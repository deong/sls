/*!
 * \file emoea.cpp
 *
 * Deb et. al.'s epsilon multiobjective evolutionary algorithm
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <climits>
#include "emoea.h"
#include "ea.h"
#include "chromosome.h"
#include "population.h"
#include "comparator.h"
#include "kvparse/kvparse.h"
#include "keywords.h"
#include "utilities.h"
#include "localsearch.h"

using namespace std;

/*!
 * \brief constructor
 */
template <typename Encoding>
emoea_chromosome<Encoding>::emoea_chromosome() :
	chromosome<Encoding>::chromosome()
{
}

/*!
 * \brief constructor
 */
template <typename Encoding>
emoea_chromosome<Encoding>::emoea_chromosome(const typename Encoding::ProblemType* prob) :
	chromosome<Encoding>::chromosome(prob)
{
	B.resize(prob->objectives());
}

/*!
 * \brief copy constructor
 */
template <typename Encoding>
emoea_chromosome<Encoding>::emoea_chromosome(const emoea_chromosome& that) :
	chromosome<Encoding>::chromosome(that)
{
	B = that.B;
	index = that.index;
}

/*!
 * \brief assignment operator
 */
template <typename Encoding>
emoea_chromosome<Encoding>& emoea_chromosome<Encoding>::operator=(const emoea_chromosome& that)
{
	chromosome<Encoding>::operator=(that);
	B = that.B;
	index = that.index;
	return *this;
}

/*!
 * \brief destructor
 */
template <typename Encoding>
emoea_chromosome<Encoding>::~emoea_chromosome()
{
}

/*!
 * \brief compute the box containing the chromosome
 */
template <typename Encoding>
void emoea_chromosome<Encoding>::compute_box(const vector<double>& eps)
{
	for(unsigned int i=0; i<B.size(); i++) {
		B[i] = floor(this->fitness[i]/eps[i]);
	}
}

/*!
 * \brief determine whether a chromosome epsilon-dominates another
 *
 * \bug why isn't this using the epsilon_dominance_comparator?
 * \sa epsilon_dominance_comparator
 */
template <typename Encoding>
bool emoea_chromosome<Encoding>::eps_dominates(const emoea_chromosome<Encoding>& other) const
{
	bool better = false;
	for(unsigned int i=0; i<B.size(); i++) {
		if(B[i] < other.B[i]) {
			better = true;
		} else if(B[i] > other.B[i]) {
			return false;
		}
	}
	return better;
}

/*!
 * \brief constructor
 */
template <typename Encoding>
emoea<Encoding>::emoea() :
	ea<emoea_chromosome,Encoding>::ea(),
	m_cross_op(0),
	m_cross_rate(1.0),
	m_mut_op(0),
	m_hc(0)
{
}

/*!
 * \brief destructor
 */
template <typename Encoding>
emoea<Encoding>::~emoea()
{
	if(m_hc) {
		delete m_hc;
	}
	delete m_cross_op;
	delete m_mut_op;
}

/*!
 * \brief initialize the emoea components
 */
template <typename Encoding>
void emoea<Encoding>::initialize()
{
	// initialize the comparator objects
	ea<emoea_chromosome,Encoding>::initialize();

	// initialize the epsilon dominance comparator
	m_eps_comp.initialize();

	// read in the same epsilon for use in the archive acceptance procedure
	kvparse::parameter_value(keywords::EPSILON, m_epsilon, true);

	// set up the crossover and mutation operators
	m_cross_op = crossover_operator_factory<emoea_chromosome,Encoding>::construct();
	kvparse::parameter_value(keywords::CROSSOVER_RATE, m_cross_rate, true);
	m_mut_op = mutation_operator_factory<emoea_chromosome,Encoding>::construct();

	// configure the hill climber
	if(kvparse::keyword_exists(keywords::LOCAL_SEARCH)) {
		local_search_factory<emoea_chromosome,Encoding> lsf;
		lsf.set_prefix("ls_");
		m_hc=lsf.construct();

		diversification = true;
		kvparse::parameter_value(keywords::DIVERSIFICATION, diversification, false);

		if(diversification) {
			// figure out how many iterations to run for the second phase
			// of local search
			m_ls_iter = INT_MAX;
			kvparse::parameter_value(keywords::LOCAL_SEARCH_ITERATIONS, m_ls_iter, true);
			ph2iter = static_cast<unsigned int>(sqrt(double(m_ls_iter)));
			kvparse::parameter_value(keywords::SECOND_PHASE_ITERATIONS, ph2iter, false);
		}
	}

	// initialize the tournament comparator
	m_pop_sel.initialize();

	// compute the boxes for the initial population
	for(unsigned int i=0; i<this->m_population.size(); i++) {
		this->m_population[i].compute_box(m_epsilon);
	}
}

/*!
 * \brief check a chromosome for acceptance into a population
 */
template <typename Encoding>
void emoea<Encoding>::population_acceptance(const emoea_chromosome<Encoding>& chr)
{
	// check to see if the solution dominates any population members
	// or is dominated by any
	vector<int> ties;
	for(unsigned int i=0; i<this->m_population.size(); i++) {
		int res = m_pareto_comp.compare(chr,this->m_population[i]);
		if(res == -1) {
			ties.push_back(i);
		} else if(res == 1) {
			// if a population member dominates the solution,
			// the solution is not added
			return;
		}
	}

	if(ties.size() > 0) {
		mtrandom mt;
		int to_replace = mt.random(0, ties.size());
		this->m_population[ties[to_replace]] = chr;
	} else {
		// if the population and the solution are mutually nondominated,
		// then the solution replaces a randomly selected population member
		mtrandom mt;
		this->m_population[mt.random(this->m_population.size())] = chr;
	}
}

/*!
 * \brief check a chromosome for acceptance into the archive
 */
template <typename Encoding>
void emoea<Encoding>::archive_acceptance(const emoea_chromosome<Encoding>& chr)
{
	int same_box = -1;
	for(unsigned int i=0; i<m_archive.size(); i++) {
		// if any archive member eps_dominates chr, then chr is rejected
		if(m_archive[i].eps_dominates(chr)) {
			return;
		}
		// if chr eps_dominates any archive member, it replaces the member
		else if(chr.eps_dominates(m_archive[i])) {
			m_archive[i] = chr;
			for(unsigned int j=m_archive.size()-1; j>i; j--) {
				if(chr.eps_dominates(m_archive[j])) {
					m_archive.remove_at(j);
				}
			}
			return;
		}
		// remember if the child shares a box with an archive member
		else if(chr.B == m_archive[i].B) {
			same_box = i;
			break;
		}
	}

	// if we get here, then the child and archive are mutually non-eps-dominated
	if(same_box != -1) {
		// if the child shares a box with an archive member, it replaces the archive
		// member if it dominates the member or if it is closer to their B vector
		if(m_pareto_comp.compare(chr,m_archive[same_box]) == -1) {
			m_archive[same_box] = chr;
		} else {
			double d1 = 0;
			double d2 = 0;
			for(unsigned int i=0; i<chr.fitness.size(); i++) {
				d1 += (chr.fitness[i]-chr.B[i]) * (chr.fitness[i]-chr.B[i]);
				d2 += (m_archive[same_box].fitness[i] - m_archive[same_box].B[i]) *
				      (m_archive[same_box].fitness[i] - m_archive[same_box].B[i]);
			}
			d1 = sqrt(d1);
			d2 = sqrt(d2);
			if(d1 < d2) {
				m_archive[same_box] = chr;
			}
		}
	} else {
		// if we are here, the archive doesn't dominate the child and the
		// child lies in a new box, so we add it to the archive
		m_archive.add(chr);
	}
}

/*!
 * \brief intensify the search using a targeted local search run
 */
template <typename Encoding>
void emoea<Encoding>::intensify(emoea_chromosome<Encoding>& chr)
{
	scalarizing_comparator<emoea_chromosome,Encoding> comp;
	comp.weights.resize(chr.fitness.size());
	comp.weights[0] = double(this->m_population.size()-chr.index-1)/(this->m_population.size()-1);
	comp.weights[1] = 1.0 - comp.weights[0];

	m_hc->improve(chr, &comp, this->m_fitfunc);
	chr.compute_box(m_epsilon);

	// put the new point into the population and archive (if necessary)
	population_acceptance(chr);
	archive_acceptance(chr);

	used[chr.index] = 1;
}

/*!
 * \brief diversify the search using "randomly" seeded local search runs
 */
template <typename Encoding>
void emoea<Encoding>::diversify(emoea_chromosome<Encoding>& chr)
{
	scalarizing_comparator<emoea_chromosome,Encoding> comp;
	comp.weights.resize(chr.fitness.size());

	for(int offset=-1; offset<=1; offset+=2) {
		emoea_chromosome<Encoding> curr(chr);

		// set the new index for the solution
		curr.index = min(max(0, int(chr.index+offset)),int(this->m_population.size()-1));
		if(used[curr.index] == 0) {
			comp.weights[0] = double(this->m_population.size()-curr.index-1)/(this->m_population.size()-1);
			comp.weights[1] = 1.0 - comp.weights[0];

			// run shorter ts runs
			m_hc->improve(curr, &comp, this->m_fitfunc);
			curr.compute_box(m_epsilon);

			// put the new point into the population and archive (if necessary)
			population_acceptance(curr);
			archive_acceptance(curr);

			used[curr.index] = 1;
		}
	}
}

/*!
 * \brief run the epsilon MOEA
 */
template <typename Encoding>
void emoea<Encoding>::run()
{
	mtrandom mt;
	vector<int> order = mt.permutation(this->m_population.size());

	used.resize(this->m_population.size());
	used.assign(this->m_population.size(), 0);

	if(m_hc) {
		// optimize initial population
		for(index = 0; index<this->m_population.size(); index++) {
			this->m_population[index].index = index;
			intensify(this->m_population[index]);
		}
	}
	m_archive.construct_front(this->m_population);
	this->generation_completed(m_archive);

	unsigned int count = 0;
	while(!this->terminate()) {
		index = order[count];
		if(used[index] == 0) {
			run_one_generation();
		}

		if(++count % this->m_population.size() == 0) {
			count = 0;
			order = mt.permutation(this->m_population.size());
			used.assign(this->m_population.size(), 0);
			this->generation_completed(m_archive);
		}
	}
	cout << m_archive << endl;

	this->compute_metrics();
	this->report_metrics(cout);
}

/*!
 * \brief run one generation of the epsilon MOEA
 */
template <typename Encoding>
void emoea<Encoding>::run_one_generation()
{
	mtrandom mt;

	m_pop_sel.set_population(&this->m_population);
	m_arch_sel.set_population(&m_archive);

	// select a parent from the population
	emoea_chromosome<Encoding> p1 = m_arch_sel.select_parent();

	// and one from the archive
	emoea_chromosome<Encoding> p2 = m_pop_sel.select_parent();

	// perform crossover, mutation
	// (discard c2 as only a single offspring required)
	emoea_chromosome<Encoding> c1 = p1;
	emoea_chromosome<Encoding> c2 = p2;

	if(mt.random() < m_cross_rate) {
		m_cross_op->crossover(p1, p2, c1, c2);
	}

	// randomly choose between c1 and c2
	if(mt.random() < 0.5) {
		c1 = c2;
	}

	// mutate the offspring
	m_mut_op->mutate(c1);

	// evaluate the offspring
	if(this->m_repair) {
		this->m_repair->repair(c1,this->m_fitfunc);
	}
	c1.evaluate(this->m_fitfunc);
	c1.compute_box(m_epsilon);
	this->chromosome_evaluated(c1);

	// improve the offspring
	if(m_hc) {
		// intensify the search toward the pareto front using a
		// particular weight vector
		c1.index = index;
		if(used[c1.index] == 0) {
			intensify(c1);

			if(diversification) {
				diversify(c1);
			}
		}
	} else {
		// include offspring in population
		population_acceptance(c1);

		// include offspring in archive
		archive_acceptance(c1);
	}
}

