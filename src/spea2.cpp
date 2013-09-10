/*!
 * \file spea2.cpp
 *
 * Zitzler et.al.'s Strength Pareto Evolutionary Algorithm 2
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <cfloat>
#include <fstream>
#include <iterator>
#include <iomanip>
#include "spea2.h"
#include "chromosome.h"
#include "population.h"
#include "selection.h"
#include "crossover.h"
#include "mutation.h"
#include "pfront.h"
#include "strategy.h"
#include "localsearch.h"
#include "kvparse/kvparse.h"
#include "keywords.h"

using namespace std;

/*!
 * \brief constructor
 */
template <typename Encoding>
spea2_chromosome<Encoding>::spea2_chromosome() :
	chromosome<Encoding>::chromosome(),
	strength(0),
	spea2_fitness(0)
{
}

/*!
 * \brief constructor
 */
template <typename Encoding>
spea2_chromosome<Encoding>::spea2_chromosome(const typename Encoding::ProblemType* prob) :
	chromosome<Encoding>::chromosome(prob),
	strength(0),
	spea2_fitness(0)
{
}

/*!
 * \brief copy constructor
 */
template <typename Encoding>
spea2_chromosome<Encoding>::spea2_chromosome(const spea2_chromosome& that) :
	chromosome<Encoding>::chromosome(that)
{
	strength = that.strength;
	spea2_fitness = that.spea2_fitness;
	neighbors = that.neighbors;
}

/*!
 * \brief assignment operator
 */
template <typename Encoding>
spea2_chromosome<Encoding>& spea2_chromosome<Encoding>::operator=(const spea2_chromosome& that)
{
	chromosome<Encoding>::operator=(that);
	strength = that.strength;
	spea2_fitness = that.spea2_fitness;
	neighbors = that.neighbors;
	return *this;
}

/*!
 * \brief destructor
 */
template <typename Encoding>
spea2_chromosome<Encoding>::~spea2_chromosome()
{
}

/*!
 * \brief compute distance to another chromosome
 */
template <typename Encoding>
double spea2_chromosome<Encoding>::distance_to(const spea2_chromosome& other) const
{
	double d = 0.0;
	for(unsigned int i=0; i<this->fitness.size(); i++) {
		d += (this->fitness[i] - other.fitness[i]) * (this->fitness[i] - other.fitness[i]);
	}
	return sqrt(d);
}

/*!
 * \brief initialize the neighbor list
 */
template <typename Encoding>
void spea2_chromosome<Encoding>::init_neighbors(unsigned int size)
{
	neighbors.clear();
	neighbors.resize(size);
}

/*!
 * \brief sort the neighbor list by distance
 */
template <typename Encoding>
void spea2_chromosome<Encoding>::sort_neighbors()
{
	sort(neighbors.begin(),neighbors.end(),neighbor_comparator());
}

/*!
 * \brief remove a neighbor from the list
 */
template <typename Encoding>
void spea2_chromosome<Encoding>::remove_neighbor(int dest)
{
	for(unsigned int i=0; i<neighbors.size(); i++) {
		if(neighbors[i].dest == dest) {
			neighbors.erase(neighbors.begin()+i);
			return;
		}
	}
}

/*!
 * \brief constructor
 */
template <typename Encoding>
spea2_comparator<Encoding>::spea2_comparator()
{
}

/*!
 * \brief destructor
 */
template <typename Encoding>
spea2_comparator<Encoding>::~spea2_comparator()
{
}

/*!
 * \brief compare chromosomes by spea2 fitness values
 */
template <typename Encoding>
inline int spea2_comparator<Encoding>::compare(const spea2_chromosome<Encoding>& c1,
        const spea2_chromosome<Encoding>& c2) const
{
	if(c1.spea2_fitness < c2.spea2_fitness) {
		return -1;
	} else if(c1.spea2_fitness > c2.spea2_fitness) {
		return 1;
	} else {
		return 0;
	}
}

/*!
 * \brief constructor
 */
template <typename Encoding>
spea2<Encoding>::spea2() :
	ea<spea2_chromosome,Encoding>::ea(),
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
spea2<Encoding>::~spea2()
{
	if(m_hc) {
		delete m_hc;
	}
	delete m_sel_scheme;
	delete m_cross_op;
	delete m_mut_op;
}

/*!
 * \brief compute the strength of all population and archive members
 */
template <typename Encoding>
void spea2<Encoding>::compute_strength()
{
	pareto_dominance_comparator<spea2_chromosome,Encoding> pdc;

	int dominates;
	for(unsigned int i=0; i<m_archive.size(); i++) {
		dominates = 0;
		for(unsigned int j=0; j<this->m_population.size(); j++) {
			if(pdc.compare(m_archive[i],this->m_population[j]) == -1) {
				dominates++;
			}
		}

		for(unsigned int j=0; j<m_archive.size(); j++) {
			if(pdc.compare(m_archive[i],m_archive[j]) == -1) {
				dominates++;
			}
		}
		m_archive[i].strength = dominates;
	}

	for(unsigned int i=0; i<this->m_population.size(); i++) {
		dominates = 0;
		for(unsigned int j=0; j<this->m_population.size(); j++) {
			if(pdc.compare(this->m_population[i],this->m_population[j]) == -1) {
				dominates++;
			}
		}

		for(unsigned int j=0; j<m_archive.size(); j++) {
			if(pdc.compare(this->m_population[i],m_archive[j]) == -1) {
				dominates++;
			}
		}
		this->m_population[i].strength = dominates;
	}
}

/*!
 * \brief compute the fitness of all population and archive members
 */
template <typename Encoding>
void spea2<Encoding>::compute_fitness()
{
	pareto_dominance_comparator<spea2_chromosome,Encoding> pdc;

	for(unsigned int i=0; i<m_archive.size(); i++) {
		m_archive[i].spea2_fitness = 0;
		for(unsigned int j=0; j<this->m_population.size(); j++) {
			if(pdc.compare(this->m_population[j],m_archive[i]) == -1) {
				m_archive[i].spea2_fitness += this->m_population[j].strength;
			}
		}

		for(unsigned int j=0; j<m_archive.size(); j++) {
			if(pdc.compare(m_archive[j],m_archive[i]) == -1) {
				m_archive[i].spea2_fitness += m_archive[j].strength;
			}
		}
	}

	for(unsigned int i=0; i<this->m_population.size(); i++) {
		this->m_population[i].spea2_fitness = 0;
		for(unsigned int j=0; j<this->m_population.size(); j++) {
			if(pdc.compare(this->m_population[j],this->m_population[i]) == -1) {
				this->m_population[i].spea2_fitness += this->m_population[j].strength;
			}
		}

		for(unsigned int j=0; j<m_archive.size(); j++) {
			if(pdc.compare(m_archive[j],this->m_population[i]) == -1) {
				this->m_population[i].spea2_fitness += m_archive[j].strength;
			}
		}
	}
}

/*!
 * \brief compute the density measure for each chromosome
 */
template <typename Encoding>
void spea2<Encoding>::compute_density()
{
	unsigned int n = m_archive.size() + this->m_population.size();

	// build up m_all
	m_all.clear();
	m_all.resize(n);
	for(unsigned int i=0; i<this->m_population.size(); i++) {
		this->m_population[i].init_neighbors(n-1);
		m_all.add(this->m_population[i]);
	}
	for(unsigned int i=0; i<m_archive.size(); i++) {
		m_archive[i].init_neighbors(n-1);
		m_all.add(m_archive[i]);
	}

	for(unsigned int i=0; i<n-1; i++) {
		for(unsigned int j=i+1; j<n; j++) {
			double d = m_all[i].distance_to(m_all[j]);

			neighbor n1;
			n1.source = i;
			n1.dest = j;
			n1.distance = d;
			m_all[i].neighbors[j-1] = n1;

			neighbor n2;
			n2.source = j;
			n2.dest = i;
			n2.distance = d;
			m_all[j].neighbors[i] = n2;
		}
	}

	for(unsigned int i=0; i<n; i++) {
		m_all[i].sort_neighbors();
		double sigma_k = m_all[i].neighbors[m_spea2_k].distance;
		double di = 1.0 / (sigma_k + 2.0);
		m_all[i].spea2_fitness += di;
	}
}

/*!
 * \brief construct the next archive from the population and current archive
 */
template <typename Encoding>
void spea2<Encoding>::construct_archive()
{
	compute_strength();
	compute_fitness();
	compute_density();

	// sort the union of archive and population
	m_all.sort(&m_comp);

	// add all the nondominated solutions to the new archive
	m_archive.clear();
	m_archive.resize(m_archive_size);
	unsigned int k = 0;
	while(k < m_all.size() && m_all[k].spea2_fitness < 1) {
		m_archive.add(m_all[k++]);
	}

	// count the number of items over/under the required number
	int numshort = m_archive_size - m_archive.size();

	// if the archive is too small, add the best dominated individuals
	if(numshort > 0) {
		while(k < min(m_all.size(),m_archive_size)) {
			m_archive.add(m_all[k++]);
		}
	}
	// otherwise, start pruning the archive by density
	else {
		// remove the neighbors corresponding to points not taken in
		// the archive
		for(unsigned int i=0; i<m_archive.size(); i++) {
			int tmp = k;
			for(int j=int(m_archive[i].neighbors.size()-1); j>=0; j--) {
				if(m_archive[i].neighbors[j].dest == m_all[tmp].neighbors[0].source) {
					m_archive[i].neighbors.erase(m_archive[i].neighbors.begin()+j);
					break;
				}
			}
			tmp++;
		}

		truncate_archive();
	}
}

/*!
 * \brief reduce archive size to acceptable levels
 */
template <typename Encoding>
void spea2<Encoding>::truncate_archive()
{
	int excess = m_archive.size() - m_archive_size;
	while(excess > 0) {
		// find the individual whose nearst neighbor is closest
		double min_dist = DBL_MAX;
		vector<int> candidates;
		for(unsigned int i=0; i<m_archive.size(); i++) {
			double d = m_archive[i].neighbors[0].distance;
			if(d < min_dist) {
				min_dist = d;
				candidates.clear();
				candidates.push_back(i);
			} else if(fabs(d-min_dist) < 1e-9) {
				candidates.push_back(i);
			}
		}

		// now candidates contains the indices of all individuals
		// with equidistant nearest neighbors
		int nTies = candidates.size();
		unsigned int neighbor = 1;
		while(nTies > 1 && neighbor < m_archive[0].neighbors.size()) {
			double min_dist = DBL_MAX;

			// keep a vector of points to maintain for another iteration
			vector<int> still_tied;

			// find all the points that are still tied after considering the next neighbor
			for(int i=candidates.size()-1; i>=0; i--) {
				// d is the distance to the current neighbor for the last candidate point
				double d = m_archive[candidates[i]].neighbors[neighbor].distance;

				if(fabs(d-min_dist) < 1e-9) {
					still_tied.push_back(candidates[i]);
				} else if(d < min_dist) {
					min_dist = d;
					still_tied.clear();
					still_tied.push_back(candidates[i]);
				}
			}

			nTies = still_tied.size();
			candidates = still_tied;

			// go on to the next neighbor
			neighbor++;
		}

		// at this point, we can remove the remaining candidate
		int removed_neighbor = m_archive[candidates[0]].neighbors[0].source;
		m_archive.remove_at(candidates[0]);
		for(unsigned int i=0; i<m_archive.size(); i++) {
			m_archive[i].remove_neighbor(removed_neighbor);
		}

		// we now have one fewer excess point in the archive
		excess--;
	}
}

/*!
 * \brief initialize the algorithm components
 */
template <typename Encoding>
void spea2<Encoding>::initialize()
{
	// call the base class initializer
	ea<spea2_chromosome,Encoding>::initialize();

	// set up the selection scheme
	m_sel_scheme = new tournament_selection<spea2_chromosome,Encoding>(&m_comp);

	// get the crossover operator parameters
	m_cross_op = crossover_operator_factory<spea2_chromosome,Encoding>::construct();
	kvparse::parameter_value(keywords::CROSSOVER_RATE, m_cross_rate, true);

	// get the mutation operator parameters
	m_mut_op = mutation_operator_factory<spea2_chromosome,Encoding>::construct();

	// configure the hill climber
	if(kvparse::keyword_exists("embedded_" + keywords::LOCAL_SEARCH)) {
		local_search_factory<spea2_chromosome,Encoding> lsf;
		lsf.set_prefix("embedded_");
		m_hc = lsf.construct();

		strat = strategy_factory::construct();
		if(strat == STRATEGY_RANDOM) {
			kvparse::parameter_value(keywords::HC_RATE, hc_rate, true);
		}
		this->optimize_population(this->m_population);
	}

	// get the archive size
	m_archive_size = this->m_population.size();
	kvparse::parameter_value(keywords::ARCHIVE_SIZE, m_archive_size, false);

	// get the density estimation parameter
	m_spea2_k = static_cast<int>(sqrt(double(this->m_population.size() + m_archive_size)));
	kvparse::parameter_value(keywords::SPEA2_K, m_spea2_k, false);
}

/*!
 * \brief run the spea2 algorithm
 */
template <typename Encoding>
void spea2<Encoding>::run()
{
	// build up the initial archive
	construct_archive();
	pareto_front<spea2_chromosome,Encoding> tmp(m_archive);

	while(!this->terminate()) {
		run_one_generation();
		pareto_front<spea2_chromosome,Encoding> front(m_archive);
		this->generation_completed(front);
	}

	pareto_front<spea2_chromosome,Encoding> front(m_archive);
	cout << front << endl;

	this->compute_metrics();
	this->report_metrics(cout);
}

/*!
 * \brief run one generation of the spea2 algorithm
 */
template <typename Encoding>
void spea2<Encoding>::run_one_generation()
{
	mtrandom mt;

	population<spea2_chromosome,Encoding> offspring(this->m_population.size());
	m_sel_scheme->set_population(&m_archive);
	for(unsigned int i=0; i<this->m_population.size()/2; i++) {
		// select two parents
		spea2_chromosome<Encoding> p1 = m_sel_scheme->select_parent();
		spea2_chromosome<Encoding> p2 = m_sel_scheme->select_parent();

		// perform crossover with some probability
		spea2_chromosome<Encoding> c1 = p1;
		spea2_chromosome<Encoding> c2 = p2;
		if(mt.random() < this->m_cross_rate) {
			this->m_cross_op->crossover(p1, p2, c1, c2);
		}

		// mutate the offspring
		this->m_mut_op->mutate(c1);
		this->m_mut_op->mutate(c2);

		// evaluate the offspring
		if(this->m_repair) {
			this->m_repair->repair(c1,this->m_fitfunc);
		}
		c1.evaluate(this->m_fitfunc);
		this->chromosome_evaluated(c1);

		if(this->m_repair) {
			this->m_repair->repair(c2,this->m_fitfunc);
		}
		c2.evaluate(this->m_fitfunc);
		this->chromosome_evaluated(c2);

		offspring.add(c1);
		offspring.add(c2);
	}

	// improve the offspring population
	if(m_hc) {
		optimize_population(offspring);
	}

	this->m_population = offspring;
	construct_archive();
}

/*!
 * \brief optimize entire population via local search method
 *
 * \todo change the way this whole thing works...this is pretty grody
 */
template <typename Encoding>
void spea2<Encoding>::optimize_population(population<spea2_chromosome,Encoding>& pop)
{
	scalarizing_comparator<spea2_chromosome,Encoding> comp;
	comp.weights.resize(pop[0].fitness.size());
	if(strat == STRATEGY_ALL) {
		for(unsigned int i=0; i<pop.size(); i++) {
			comp.weights[0] = double(pop.size()-i-1)/(pop.size()-1);
			comp.weights[1] = 1.0 - comp.weights[0];
			m_hc->improve(pop[i], &comp, this->m_fitfunc);
		}
	} else if(strat == STRATEGY_BEST) {
		mtrandom mt;
		comp.weights[0] = mt.random();
		comp.weights[1] = 1.0 - comp.weights[0];
		int bestidx = 0;
		for(unsigned int i=1; i<pop.size(); i++) {
			if(comp.compare(pop[i],pop[bestidx]) == -1) {
				bestidx = i;
			}
		}
		m_hc->improve(pop[bestidx], &comp, this->m_fitfunc);
	} else if(strat == STRATEGY_WORST) {
		mtrandom mt;
		comp.weights[0] = mt.random();
		comp.weights[1] = 1.0 - comp.weights[0];
		int worstidx = 0;
		for(unsigned int i=1; i<pop.size(); i++) {
			if(comp.compare(pop[i],pop[worstidx]) == 1) {
				worstidx = i;
			}
		}
		m_hc->improve(pop[worstidx], &comp, this->m_fitfunc);
	} else if(strat == STRATEGY_RANDOM) {
		mtrandom mt;
		comp.weights[0] = mt.random();
		comp.weights[1] = 1.0 - comp.weights[0];
		int idx = mt.random(0,pop.size());
		m_hc->improve(pop[idx], &comp, this->m_fitfunc);
	} else if(strat == STRATEGY_NONE) {
		return;
	}
}
