/*!
 * \file nsga2.cpp
 *
 * implementation of Deb's nondominated sorting genetic algorithm, version 2
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <cmath>
#include <cfloat>
#include <climits>
#include "nsga2.h"
#include "chromosome.h"
#include "population.h"
#include "pfront.h"
#include "crossover.h"
#include "mutation.h"
#include "comparator.h"
#include "selection.h"
#include "mtrandom.h"
#include "localsearch.h"
#include "kvparse/kvparse.h"
#include "keywords.h"

using namespace std;

/*!
 * \brief constructor
 */
template <typename Encoding>
nsga2_chromosome<Encoding>::nsga2_chromosome() :
	chromosome<Encoding>::chromosome(),
	rank(0),
	crowding_distance(0)
{
}

/*!
 * \brief constructor
 */
template <typename Encoding>
nsga2_chromosome<Encoding>::nsga2_chromosome(const typename Encoding::ProblemType* prob) :
	chromosome<Encoding>::chromosome(prob),
	rank(0),
	crowding_distance(0)
{
}

/*!
 * \brief copy constructor
 */
template <typename Encoding>
nsga2_chromosome<Encoding>::nsga2_chromosome(const nsga2_chromosome& that) :
	chromosome<Encoding>::chromosome(that)
{
	rank = that.rank;
	crowding_distance = that.crowding_distance;
}

/*!
 * \brief assignment operator
 */
template <typename Encoding>
nsga2_chromosome<Encoding>& nsga2_chromosome<Encoding>::operator=(const nsga2_chromosome& that)
{
	chromosome<Encoding>::operator=(that);
	rank = that.rank;
	crowding_distance = that.crowding_distance;
	return *this;
}

/*!
 * \brief destructor
 */
template <typename Encoding>
nsga2_chromosome<Encoding>::~nsga2_chromosome()
{
}

/*!
 * \brief constructor
 */
template <typename Encoding>
crowded_comparator<Encoding>::crowded_comparator()
{
}

/*!
 * \brief destructor
 */
template <typename Encoding>
crowded_comparator<Encoding>::~crowded_comparator()
{
}

/*!
 * \brief compare chromosomes according to rank and crowding distance
 *
 * A is better than B if A has lower rank or if both have the same rank
 * and A has a lower crowding distance than B.
 */
template <typename Encoding>
int crowded_comparator<Encoding>::compare(const nsga2_chromosome<Encoding>& c1,
        const nsga2_chromosome<Encoding>& c2) const
{
	if(c1.rank < c2.rank) {
		return -1;
	} else if (c1.rank > c2.rank) {
		return 1;
	} else {
		if(c1.crowding_distance > c2.crowding_distance) {
			return -1;
		} else if (c1.crowding_distance < c2.crowding_distance) {
			return 1;
		} else {
			return 0;
		}
	}
}

/*!
 * \brief constructor
 */
template <typename Encoding>
nsga2<Encoding>::nsga2() :
	m_cross_op(0),
	m_cross_rate(1.0),
	m_mut_op(0),
	m_ls_op(0)
{
}

/*!
 * \brief destructor
 */
template <typename Encoding>
nsga2<Encoding>::~nsga2()
{
	delete m_ls_op;
	delete m_cross_op;
	delete m_mut_op;
	delete m_sel_scheme;
}

/*!
 * \brief initialize the algorithm components
 */
template <typename Encoding>
void nsga2<Encoding>::initialize()
{
	ea<nsga2_chromosome,Encoding>::initialize();

	// get the crossover operator parameters
	m_cross_op = crossover_operator_factory<nsga2_chromosome,Encoding>::construct();
	kvparse::parameter_value(keywords::CROSSOVER_RATE, m_cross_rate, true);

	// get the mutation operator parameters
	m_mut_op = mutation_operator_factory<nsga2_chromosome,Encoding>::construct();

	// initialize the selection scheme
	m_sel_scheme = new tournament_selection<nsga2_chromosome,Encoding>(&m_comp);

	// initialize the local search operator
	if(kvparse::keyword_exists("embedded_" + keywords::LOCAL_SEARCH)) {
		local_search_factory<nsga2_chromosome,Encoding> lsf;
		lsf.set_prefix("embedded_");
		m_ls_op=lsf.construct();

		strat = strategy_factory::construct();
		if(strat == STRATEGY_RANDOM) {
			kvparse::parameter_value(keywords::HC_RATE, hc_rate, true);
		}
		this->optimize_population(this->m_population);
	}
}

/*!
 * \brief find all nondominated points in the current population
 */
template <typename Encoding>
void nsga2<Encoding>::find_nondominated_front(population<nsga2_chromosome,Encoding>& pop,
        pareto_front<nsga2_chromosome,Encoding>& front)
{
	// make sure the front is empty initially
	front.clear();

	// take the first member of the population
	front.add(pop[0]);

	// for every other member of the population, check to see if the population
	// member is dominated by anyone in the front before adding it, and check to
	// see if it dominates any elements of the front and delete all such elements
	for(unsigned int i=1; i<pop.size(); i++) {
		front.add(pop[i]);
	}
}

/*!
 * \brief sort the population into non-overlapping fronts
 */
template <typename Encoding>
void nsga2<Encoding>::nondominated_sort(population<nsga2_chromosome,Encoding>& pop,
                                        vector<pareto_front<nsga2_chromosome,Encoding> >& fronts)
{
	// clear any previous sorting
	fronts.clear();

	// sort the population into separate fronts
	unsigned int fnum = 0;
	while(pop.size() > 0) {
		pareto_front<nsga2_chromosome,Encoding> curr;
		find_nondominated_front(pop, curr);
		for(unsigned int i=0; i<curr.size(); i++) {
			curr[i].rank = fnum;
		}
		fronts.push_back(curr);
		pop -= curr;
		fnum++;
	}
}

/*!
 * \brief compute the crowding distance for each chromosome in the given front
 */
template <typename Encoding>
void nsga2<Encoding>::compute_crowding_distances(pareto_front<nsga2_chromosome,Encoding>& pop)
{
	for(unsigned int i=0; i<pop.size(); i++) {
		pop[i].crowding_distance = 0;
	}

	unsigned int objectives = pop[0].fitness.size();
	for(unsigned int m=0; m<objectives; m++) {
		single_objective_comparator<nsga2_chromosome,Encoding> socomp(m);
		pop.sort(&socomp);

		// always prefer boundary points in each objective
		pop[0].crowding_distance = DBL_MAX;
		pop[pop.size()-1].crowding_distance = DBL_MAX;

		for(unsigned int i=1; i<pop.size()-1; i++) {
			pop[i].crowding_distance =
			    pop[i].crowding_distance + fabs((double)(pop[i+1].fitness[m] - pop[i-1].fitness[m]));
		}
	}
}

/*!
 * \brief run the nsga2 algorithm
 */
template <typename Encoding>
void nsga2<Encoding>::run()
{
	while(!this->terminate()) {
		run_one_generation();
		this->generation_completed(this->m_population);
	}
	pareto_front<nsga2_chromosome,Encoding> pfront(this->m_population);
	cout << pfront << endl;

	this->compute_metrics();
	this->report_metrics(cout);
}

/*!
 * \brief run one generation of the nsga2 algorithm
 */
template <typename Encoding>
void nsga2<Encoding>::run_one_generation()
{
	mtrandom mt;

	// first, we perform genetic operators on the current population
	// producing an offspring population
	m_sel_scheme->set_population(&this->m_population);
	for(unsigned int i=0; i<this->m_population.size()/2; i++) {
		// select two parents
		nsga2_chromosome<Encoding> p1 = m_sel_scheme->select_parent();
		nsga2_chromosome<Encoding> p2 = m_sel_scheme->select_parent();

		// perform crossover with some probability
		nsga2_chromosome<Encoding> c1 = p1;
		nsga2_chromosome<Encoding> c2 = p2;
		if(mt.random() < this->m_cross_rate) {
			this->m_cross_op->crossover(p1, p2, c1, c2);
		}

		// mutate the offspring
		this->m_mut_op->mutate(c1);
		this->m_mut_op->mutate(c2);

		// repair offspring if necessary
		if(this->m_repair) {
			this->m_repair->repair(c1,this->m_fitfunc);
			this->m_repair->repair(c2,this->m_fitfunc);
		}

		// evaluate the offspring
		c1.evaluate(this->m_fitfunc);
		c2.evaluate(this->m_fitfunc);
		this->chromosome_evaluated(c1);
		this->chromosome_evaluated(c2);

		m_offspring.add(c1);
		m_offspring.add(c2);
	}

	// improve the offspring population
	if(m_ls_op) {
		optimize_population(m_offspring);
	}

	// now we must generate a new parent population from the previous population
	// and the newly generated offspring population

	// sort the union of parents and offspring into domination ranks
	population<nsga2_chromosome,Encoding> Rt = this->m_population + m_offspring;
	m_offspring.clear();
	vector<pareto_front<nsga2_chromosome,Encoding> > F;
	nondominated_sort(Rt, F);

	// empty the population (and remember the size)
	unsigned int popsize = this->m_population.size();
	this->m_population.clear();

	// take as many whole nondominated fronts as possible
	unsigned int fi = 0;
	while(this->m_population.size() + F[fi].size() < popsize) {
		compute_crowding_distances(F[fi]);
		this->m_population += F[fi++];
	}

	// if the current front cannot be taken in its entirety, use the
	// crowding distance to decide which elements to take
	compute_crowding_distances(F[fi]);
	F[fi].sort(&m_comp);
	unsigned int k = 0;
	while(this->m_population.size() < popsize) {
		this->m_population.add(F[fi][k++]);
	}
}

/*!
 * \brief optimize entire population via local search method
 *
 * \todo change the way this whole thing works...this is pretty grody
 */
template <typename Encoding>
void nsga2<Encoding>::optimize_population(population<nsga2_chromosome,Encoding>& pop)
{
	scalarizing_comparator<nsga2_chromosome,Encoding> comp;
	comp.weights.resize(pop[0].fitness.size());
	if(strat == STRATEGY_ALL) {
		for(unsigned int i=0; i<pop.size(); i++) {
			comp.weights[0] = double(pop.size()-i-1)/(pop.size()-1);
			comp.weights[1] = 1.0 - comp.weights[0];
			m_ls_op->improve(pop[i], &comp, this->m_fitfunc);
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
		m_ls_op->improve(pop[bestidx], &comp, this->m_fitfunc);
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
		m_ls_op->improve(pop[worstidx], &comp, this->m_fitfunc);
	} else if(strat == STRATEGY_RANDOM) {
		mtrandom mt;
		comp.weights[0] = mt.random();
		comp.weights[1] = 1.0 - comp.weights[0];
		int idx = mt.random(0,pop.size());
		m_ls_op->improve(pop[idx], &comp, this->m_fitfunc);
	} else if(strat == STRATEGY_NONE) {
		return;
	}
}
