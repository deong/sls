/*!
 * \file nsga2.h
 *
 * implementation of Deb's nondominated sorting genetic algorithm, version 2
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _NSGA2_H_
#define _NSGA2_H_

#include "chromosome.h"
#include "comparator.h"
#include "ea.h"
#include "population.h"
#include "pfront.h"
#include "crossover.h"
#include "mutation.h"
#include "localsearch.h"

using namespace std;

/*!
 * \class nsga2_chromosome
 *
 * adds crowding distance and rank to standard chromosomes
 */
template <typename Encoding>
class nsga2_chromosome : public chromosome<Encoding>
{
public:
	// nondomination rank
	unsigned int rank;

	// crowding distance
	double crowding_distance;

public:
	nsga2_chromosome();
	nsga2_chromosome(const typename Encoding::ProblemType* prob);
	nsga2_chromosome(const nsga2_chromosome& that);
	nsga2_chromosome& operator=(const nsga2_chromosome& that);
	virtual ~nsga2_chromosome();
};

/*!
 * \class crowded_comparator
 *
 * used to compare two nsga2 chromosomes
 */
template <typename Encoding>
class crowded_comparator : public comparator<nsga2_chromosome,Encoding>
{
public:
	crowded_comparator();
	~crowded_comparator();

	virtual inline int compare(const nsga2_chromosome<Encoding>& c1,
	                           const nsga2_chromosome<Encoding>& c2) const;
};

/*!
 * \class nsga2
 *
 * Deb et. al.'s Nondominated Sorting Genetic Algorithm
 */
template <typename Encoding>
class nsga2 : public ea<nsga2_chromosome,Encoding>
{
private:
	// disable copying of heavyweight class
	nsga2(const nsga2& that);
	nsga2& operator=(const nsga2& that);

	// simplifies implementation to keep offspring population
	// as a class member so it is shared between methods
	population<nsga2_chromosome,Encoding> m_offspring;

protected:
	// compare chromosomes using crowding distance comparator
	crowded_comparator<Encoding> m_comp;

	// tournament selection (using crowding distance comparator
	tournament_selection<nsga2_chromosome,Encoding>* m_sel_scheme;

	// allow whatever crossover operator we see fit
	crossover_operator<nsga2_chromosome,Encoding>* m_cross_op;
	double m_cross_rate;

	// allow whatever mutation operator we see fit
	mutation_operator<nsga2_chromosome,Encoding>* m_mut_op;

	// local search operator to improve individuals
	local_search<nsga2_chromosome,Encoding>* m_ls_op;
	unsigned int m_ls_iter;

	// compute the crowding distances
	void compute_crowding_distances(pareto_front<nsga2_chromosome,Encoding>& pop);

	// find the nondominated front
	void find_nondominated_front(population<nsga2_chromosome,Encoding>& pop,
	                             pareto_front<nsga2_chromosome,Encoding>& front);

	// sort the population into nondomination ranks
	void nondominated_sort(population<nsga2_chromosome,Encoding>& pop,
	                       vector<pareto_front<nsga2_chromosome,Encoding> >& fronts);

public:
	// ctors and dtor
	nsga2();
	virtual ~nsga2();

	// initialize the nsga2 algorithm
	virtual void initialize();

	// run the nsga2 algorithm
	virtual void run();

	// run a single generation of the nsga2 algorithm
	virtual void run_one_generation();
};

#include "nsga2.cpp"

#endif

