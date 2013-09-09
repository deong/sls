/*!
 * \file spea2.h
 *
 * Zitzler et.al.'s Strength Pareto Evolutionary Algorithm 2
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _SPEA2_H_
#define _SPEA2_H_

#include <vector>
#include "chromosome.h"
#include "comparator.h"
#include "selection.h"
#include "crossover.h"
#include "mutation.h"
#include "repair.h"
#include "problems.h"
#include "population.h"
#include "localsearch.h"
#include "strategy.h"

using namespace std;

typedef struct
{
    int source;
    int dest;
    double distance;
} neighbor;

/*!
 * \class neighbor_comparator
 *
 * order neighbors by increasing order of distance
 */
class neighbor_comparator
{
public:
    bool operator()(const neighbor& n1, const neighbor& n2)
    {
        return n1.distance < n2.distance;
    }
};

/*!
 * \class spea2_chromosome
 */
template <typename Encoding>
class spea2_chromosome : public chromosome<Encoding>
{
public:
    double strength;
    double spea2_fitness;
    vector<neighbor> neighbors;

public:
    spea2_chromosome();
    spea2_chromosome(const typename Encoding::ProblemType* prob);
    spea2_chromosome(const spea2_chromosome& that);
    spea2_chromosome& operator=(const spea2_chromosome& that);
    virtual ~spea2_chromosome();

    virtual double distance_to(const spea2_chromosome& other) const;
    virtual void init_neighbors(unsigned int size);
    virtual void sort_neighbors();
    virtual void remove_neighbor(int index);
};

/*!
 * \class spea2_comparator
 */
template <typename Encoding>
class spea2_comparator : public comparator<spea2_chromosome,Encoding>
{
public:
    spea2_comparator();
    virtual ~spea2_comparator();
    
    virtual inline int compare(const spea2_chromosome<Encoding>& c1,
                               const spea2_chromosome<Encoding>& c2) const;
};

/*!
 * \class spea2
 *
 * The Strength Pareto Evolutionary Algorithm
 */
template <typename Encoding>
class spea2 : public ea<spea2_chromosome,Encoding>
{
private:
    // disable copying
    spea2(const spea2& that);
    spea2& operator=(const spea2& that);

protected:
    // store nondominated solutions an archive
    population<spea2_chromosome,Encoding> m_archive;
    population<spea2_chromosome,Encoding> m_all;
    unsigned int m_archive_size;

    // which neighbor should be used for density estimation
    unsigned int m_spea2_k;
    
    // tournament selection using spea2_fitness
    tournament_selection<spea2_chromosome,Encoding>* m_sel_scheme;
    spea2_comparator<Encoding> m_comp;
    
    // allow any crossover
    crossover_operator<spea2_chromosome,Encoding>* m_cross_op;
    double m_cross_rate;

    // allow any mutation operator
    mutation_operator<spea2_chromosome,Encoding>* m_mut_op;

    // local search operator
    strategy strat;
    double hc_rate;
    local_search<spea2_chromosome,Encoding>* m_hc;

    // compute the strength of each nondominated member
    void compute_strength();

    // compute the fitness of dominated individuals
    void compute_fitness();

    // compute the density estimates for each individual
    void compute_density();

    // construct the archive
    void construct_archive();
    void truncate_archive();

    // subject some individuals to local search
    void optimize_population(population<spea2_chromosome,Encoding>& pop);
    
public:
    spea2();
    virtual ~spea2();

    // initialize the algorithm
    virtual void initialize();

    // run the spea2 algorithm
    virtual void run();
    virtual void run_one_generation();
};

#include "spea2.cpp"

#endif


    
