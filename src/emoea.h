/*!
 * \file emoea.h
 *
 * Deb et. al.'s epsilon multiobjective evolutionary algorithm
 *
 * Deon Garrett
 * deong@acm.org
 */

#ifndef _EMOEA_H_
#define _EMOEA_H_

#include "ea.h"
#include "chromosome.h"
#include "encoding.h"
#include "comparator.h"
#include "population.h"
#include "localsearch.h"

using namespace std;

/*!
 * \class emoea_chromosome
 *
 * augments standard chromosome with grid information
 *
 * \author deong
 * \date 05/09/2007
 */
template <typename Encoding>
class emoea_chromosome : public chromosome<Encoding>
{
public:
    vector<double> B;
    unsigned int index;
  
public:
    emoea_chromosome();
    emoea_chromosome(const typename Encoding::ProblemType* p);
    emoea_chromosome(const emoea_chromosome& that);
    emoea_chromosome& operator=(const emoea_chromosome& that);
    virtual ~emoea_chromosome();

    void compute_box(const vector<double>& eps);
    bool eps_dominates(const emoea_chromosome<Encoding>& other) const;
};

/*!
 * \class emoea
 *
 * Epsilon Multiobjective Evolutionary Algorithm
 *
 * \author deong
 * \date 05/09/2007
 */
template <typename Encoding>
class emoea : public ea<emoea_chromosome,Encoding>
{
private:
    // disable copying of heavyweight classes
    emoea(const emoea& that);
    emoea& operator=(const emoea& that);
    
protected:
    // archive of nondominated solutions
    pareto_front<emoea_chromosome,Encoding> m_archive;

    // genetic operators
    crossover_operator<emoea_chromosome,Encoding>* m_cross_op;
    double m_cross_rate;
    mutation_operator<emoea_chromosome,Encoding>* m_mut_op;

    // hill climbing operator
    local_search<emoea_chromosome,Encoding>* m_hc;
    next_descent<emoea_chromosome,Encoding>* m_hc2;
    
    // epsilon used for archive acceptance procedure
    vector<double> m_epsilon;
    
    // tournament selection with pareto dominance used to select
    // parent from population; random selection used to select
    // parent from archive
    tournament_selection<emoea_chromosome,Encoding> m_pop_sel;
    random_selection<emoea_chromosome,Encoding> m_arch_sel;
    
    // use the epsilon dominance comparator to update the archive
    epsilon_dominance_comparator<emoea_chromosome,Encoding> m_eps_comp;

    // use a conventional domination check for population acceptance
    pareto_dominance_comparator<emoea_chromosome,Encoding> m_pareto_comp;

    // keep track of the current "number" of each offspring (0 to popsize-1)
    // and the associated comparator 
    unsigned int index;

    // determines whether or not to use the diversification method
    bool diversification;
    unsigned int m_ls_iter;
    unsigned int ph2iter;
    vector<int> used;
    
    // helper methods
    void population_acceptance(const emoea_chromosome<Encoding>& chr);
    void archive_acceptance(const emoea_chromosome<Encoding>& chr);
    void intensify(emoea_chromosome<Encoding>& chr);
    void diversify(emoea_chromosome<Encoding>& chr);
    
public:
    emoea();
    virtual ~emoea();

    virtual void initialize();
    virtual void run();
    virtual void run_one_generation();
};

#include "emoea.cpp"

#endif
