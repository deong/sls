/*!
 * \file simplega.h
 *
 * Simple Genetic Algorithm (generational ga)
 *
 * Deon Garrett
 * deong@acm.org
 */

#ifndef _SIMPLEGA_H_
#define _SIMPLEGA_H_

#include "ea.h"
#include "selection.h"
#include "replacement.h"
#include "crossover.h"
#include "mutation.h"

/*!
 * \class simple_ga
 *
 * \author deong
 * \date 05/11/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class simple_ga : public ea<Chromosome,Encoding>
{
private:
    // disable copying of heavyweight class
    simple_ga(const simple_ga& that);
    simple_ga& operator=(const simple_ga& that);

protected:
    // utilize the ranking selection scheme
    selection_scheme<Chromosome,Encoding>* m_sel_scheme;

    // utilize a replace-worst environment selection scheme
    replacement_scheme<Chromosome,Encoding>* m_rep_scheme;

    // allow whatever crossover operator we see fit
    crossover_operator<Chromosome,Encoding>* m_cross_op;
    double m_cross_rate;
    
    // allow whatever mutation operator we see fit
    mutation_operator<Chromosome,Encoding>* m_mut_op;

    // define a comparator operator
    comparator<Chromosome,Encoding>* m_comp;
    
public:
    // ctors and dtor
    simple_ga();
    virtual ~simple_ga();

    // initialize the simple_ga algorithm
    virtual void initialize();

    // run the simple_ga algorithm
    virtual void run();

    // run a single generation of the simple_ga algorithm
    virtual void run_one_generation();
};

#include "simplega.cpp"

#endif 
