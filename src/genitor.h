/*!
 * \file genitor.h
 *
 * Whitley's steady state ga which first advocated the use of
 * ranking for parent selection
 *
 * Deon Garrett
 * deong@acm.org
 */

#ifndef _GENITOR_H_
#define _GENITOR_H_

#include "ea.h"
#include "selection.h"
#include "replacement.h"
#include "crossover.h"
#include "mutation.h"

/*!
 * \class genitor
 *
 * Darrell Whitley's Genitor algorithm
 *
 * \author deong
 * \date 05/09/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class genitor : public ea<Chromosome,Encoding>
{
private:
    // disable copying of heavyweight class
    genitor(const genitor& that);
    genitor& operator=(const genitor& that);

protected:
    // utilize the ranking selection scheme
    ranking_selection<Chromosome,Encoding>* m_sel_scheme;

    // utilize a replace-worst environment selection scheme
    replace_worst<Chromosome,Encoding>* m_rep_scheme;

    // allow whatever crossover operator we see fit
    crossover_operator<Chromosome,Encoding>* m_cross_op;
    double m_cross_rate;

    // allow whatever mutation operator we see fit
    mutation_operator<Chromosome,Encoding>* m_mut_op;

    // define a comparator
    comparator<Chromosome,Encoding>* m_comp;
    
public:
    // ctors and dtor
    genitor();
    virtual ~genitor();

    // initialize the genitor algorithm
    virtual void initialize();

    // run the genitor algorithm
    virtual void run();

    // run a single generation of the genitor algorithm
    virtual void run_one_generation();
};

#include "genitor.cpp"

#endif 
