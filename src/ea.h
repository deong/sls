/*!
 * \file ea.h
 *
 * all evolutionary algorithms share certain aspects in common.  each
 * maintains a population of candidate chromosomes, using defined genetic
 * operators to broaden the search and selection based on fitness to focus
 * the search.  this class provides the basic infrastructure for all sorts
 * of evolutionary algorithms.
 *
 * Deon Garrett
 * deong@acm.org
 */

#ifndef _EA_H_
#define _EA_H_

#include <iostream>
#include <list>
#include "sls.h"
#include "chromosome.h"
#include "population.h"
#include "convergence.h"

using namespace std;

/*!
 * \class ea
 *
 * base class for evolutionary algorithms
 *
 * \author deong
 * \date 05/09/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class ea : virtual public sls<Chromosome,Encoding>
{
private:
    // disable copying of heavyweight class
    ea(const ea& that);
    ea& operator=(const ea& that);
    
protected:
    // store the individuals in the population
    population<Chromosome,Encoding> m_population;

    // maintain a list of convergence criteria
    list<convergence<Chromosome,Encoding>*> m_convergence;
    
public:
    ea();
    virtual ~ea();

    // update the performance metrics and termination criteria
    virtual void chromosome_evaluated(const Chromosome<Encoding>& sol);
    virtual void generation_completed(const population<Chromosome,Encoding>& pop);

    // determine whether or not the algorithm has converged
    virtual bool converged() const;
    
    // set up the ea parameters
    virtual void initialize();

    // run the evolutionary algorithm
    virtual void run() = 0;

    // run one generation of the evolutionary algorithm
    virtual void run_one_generation() = 0;
};

#include "ea.cpp"

#endif    
