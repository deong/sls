/*!
 * \file simanneal.h
 *
 * Simulated Annealing algorithm
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _SIMANNEAL_H_
#define _SIMANNEAL_H_

#include "localsearch.h"
#include "chromosome.h"
#include "comparator.h"
#include "problems.h"
#include "metrics.h"
#include "terminators.h"

/*!
 * \class simulated_annealing
 */
template <template <typename> class Chromosome, typename Encoding>
class simulated_annealing : public local_search<Chromosome,Encoding>
{
private:
    // disable copy
    simulated_annealing(const simulated_annealing& that);
    simulated_annealing& operator=(const simulated_annealing& that);

public:
    simulated_annealing();
    virtual ~simulated_annealing();
    virtual void initialize();
    virtual void reset();
    virtual void improve(Chromosome<Encoding>& chr,
			 comparator<Chromosome,Encoding>* comp,
			 const typename Encoding::ProblemType* prob);

private:
    bool accept_move(double delta) const;
    
private:
    double m_alpha;
    double m_temp;
    unsigned int m_schedule;
};

#include "simanneal.cpp"

#endif
