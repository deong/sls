/*!
 * \file pfront.h
 *
 * class representing a pareto front
 *
 * Deon Garrett
 * deong@acm.org
 */

#ifndef _PFRONT_H_
#define _PFRONT_H_

#include <iostream>
#include "population.h"
#include "chromosome.h"
#include "problems.h"

using namespace std;

template <template <typename> class Chromosome, typename Encoding> class pareto_front;
template <template <typename> class Chromosome, typename Encoding>
ostream& operator<<(ostream& s, const pareto_front<Chromosome,Encoding>& p);

/*!
 * \class pareto_front
 *
 * maintains a set of mutually nondominated solutions
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class pareto_front : public population<Chromosome,Encoding>
{
public:
    pareto_front();
    pareto_front(unsigned int sz);
    pareto_front(const pareto_front& pf);
    pareto_front(const population<Chromosome,Encoding>& pop);
    virtual ~pareto_front();

    // overload the add method to respect pareto dominance
    virtual bool add(const Chromosome<Encoding>& chr);

    void construct_front(const population<Chromosome,Encoding>& pop);
    void construct_front(string filename, const typename Encoding::ProblemType* prob);
    
    friend ostream& operator<<<>(ostream& ostr, const population<Chromosome,Encoding>& p);
};

#include "pfront.cpp"

#endif

