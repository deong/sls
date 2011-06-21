/*!
 * \file sls.h
 *
 * base class for stochastic local search algorithms
 *
 * Deon Garrett
 * deong@acm.org
 */

#ifndef _SLS_H_
#define _SLS_H_

#include <list>
#include "chromosome.h"
#include "comparator.h"
#include "metrics.h"
#include "terminators.h"
#include "problems.h"
#include "repair.h"

using namespace std;

/*!
 * \class sls
 *
 * base class for all stochastic local search algorithms
 *
 * \author deong
 * \date 05/11/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class sls 
{
private:
    // disable copying of heavyweight class
    sls(const sls& that);
    sls& operator=(const sls& that);

protected:
    typename Encoding::ProblemType* m_fitfunc;
    repair_operator<Chromosome,Encoding>* m_repair;
    static list<metric<Chromosome,Encoding>*> m_metrics;
    static list<terminator<Chromosome,Encoding>*> m_terminators;
    
public:
    sls();
    virtual ~sls();
    static void chromosome_evaluated(const Chromosome<Encoding>& sol);
    static void generation_completed();
    static void generation_completed(const population<Chromosome,Encoding>& pop);
    static void compute_metrics();
    static void report_metrics(ostream& ostr);
    static bool terminate();
    virtual void initialize();
    virtual void run() = 0;
};

#include "sls.cpp"

#endif
