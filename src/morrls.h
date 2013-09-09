/*!
 * \file morrls.h
 *
 * multiobjective random restarts local search
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _MORRLS_H_
#define _MORRLS_H_

#include "sls.h"
#include "localsearch.h"
#include "pfront.h"
#include "comparator.h"

using namespace std;

/*!
 * \class morrls
 *
 * multiobjective random restarts local search
 */
template <template <typename> class Chromosome, typename Encoding>
class morrls : virtual public sls<Chromosome,Encoding>
{
protected:
    unsigned int mu;
    unsigned int m_iter;
    local_search<Chromosome,Encoding>* m_hc;
    scalarizing_comparator<Chromosome,Encoding> m_comp;
    pareto_front<Chromosome,Encoding> m_front;
    
private:
    // disable copying
    morrls(const morrls& that);
    morrls& operator=(const morrls& that);
    
public:
    morrls();
    virtual ~morrls();

    virtual void generation_completed();
    virtual void initialize();
    virtual void run();
};

#include "morrls.cpp"

#endif
