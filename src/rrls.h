/*!
 * \file rrls.h
 *
 * random restarts local search
 *
 * Deon Garrett
 * deong@acm.org
 */

#ifndef _RRLS_H_
#define _RRLS_H_

#include "sls.h"
#include "localsearch.h"

/*!
 * \class rrls
 *
 * random restarts local search
 *
 * \author deong
 * \date 05/11/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class rrls : virtual public sls<Chromosome,Encoding>
{
protected:
    local_search<Chromosome,Encoding>* m_hc;
    comparator<Chromosome,Encoding>* m_comp;
    
private:
    // disable copying
    rrls(const rrls& that);
    rrls& operator=(const rrls& that);
    
public:
    rrls();
    virtual ~rrls();

    virtual void initialize();
    virtual void run();
};

#include "rrls.cpp"

#endif
