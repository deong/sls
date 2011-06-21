#ifndef _SANDM_H_
#define _SANDM_H_

/**
 * \file sandm.h
 *
 * \author deong@acm.org
 * \date 01/17/2008
 */

#include "sls.h"
#include "pfront.h"
#include "mutation.h"

/**
 * \class sandm
 * \brief sample and mutate search algorithm
 *
 * \author deong@acm.org
 * \date 01/17/2008
 */
template <template <typename> class Chromosome, typename Encoding>
class sandm : public sls<Chromosome,Encoding>
{
public:
    sandm();
    virtual ~sandm();
    virtual void initialize();
    virtual void run();

protected:
    unsigned int m_sample_size;
    pareto_front<Chromosome,Encoding> m_front;
    mutation_operator<Chromosome,Encoding>* m_mut_op;
};

#include "sandm.cpp"

#endif
