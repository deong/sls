/**
 * \file sandm.h
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _SANDM_H_
#define _SANDM_H_

#include "sls.h"
#include "pfront.h"
#include "mutation.h"

/**
 * \class sandm
 * \brief sample and mutate search algorithm
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
