/*!
 * \file convergence.h
 *
 * classes which are used to determine whether an evolutionary
 * algorithm has converged
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _CONVERGENCE_H_
#define _CONVERGENCE_H_

#include <list>
#include "chromosome.h"

/*!
 * \class convergence
 *
 * monitor the convergence behavior of a search algorithm
 */
template <template <typename> class Chromosome, typename Encoding>
class convergence
{
private:
	// disable copying of convergence measures
	convergence(const convergence& that);
	convergence& operator=(const convergence& that);

public:
	convergence();
	virtual ~convergence();

	virtual void chromosome_evaluated(const Chromosome<Encoding>& sol);
	virtual void generation_completed(const population<Chromosome,Encoding>& pop);

	virtual bool converged() const = 0;
};

/*!
 * \class convergence_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class convergence_factory
{
public:
	static list<convergence<Chromosome,Encoding>*> construct();
};

#include "convergence.cpp"

#endif
