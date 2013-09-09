/*!
 * \file chc.h
 *
 * Simple Genetic Algorithm (generational ga)
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _CHC_H_
#define _CHC_H_

#include "ea.h"
#include "selection.h"
#include "replacement.h"
#include "crossover.h"
#include "mutation.h"
#include "comparator.h"
#include "neighborhood.h"

/*!
 * \class chc
 *
 * Eshelman's CHC genetic algorithm
 */
template <template <typename> class Chromosome, typename Encoding>
class chc : public ea<Chromosome,Encoding>
{
public:
	chc();
	virtual ~chc();
	virtual void initialize();
	virtual void run();
	virtual void run_one_generation();
	void diverge();

protected:
	random_selection<Chromosome,Encoding>* m_sel_scheme;
	truncation_replacement<Chromosome,Encoding>* m_rep_scheme;
	hux_crossover<Chromosome,Encoding>* m_cross_op;
	comparator<Chromosome,Encoding>* m_comp;
	neighborhood<Chromosome,Encoding>* m_distmeas;
	double m_difference_threshold;
	double m_divergence_rate;

private:
	// disable copying of heavyweight class
	chc(const chc& that);
	chc& operator=(const chc& that);
};

#include "chc.cpp"

#endif
