/**
 * @author Deon Garrett <jdgarrett@gmail.com>
 * @date 02/05/2014
 * 
 * Copyright (c) 2014 IIIM
 */

#ifndef _DE_H_
#define _DE_H_

/**
 * @file de.h
 */

#include "ea.h"
#include "chromosome.h"
#include "encoding.h"
#include "comparator.h"
#include "selection.h"
#include "repair.h"

template <template <typename> class Chromosome>
class differential_evolution : public ea<Chromosome, real_encoding> 
{
private:
    // disable copying of heavyweight class
	differential_evolution(const differential_evolution& that);
	differential_evolution& operator=(const differential_evolution& that);

protected:
	// crossover rate
	double m_cross_rate;
	
	// force constant determining contribution of A vs B-C in crossover
	double m_force;
	
	// define a comparator
	comparator<Chromosome,real_encoding>* m_comp;

	// DE uses random selection
	random_selection<Chromosome,real_encoding>* m_sel_scheme;
	
	// keep the parameters in bounds
	clip_bounds_repair<Chromosome,real_encoding> m_repair;
	
public:
	differential_evolution();
	virtual ~differential_evolution();
	virtual void initialize();
	virtual void run();
	virtual void run_one_generation();
};

#include "de.cpp"

#endif
