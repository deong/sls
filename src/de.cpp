/**
 * @author Deon Garrett <jdgarrett@gmail.com>
 * @date 02/05/2014
 * 
 * Copyright (c) 2014 IIIM
 */

/**
 * @file de.cpp
 */

#include <iostream>
#include "kvparse/kvparse.h"
#include "keywords.h"

using namespace std;

template <template <typename> class Chromosome>
differential_evolution<Chromosome>::differential_evolution() :
	m_cross_rate(1.0),
	m_force(1.0),
	m_comp(0),
	m_sel_scheme(0)
{
}

template <template <typename> class Chromosome>
differential_evolution<Chromosome>::~differential_evolution()
{
	if(m_comp) {
		delete m_comp;
		m_comp = 0;
	}
	if(m_sel_scheme) {
		delete m_sel_scheme;
		m_sel_scheme = 0;
	}
}

/*!
 * \brief initialize the algorithm components
 */
template <template <typename> class Chromosome>
void differential_evolution<Chromosome>::initialize()
{
	// initialize the base class members
	ea<Chromosome,real_encoding>::initialize();

	// initialize the comparator
	comparator_factory<Chromosome,real_encoding> cf;
	m_comp = cf.construct();

	// set up the crossover rate
	kvparse::parameter_value(keywords::CROSSOVER_RATE, m_cross_rate, false);

	// set up the force paramater
	kvparse::parameter_value(keywords::DE_FORCE, m_force, false);
	
	// create the selection operator
	m_sel_scheme = new random_selection<Chromosome, real_encoding>();
}

/*!
 * \brief run the differential evolution algorithm
 */
template <template <typename> class Chromosome>
void differential_evolution<Chromosome>::run()
{
	while(!this->terminate()) {
		this->run_one_generation();
		this->generation_completed(this->m_population);
	}
	this->compute_metrics();
	this->report_metrics(cout);
}

/*!
 * \brief run one generation of the differential_evolution algorithm
 */
template <template <typename> class Chromosome>
void differential_evolution<Chromosome>::run_one_generation()
{
	mtrandom mt;
	m_sel_scheme->set_population(&this->m_population);

	for(unsigned int i=0; i<this->m_population.size(); ++i) 
	{
		// individual to alter
		Chromosome<real_encoding> x = this->m_population[i];
		
		// select three parents
		Chromosome<real_encoding> a = m_sel_scheme->select_parent();
		Chromosome<real_encoding> b = m_sel_scheme->select_parent();
		Chromosome<real_encoding> c = m_sel_scheme->select_parent();
	
		// select an index at random
		unsigned int R = mt.random(0, a.length());
	
		for(unsigned int j=0; j<a.length(); ++j) {
			if(mt.random() < m_cross_rate || j == R) {
				x[j] = a[j] + m_force * (b[j] - c[j]);
			}
		}
		m_repair.repair(x, this->m_fitfunc);
		
		// evaluate the offspring
		x.evaluate(this->m_fitfunc);
		this->chromosome_evaluated(x);

		if(m_comp->compare(x, this->m_population[i]) == -1) {
			this->m_population[i] = x;
		}
	}
}
