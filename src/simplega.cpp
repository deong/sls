/*!
 * \file simplega.cpp
 *
 * Simple Genetic Algorithm uses a specified selection scheme,
 * crossover operator, mutation operator, and generational or
 * elitist replacement.
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <iostream>
#include "simplega.h"
#include "chromosome.h"
#include "encoding.h"
#include "population.h"
#include "crossover.h"
#include "mutation.h"
#include "selection.h"
#include "replacement.h"
#include "mtrandom.h"
#include "configuration.h"
#include "keywords.h"
#include "utilities.h"

using namespace std;

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
simple_ga<Chromosome,Encoding>::simple_ga() :
	m_sel_scheme(0),
	m_rep_scheme(0),
	m_cross_op(0),
	m_cross_rate(1.0),
	m_mut_op(0)
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
simple_ga<Chromosome,Encoding>::~simple_ga()
{
	delete m_sel_scheme;
	delete m_rep_scheme;
	delete m_cross_op;
	delete m_mut_op;
	delete m_comp;
}

/*!
 * \brief initialize the simple ga components
 */
template <template <typename> class Chromosome, typename Encoding>
void simple_ga<Chromosome,Encoding>::initialize()
{
	// initialize the base class members
	ea<Chromosome,Encoding>::initialize();

	// initialize the comparator
	comparator_factory<Chromosome,Encoding> cf;
	m_comp = cf.construct();

	// sort the population by fitness
	this->m_population.sort(m_comp);

	// configure the selection and replacement schemes
	m_sel_scheme = selection_scheme_factory<Chromosome,Encoding>::construct();
	m_rep_scheme = replacement_scheme_factory<Chromosome,Encoding>::construct();

	// figure out which crossover operator to use
	m_cross_op = crossover_operator_factory<Chromosome,Encoding>::construct();

	// set up the crossover rate
	configuration::double_parameter(keywords::CROSSOVER_RATE, m_cross_rate, false);

	// figure out which mutation operator to use
	m_mut_op = mutation_operator_factory<Chromosome,Encoding>::construct();
}

/*!
 * \brief run the simple ga
 */
template <template <typename> class Chromosome, typename Encoding>
void simple_ga<Chromosome,Encoding>::run()
{
	while(!this->terminate()) {
		this->run_one_generation();
		this->generation_completed(this->m_population);
	}
	this->compute_metrics();
	this->report_metrics(cout);
}

/*!
 * \brief run one generation of the simple ga
 */
template <template <typename> class Chromosome, typename Encoding>
void simple_ga<Chromosome,Encoding>::run_one_generation()
{
	mtrandom mt;
	population<Chromosome,Encoding> offspring(this->m_population.size());

	m_sel_scheme->set_population(&this->m_population);

	for(unsigned int i=0; i<this->m_population.size()/2; i++) {
		// select two parents
		Chromosome<Encoding> p1 = m_sel_scheme->select_parent();
		Chromosome<Encoding> p2 = m_sel_scheme->select_parent();

		// perform crossover with some probability
		Chromosome<Encoding> c1 = p1;
		Chromosome<Encoding> c2 = p2;
		if(mt.random() < m_cross_rate) {
			m_cross_op->crossover(p1, p2, c1, c2);
		}

		// mutate the offspring
		m_mut_op->mutate(c1);
		m_mut_op->mutate(c2);

		// evaluate the offspring
		c1.evaluate(this->m_fitfunc);
		c2.evaluate(this->m_fitfunc);
		this->chromosome_evaluated(c1);
		this->chromosome_evaluated(c2);

		offspring.add(c1);
		offspring.add(c2);
	}

	// merge the offspring into the population according to the specified
	// replacement scheme
	m_rep_scheme->merge_populations(this->m_population, offspring, this->m_population);
}
