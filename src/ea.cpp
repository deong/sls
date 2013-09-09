/*!
 * \file ea.cpp
 *
 * all evolutionary algorithms share certain aspects in common.  each
 * maintains a population of candidate chromosomes, using defined genetic
 * operators to broaden the search and selection based on fitness to focus
 * the search.  this class provides the basic infrastructure for all sorts
 * of evolutionary algorithms.
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <iostream>
#include <list>
#include "ea.h"
#include "sls.h"
#include "problems.h"
#include "chromosome.h"
#include "population.h"
#include "encoding.h"
#include "metrics.h"
#include "terminators.h"
#include "convergence.h"
#include "configuration.h"
#include "keywords.h"

using namespace std;

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
ea<Chromosome,Encoding>::ea()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
ea<Chromosome,Encoding>::~ea()
{
	for(typename list<convergence<Chromosome,Encoding>*>::iterator it=this->m_convergence.begin();
	        it!=this->m_convergence.end();
	        it++) {
		delete *it;
	}
}

/*!
 * \brief notify of new chromosome evaluation
 */
template <template <typename> class Chromosome, typename Encoding>
void ea<Chromosome,Encoding>::chromosome_evaluated(const Chromosome<Encoding>& sol)
{
	sls<Chromosome,Encoding>::chromosome_evaluated(sol);

	for(typename list<convergence<Chromosome,Encoding>*>::iterator it=this->m_convergence.begin();
	        it!=this->m_convergence.end();
	        it++) {
		(*it)->chromosome_evaluated(sol);
	}
}

/*!
 * \brief notify of generation completion
 */
template <template <typename> class Chromosome, typename Encoding>
void ea<Chromosome,Encoding>::generation_completed(const population<Chromosome,Encoding>& pop)
{
	sls<Chromosome,Encoding>::generation_completed(pop);

	for(typename list<convergence<Chromosome,Encoding>*>::iterator it=this->m_convergence.begin();
	        it!=this->m_convergence.end();
	        it++) {
		(*it)->generation_completed(pop);
	}
}

/*!
 * \brief check whether algorithm has converged
 */
template <template <typename> class Chromosome, typename Encoding>
bool ea<Chromosome,Encoding>::converged() const
{
	for(typename list<convergence<Chromosome,Encoding>*>::const_iterator it=this->m_convergence.begin();
	        it!=this->m_convergence.end();
	        it++) {
		if((*it)->converged()) {
			return true;
		}
	}
	return false;
}

/*!
 * \brief initialize the generic ea components
 */
template <template <typename> class Chromosome, typename Encoding>
void ea<Chromosome,Encoding>::initialize()
{
	// initialize the parent class
	sls<Chromosome,Encoding>::initialize();

	// get and construct the convergence criteria
	this->m_convergence = convergence_factory<Chromosome,Encoding>::construct();

	// read in the population size
	int popsize;
	configuration::integer_parameter(keywords::POPULATION_SIZE, popsize, true);

	// construct the population
	for(int i=0; i<popsize; i++) {
		Chromosome<Encoding> sol(this->m_fitfunc);
		sol.randomize();
		if(this->m_repair) {
			this->m_repair->repair(sol,this->m_fitfunc);
		}
		sol.evaluate(this->m_fitfunc);
		chromosome_evaluated(sol);
		this->m_population.add(sol);
	}
	generation_completed(this->m_population);
}
