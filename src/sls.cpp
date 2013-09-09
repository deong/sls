/*!
 * \file sls.cpp
 *
 * base class for all stochastic local search algorithms
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <list>
#include "sls.h"
#include "encoding.h"
#include "problems.h"
#include "metrics.h"
#include "configuration.h"
#include "keywords.h"
#include "utilities.h"

using namespace std;

template <template <typename> class Chromosome, typename Encoding>
list<metric<Chromosome,Encoding>*> sls<Chromosome,Encoding>::m_metrics;

template <template <typename> class Chromosome, typename Encoding>
list<terminator<Chromosome,Encoding>*> sls<Chromosome,Encoding>::m_terminators;

template <template <typename> class Chromosome, typename Encoding>
bool sls<Chromosome,Encoding>::m_print_every_generation;

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
sls<Chromosome,Encoding>::sls() :
	m_fitfunc(0),
	m_repair(0)
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
sls<Chromosome,Encoding>::~sls()
{
	if(this->m_fitfunc) {
		delete this->m_fitfunc;
	}
	if(this->m_repair) {
		delete this->m_repair;
	}

	Encoding::clear_parameters();
	for(typename list<terminator<Chromosome,Encoding>*>::iterator it=this->m_terminators.begin();
	        it!=this->m_terminators.end();
	        it++) {
		delete *it;
	}
	for(typename list<metric<Chromosome,Encoding>*>::iterator it=this->m_metrics.begin();
	        it!=this->m_metrics.end();
	        it++) {
		delete *it;
	}
}

/*!
 * \brief notify metrics and terminators of new chromosome evaluation
 */
template <template <typename> class Chromosome, typename Encoding>
void sls<Chromosome,Encoding>::chromosome_evaluated(const Chromosome<Encoding>& chr)
{
	for(typename list<terminator<Chromosome,Encoding>*>::iterator it=m_terminators.begin();
	        it!=m_terminators.end();
	        it++) {
		(*it)->chromosome_evaluated(chr);
	}
	for(typename list<metric<Chromosome,Encoding>*>::iterator it=m_metrics.begin();
	        it!=m_metrics.end();
	        it++) {
		(*it)->chromosome_evaluated(chr);
	}
}

/*!
 * \brief notify metrics and terminators of generation completion
 */
template <template <typename> class Chromosome, typename Encoding>
void sls<Chromosome,Encoding>::generation_completed()
{
	for(typename list<terminator<Chromosome,Encoding>*>::iterator it=m_terminators.begin();
	        it!=m_terminators.end();
	        it++) {
		(*it)->generation_completed();
	}
	for(typename list<metric<Chromosome,Encoding>*>::iterator it=m_metrics.begin();
	        it!=m_metrics.end();
	        it++) {
		(*it)->generation_completed();
	}
	if(m_print_every_generation) {
		report_metrics(cout);
		cout << endl;
	}
}

/*!
 * \brief notify metrics and terminators of new generation completion
 */
template <template <typename> class Chromosome, typename Encoding>
void sls<Chromosome,Encoding>::generation_completed(const population<Chromosome,Encoding>& pop)
{
	for(typename list<terminator<Chromosome,Encoding>*>::iterator it=m_terminators.begin();
	        it!=m_terminators.end();
	        it++) {
		(*it)->generation_completed(pop);
	}
	for(typename list<metric<Chromosome,Encoding>*>::iterator it=m_metrics.begin();
	        it!=m_metrics.end();
	        it++) {
		(*it)->generation_completed(pop);
	}
	if(m_print_every_generation) {
		report_metrics(cout);
		cout << endl;
	}
}

/*!
 * \brief determine whether or not to terminate the algorithm
 */
template <template <typename> class Chromosome, typename Encoding>
bool sls<Chromosome,Encoding>::terminate()
{
	for(typename list<terminator<Chromosome,Encoding>*>::const_iterator it=m_terminators.begin();
	        it != m_terminators.end();
	        it++) {
		if((*it)->terminate()) {
			return true;
		}
	}
	return false;
}

/*!
 * \brief compute all performance metrics
 */
template <template <typename> class Chromosome, typename Encoding>
void sls<Chromosome,Encoding>::compute_metrics()
{
	for(typename list<metric<Chromosome,Encoding>*>::iterator it=m_metrics.begin();
	        it != m_metrics.end();
	        it++) {
		(*it)->compute();
	}
}

/*!
 * \brief print out all performance metrics
 */
template <template <typename> class Chromosome, typename Encoding>
void sls<Chromosome,Encoding>::report_metrics(ostream& ostr)
{
	for(typename list<metric<Chromosome,Encoding>*>::const_iterator it=m_metrics.begin();
	        it!=m_metrics.end();
	        it++) {
		(*it)->report(ostr);
	}
}

/*!
 * \brief initialize the sls components
 */
template <template <typename> class Chromosome, typename Encoding>
void sls<Chromosome,Encoding>::initialize()
{
	// read in the objective function information
	this->m_fitfunc = Encoding::ProblemFactoryType::construct();

	// initialize the encoding parameters
	Encoding::initialize_parameters(m_fitfunc);

	// get and construct the termination criteria
	terminator_factory<Chromosome,Encoding> tf;
	m_terminators = tf.construct();

	// initialize the performance metrics
	metric_factory<Chromosome,Encoding> mf;
	m_metrics = mf.construct();

	// configure the repair operator
	if(configuration::keyword_exists(keywords::REPAIR_OPERATOR)) {
		repair_factory<Chromosome,Encoding> rof;
		this->m_repair=rof.construct();
	}

	// verbose logging
	m_print_every_generation = false;
	configuration::boolean_parameter(keywords::PRINT_EVERY_GENERATION, m_print_every_generation, false);
}
