/*!
 * \file sls.cpp
 *
 * base class for all stochastic local search algorithms
 *
 * Deon Garrett
 * deong@acm.org
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

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/11/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
sls<Chromosome,Encoding>::sls() :
    m_fitfunc(0),
    m_repair(0)
{
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/11/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
sls<Chromosome,Encoding>::~sls()
{
    if(this->m_fitfunc)
		delete this->m_fitfunc;
    if(this->m_repair)
		delete this->m_repair;
    
    Encoding::clear_parameters();
    for(typename list<terminator<Chromosome,Encoding>*>::iterator it=this->m_terminators.begin();
        it!=this->m_terminators.end();
        it++)
    {
        delete *it;
    }
    for(typename list<metric<Chromosome,Encoding>*>::iterator it=this->m_metrics.begin();
        it!=this->m_metrics.end();
        it++)
    {
        delete *it;
    }
}

/*!
 * \brief notify metrics and terminators of new chromosome evaluation
 *
 * \author deong
 * \date 05/11/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void sls<Chromosome,Encoding>::chromosome_evaluated(const Chromosome<Encoding>& chr)
{
    for(typename list<terminator<Chromosome,Encoding>*>::iterator it=m_terminators.begin();
        it!=m_terminators.end();
        it++)
    {
        (*it)->chromosome_evaluated(chr);
    }
    for(typename list<metric<Chromosome,Encoding>*>::iterator it=m_metrics.begin();
        it!=m_metrics.end();
        it++)
    {
        (*it)->chromosome_evaluated(chr);
    }
	//cout << chr << endl;
}

/*!
 * \brief notify metrics and terminators of generation completion
 *
 * \author deong
 * \date 05/11/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void sls<Chromosome,Encoding>::generation_completed()
{
    for(typename list<terminator<Chromosome,Encoding>*>::iterator it=m_terminators.begin();
        it!=m_terminators.end();
        it++)
    {
        (*it)->generation_completed();
    }
    for(typename list<metric<Chromosome,Encoding>*>::iterator it=m_metrics.begin();
        it!=m_metrics.end();
        it++)
    {
        (*it)->generation_completed();
    }
}

/*!
 * \brief notify metrics and terminators of new generation completion
 *
 * \author deong
 * \date 05/12/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void sls<Chromosome,Encoding>::generation_completed(const population<Chromosome,Encoding>& pop)
{
    for(typename list<terminator<Chromosome,Encoding>*>::iterator it=m_terminators.begin();
        it!=m_terminators.end();
        it++)
    {
        (*it)->generation_completed(pop);
    }
    for(typename list<metric<Chromosome,Encoding>*>::iterator it=m_metrics.begin();
        it!=m_metrics.end();
        it++)
    {
        (*it)->generation_completed(pop);
    }
}

/*!
 * \brief determine whether or not to terminate the algorithm
 *
 * \author deong
 * \date 06/25/2007
 */
template <template <typename> class Chromosome, typename Encoding>
bool sls<Chromosome,Encoding>::terminate()
{
    for(typename list<terminator<Chromosome,Encoding>*>::const_iterator it=m_terminators.begin();
	it != m_terminators.end();
	it++)
    {
	if((*it)->terminate())
	{
	    return true;
	}
    }
    return false;
}

/*!
 * \brief compute all performance metrics
 *
 * \author deong
 * \date 05/11/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void sls<Chromosome,Encoding>::compute_metrics()
{
    for(typename list<metric<Chromosome,Encoding>*>::iterator it=m_metrics.begin();
        it != m_metrics.end();
        it++)
    {
        (*it)->compute();
    }
}

/*!
 * \brief print out all performance metrics
 *
 * \author deong
 * \date 05/11/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void sls<Chromosome,Encoding>::report_metrics(ostream& ostr)
{
    for(typename list<metric<Chromosome,Encoding>*>::const_iterator it=m_metrics.begin();
        it!=m_metrics.end();
        it++)
    {
        (*it)->report(ostr);
    }
}

/*!
 * \brief initialize the sls components
 *
 * \author deong
 * \date 05/11/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
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
    if(configuration::keyword_exists(keywords::REPAIR_OPERATOR))
    {
		repair_factory<Chromosome,Encoding> rof;
		this->m_repair=rof.construct();
    }
}
