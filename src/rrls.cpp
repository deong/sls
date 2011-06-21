/*!
 * \file rrls.cpp
 *
 * random restarts local search
 *
 * Deon Garrett
 * deong@acm.org
 */

#include "rrls.h"
#include "sls.h"
#include "localsearch.h"

using namespace std;

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
rrls<Chromosome,Encoding>::rrls() :
    sls<Chromosome,Encoding>::sls(),
    m_hc(0),
    m_comp(0)
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
rrls<Chromosome,Encoding>::~rrls()
{
    delete m_hc;
    delete m_comp;
}

/*!
 * \brief initialize the local search components
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
void rrls<Chromosome,Encoding>::initialize()
{
    // initialize the stochastic local search components
    sls<Chromosome,Encoding>::initialize();

    // initialize the comparator
    comparator_factory<Chromosome,Encoding> cf;
    m_comp = cf.construct();
    
    // initialize the local search components
    local_search_factory<Chromosome,Encoding> lsf;
    lsf.set_prefix("ls_");
    m_hc = lsf.construct();
}

/*!
 * \brief run the random restart local search until termination
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
void rrls<Chromosome,Encoding>::run()
{
    while(!this->terminate())
    {
        Chromosome<Encoding> chr(this->m_fitfunc);
        chr.randomize();
        chr.evaluate(this->m_fitfunc);
        this->chromosome_evaluated(chr);
        m_hc->improve(chr, m_comp, this->m_fitfunc);
        this->generation_completed();
    }
    this->compute_metrics();
    this->report_metrics(cout);
}

