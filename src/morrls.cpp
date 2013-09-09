/*!
 * \file morrls.cpp
 *
 * random restarts local search
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include "morrls.h"
#include "sls.h"
#include "localsearch.h"
#include "metrics.h"

using namespace std;

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
morrls<Chromosome,Encoding>::morrls() :
	sls<Chromosome,Encoding>::sls(),
	mu(0),
	m_iter(0),
	m_hc(0)
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
morrls<Chromosome,Encoding>::~morrls()
{
	delete m_hc;
}

/*!
 * \brief notify of generation completion
 */
template <template <typename> class Chromosome, typename Encoding>
void morrls<Chromosome,Encoding>::generation_completed()
{
	sls<Chromosome,Encoding>::generation_completed(m_front);
}

/*!
 * \brief initialize local search components
 */
template <template <typename> class Chromosome, typename Encoding>
void morrls<Chromosome,Encoding>::initialize()
{
	// initialize the stochastic local search components
	sls<Chromosome,Encoding>::initialize();

	// initialize the local search components
	local_search_factory<Chromosome,Encoding> lsf;
	lsf.set_prefix("embedded_");
	m_hc = lsf.construct();

	// set up the number of different weight vectors to use
	configuration::unsigned_integer_parameter(keywords::NUM_WEIGHTS, mu, true);
}

/*!
 * \brief iteratively run the local search until termination
 */
template <template <typename> class Chromosome, typename Encoding>
void morrls<Chromosome,Encoding>::run()
{
	while(!this->terminate()) {
		for(unsigned int index=0; index<mu; index++) {
			m_comp.randomize_weights(this->m_fitfunc->objectives());
			// m_comp.weights[0] = double(mu-index-1)/(mu-1);
			// m_comp.weights[1] = 1.0 - m_comp.weights[0];
			Chromosome<Encoding> chr(this->m_fitfunc);
			chr.randomize();
			chr.evaluate(this->m_fitfunc);
			if(this->m_repair) {
				this->m_repair->repair(chr,this->m_fitfunc);
			}
			this->chromosome_evaluated(chr);
			m_hc->improve(chr, &m_comp, this->m_fitfunc);
			m_hc->reset();
			m_front.add(chr);
		}
		this->generation_completed();
	}

	cout << m_front << endl;
	this->compute_metrics();
	this->report_metrics(cout);
}

