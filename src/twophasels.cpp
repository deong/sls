/*!
 * \file twophasels.cpp
 *
 * two phase local search
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include "twophasels.h"
#include "sls.h"
#include "localsearch.h"
#include "metrics.h"

using namespace std;

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
twophasels<Chromosome,Encoding>::twophasels() :
	sls<Chromosome,Encoding>::sls(),
	m_phase_two_runs(0),
	m_hc(0),
	m_hc2(0)
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
twophasels<Chromosome,Encoding>::~twophasels()
{
	delete m_hc;
	delete m_hc2;
}

/*!
 * \brief notify of generation completion
 */
template <template <typename> class Chromosome, typename Encoding>
void twophasels<Chromosome,Encoding>::generation_completed()
{
	sls<Chromosome,Encoding>::generation_completed(m_front);
}

/*!
 * \brief initialize local search components
 */
template <template <typename> class Chromosome, typename Encoding>
void twophasels<Chromosome,Encoding>::initialize()
{
	// initialize the stochastic local search components
	sls<Chromosome,Encoding>::initialize();

	// initialize the local search components
	local_search_factory<Chromosome,Encoding> lsf1;
	lsf1.set_prefix("ls_");
	m_hc = lsf1.construct();

	local_search_factory<Chromosome,Encoding> lsf2;
	lsf2.set_prefix("ls2_");
	m_hc2 = lsf2.construct();

	// how many different weight vectors to try in phase two
	configuration::unsigned_integer_parameter(keywords::PHASE_TWO_RUNS,m_phase_two_runs,true);
}

/*!
 * \brief iteratively run the local search until termination
 */
template <template <typename> class Chromosome, typename Encoding>
void twophasels<Chromosome,Encoding>::run()
{
	while(!this->terminate()) {
		// create a new solution via random weight vector
		m_comp.randomize_weights(this->m_fitfunc->objectives());
		Chromosome<Encoding> chr(this->m_fitfunc);
		chr.randomize();
		if(this->m_repair) {
			this->m_repair->repair(chr,this->m_fitfunc);
		}
		chr.evaluate(this->m_fitfunc);
		this->chromosome_evaluated(chr);
		m_hc->improve(chr, &m_comp, this->m_fitfunc);
		m_front.add(chr);

		// randomly mutate the weight vector and start phase two
		for(unsigned int i=0; i<m_phase_two_runs; i++) {
			// update weight vector by random amount
			mtrandom mt;
			double update=mt.gaussian(0.1,0.025);
			if(mt.random()<0.5) {
				update=-update;
			}

			// make sure the update keeps the weight vector valid
			if(m_comp.weights[0]+update>1.0) {
				update=-update;
			} else if(m_comp.weights[0]+update<0.0) {
				update=-update;
			}

			m_comp.weights[0]+=update;
			m_comp.weights[1]-=update;

			Chromosome<Encoding> chr2=chr;
			m_hc2->improve(chr2,&m_comp,this->m_fitfunc);
			m_front.add(chr2);
		}
		this->generation_completed();
	}

	cout << m_front << endl;
	this->compute_metrics();
	this->report_metrics(cout);
}

