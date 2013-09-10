/*!
 * \file localsearch.cpp
 *
 * local search algorithms
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include "chromosome.h"
#include "encoding.h"
#include "problems.h"
#include "neighborhood.h"
#include "comparator.h"
#include "tabu.h"
#include "simanneal.h"
#include "metrics.h"
#include "terminators.h"
#include "kvparse/kvparse.h"
#include "keywords.h"

using namespace std;

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
local_search<Chromosome,Encoding>::local_search() :
	m_nf(0),
	m_repair(0),
	debug_generations(false)
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
local_search<Chromosome,Encoding>::~local_search()
{
	if(m_nf) {
		delete m_nf;
	}
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

	if(m_repair) {
		delete m_repair;
	}
}

/*!
 * \brief initialize local search components
 */
template <template <typename> class Chromosome, typename Encoding>
void local_search<Chromosome,Encoding>::initialize()
{
	m_nf = neighborhood_factory<Chromosome,Encoding>::construct();

	terminator_factory<Chromosome,Encoding> tf;
	tf.set_prefix(m_prefix);
	m_terminators = tf.construct();

	metric_factory<Chromosome,Encoding> mf;
	mf.set_prefix(m_prefix);
	m_metrics = mf.construct();

	repair_factory<Chromosome,Encoding> rf;
	rf.set_prefix(m_prefix);
	m_repair=rf.construct();

	kvparse::parameter_value(keywords::DEBUG_LS_GENERATIONS,debug_generations);
}

/*!
 * \brief reset the terminators and metrics for the algorithm
 */
template <template <typename> class Chromosome, typename Encoding>
void local_search<Chromosome,Encoding>::reset()
{
	for(typename list<metric<Chromosome,Encoding>*>::iterator it=this->m_metrics.begin();
	        it!=m_metrics.end();
	        it++) {
		(*it)->reset();
	}
	for(typename list<terminator<Chromosome,Encoding>*>::iterator it=this->m_terminators.begin();
	        it!=m_terminators.end();
	        it++) {
		(*it)->reset();
	}
}

/*!
 * \brief return the local search neighborhood operator
 */
template <template <typename> class Chromosome, typename Encoding>
neighborhood<Chromosome,Encoding>* local_search<Chromosome,Encoding>::get_neighborhood()
{
	return m_nf;
}

/*!
 * \brief notify of new chromosome evaluation
 */
template <template <typename> class Chromosome, typename Encoding>
void local_search<Chromosome,Encoding>::chromosome_evaluated(const Chromosome<Encoding>& chr)
{
	sls<Chromosome,Encoding>::chromosome_evaluated(chr);

	for(typename list<terminator<Chromosome,Encoding>*>::iterator it=this->m_terminators.begin();
	        it!=this->m_terminators.end();
	        it++) {
		(*it)->chromosome_evaluated(chr);
	}
	for(typename list<metric<Chromosome,Encoding>*>::iterator it=this->m_metrics.begin();
	        it!=this->m_metrics.end();
	        it++) {
		(*it)->chromosome_evaluated(chr);
	}
}

/*!
 * \brief notify of generation completion
 *
 * Note that a generation is basically defined as all moves leading
 * up to an accepted move.
 */
template <template <typename> class Chromosome, typename Encoding>
void local_search<Chromosome,Encoding>::generation_completed()
{
	//! we don't propagate generation_completed signals to the global scope
	//! as a matter of principle -- unlike evaluations, generations are not
	//! intended to be global.  each component algorithm maintains a wholely
	//! separate generation count
	//! sls<Chromosome,Encoding>::generation_completed();

	for(typename list<terminator<Chromosome,Encoding>*>::iterator it=this->m_terminators.begin();
	        it!=this->m_terminators.end();
	        it++) {
		(*it)->generation_completed();
	}
	for(typename list<metric<Chromosome,Encoding>*>::iterator it=this->m_metrics.begin();
	        it!=this->m_metrics.end();
	        it++) {
		(*it)->generation_completed();
		if(debug_generations) {
			(*it)->compute();
			(*it)->report(cout);
		}
	}
}

/**
 * \brief set the configuration prefix for the local search algorithm
 */
template <template <typename> class Chromosome, typename Encoding>
void local_search<Chromosome,Encoding>::set_prefix(const string& prefix)
{
	m_prefix=prefix;
}

/*!
 * \brief determine whether to terminate the local search procedure
 */
template <template <typename> class Chromosome, typename Encoding>
bool local_search<Chromosome,Encoding>::terminate() const
{
	if(sls<Chromosome,Encoding>::terminate()) {
		return true;
	}

	for(typename list<terminator<Chromosome,Encoding>*>::const_iterator it=this->m_terminators.begin();
	        it != m_terminators.end();
	        it++) {
		if((*it)->terminate()) {
			return true;
		}
	}

	return false;
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
next_descent<Chromosome,Encoding>::next_descent()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
next_descent<Chromosome,Encoding>::~next_descent()
{
}

/*!
 * \brief improve a given chromosome using next descent local search
 */
template <template <typename> class Chromosome, typename Encoding>
void next_descent<Chromosome,Encoding>::improve(Chromosome<Encoding>& chr,
        comparator<Chromosome,Encoding>* comp,
        const typename Encoding::ProblemType* prob)
{
	this->reset();

	bool at_local_opt = false;
	while(!at_local_opt && !this->terminate()) {
		Chromosome<Encoding> saved=chr;
		at_local_opt = true;
		this->m_nf->initialize(chr);
		while(this->m_nf->has_more_neighbors()) {
			move<Chromosome,Encoding> m=this->m_nf->next_neighbor();
			m.apply(chr);
			chr.evaluate(prob);
			if(this->m_repair) {
				this->m_repair->repair(chr,prob);
			}
			this->chromosome_evaluated(chr);
			if(comp->compare(chr,saved) == -1) {
				saved = chr;
				at_local_opt = false;
				break;
			} else {
				chr=saved;
			}
		}
	}
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
steepest_descent<Chromosome,Encoding>::steepest_descent()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
steepest_descent<Chromosome,Encoding>::~steepest_descent()
{
}

/*!
 * \brief improve a chromosome using steepest descent local search
 */
template <template <typename> class Chromosome, typename Encoding>
void steepest_descent<Chromosome,Encoding>::improve(Chromosome<Encoding>& chr,
        comparator<Chromosome,Encoding>* comp,
        const typename Encoding::ProblemType* prob)
{
	this->reset();

	Chromosome<Encoding> current_best = chr;
	bool at_local_opt = false;
	while(!at_local_opt && !this->terminate()) {
		at_local_opt = true;
		this->m_nf->initialize(chr);
		while(this->m_nf->has_more_neighbors()) {
			Chromosome<Encoding> next(chr);
			move<Chromosome,Encoding> m=this->m_nf->next_neighbor();
			m.apply(next);
			next.evaluate(prob);
			if(this->m_repair) {
				this->m_repair->repair(next,prob);
			}
			this->chromosome_evaluated(next);
			if(comp->compare(next,current_best) == -1) {
				current_best = next;
				at_local_opt = false;
			}
		}
		chr = current_best;
	}
}

/*!
 * \brief create and return an initialized local search operator
 */
template <template <typename> class Chromosome, typename Encoding>
local_search<Chromosome,Encoding>*
local_search_factory<Chromosome,Encoding>::construct()
{
	string hcname;
	kvparse::parameter_value(this->m_prefix+keywords::LOCAL_SEARCH, hcname, true);
	if(hcname == keywords::NEXT_DESCENT) {
		local_search<Chromosome,Encoding>* hc = new next_descent<Chromosome,Encoding>;
		hc->set_prefix(this->m_prefix);
		hc->initialize();
		return hc;
	} else if(hcname == keywords::STEEPEST_DESCENT) {
		local_search<Chromosome,Encoding>* hc = new steepest_descent<Chromosome,Encoding>;
		hc->set_prefix(this->m_prefix);
		hc->initialize();
		return hc;
	} else if(hcname == keywords::TABU_SEARCH) {
		local_search<Chromosome,Encoding>* hc = new tabu_search<Chromosome,Encoding>;
		hc->set_prefix(this->m_prefix);
		hc->initialize();
		return hc;
	} else if(hcname == keywords::SIMULATED_ANNEALING) {
		local_search<Chromosome,Encoding>* hc = new simulated_annealing<Chromosome,Encoding>;
		hc->set_prefix(this->m_prefix);
		hc->initialize();
		return hc;
	} else if(hcname==keywords::VARIABLE_DEPTH_SEARCH) {
//         vds<Chromosome,Encoding>* hc=new vds<Chromosome,Encoding>;
//         hc->set_prefix(this->m_prefix);
//         hc->initialize();
//         return hc;
		return 0;
	} else {
		cerr << "illegal value for parameter " << this->m_prefix+keywords::LOCAL_SEARCH << ": "
		     << hcname << " specified" << endl;
		exit(0);
		return 0;
	}
}

