/*!
 * \file simanneal.cpp
 *
 * Simulated Annealing
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <cmath>
#include "simanneal.h"
#include "localsearch.h"
#include "chromosome.h"
#include "encoding.h"
#include "problems.h"
#include "comparator.h"
#include "move.h"
#include "neighborhood.h"
#include "metrics.h"
#include "terminators.h"
#include "configuration.h"
#include "keywords.h"

using namespace std;

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
simulated_annealing<Chromosome,Encoding>::simulated_annealing()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
simulated_annealing<Chromosome,Encoding>::~simulated_annealing()
{
}

/*!
 * \brief initialize the algorithm components
 */
template <template <typename> class Chromosome, typename Encoding>
void simulated_annealing<Chromosome,Encoding>::initialize()
{
    local_search<Chromosome,Encoding>::initialize();
    double init_accept;
    double final_accept;
    configuration::double_parameter(keywords::INITIAL_50_PERCENTILE,init_accept,true);
    configuration::double_parameter(keywords::FINAL_50_PERCENTILE,final_accept,true);
    m_temp = -init_accept / log(0.5);

    // compute the appropriate alpha
    unsigned int maxeval;
    unsigned int steps;
    configuration::unsigned_integer_parameter("ls_"+keywords::MAX_EVALUATIONS,maxeval,true);
    configuration::unsigned_integer_parameter(keywords::COOLING_STEPS,steps,true);
    m_alpha=pow((final_accept/init_accept),1.0/steps);
    m_schedule=maxeval/steps;
}

/*!
 * \brief reset the algorithm parameters
 */
template <template <typename> class Chromosome, typename Encoding>
void simulated_annealing<Chromosome,Encoding>::reset()
{
    local_search<Chromosome,Encoding>::reset();
    double init_accept;
    double final_accept;
    configuration::double_parameter(keywords::INITIAL_50_PERCENTILE,init_accept,true);
    configuration::double_parameter(keywords::FINAL_50_PERCENTILE,final_accept,true);
    m_temp = -init_accept / log(0.5);

    // compute the appropriate alpha
    unsigned int maxeval;
    unsigned int steps;
    configuration::unsigned_integer_parameter("ls_"+keywords::MAX_EVALUATIONS,maxeval,true);
    configuration::unsigned_integer_parameter(keywords::COOLING_STEPS,steps,true);
    m_alpha=pow((final_accept/init_accept),1.0/steps);
    m_schedule=maxeval/steps;
}

/*!
 * \brief determine whether to accept a given move
 */
template <template <typename> class Chromosome, typename Encoding>
bool simulated_annealing<Chromosome,Encoding>::accept_move(double delta) const
{
    if(delta<0)
    {
	return true;
    }
    else
    {
	mtrandom mt;
	if(mt.random()<exp(-delta/m_temp))
	{
	    return true;
	}
	else
	{
	    return false;
	}
    }
}

/*!
 * \brief improve a chromosome using the simulated annealing algorithm
 */
template <template <typename> class Chromosome, typename Encoding>
void simulated_annealing<Chromosome,Encoding>::improve(Chromosome<Encoding>& chr,
						       comparator<Chromosome,Encoding>* comp,
						       const typename Encoding::ProblemType* prob)
{
    this->reset();
    
    // make sure the comparator is a scalarizing comparator
    scalarizing_comparator<Chromosome,Encoding>* sc;
    sc=dynamic_cast<scalarizing_comparator<Chromosome,Encoding>*>(comp);
    
    if(!sc)
    {
	error("attempt to use simulated annealing with comparator other than scalarizing");
	return;
    }

    mtrandom mt;
    unsigned int evals=0;

    //! keep track of current best known solution
    Chromosome<Encoding> best=chr;
    
    // loop until we reach maximum iterations
    while(!this->terminate())
    {
	this->m_nf->initialize(chr);
	while(this->m_nf->has_more_neighbors())
	{
	    // generate the next move
	    move<Chromosome,Encoding> m=this->m_nf->next_neighbor();
	    Chromosome<Encoding> tmp=chr;

	    // evaluate the move 
	    m.apply(tmp);
	    tmp.evaluate(prob);
	    this->chromosome_evaluated(tmp);

	    bool accepted=false;
	    double delta=sc->difference(tmp,chr);
	    if(this->accept_move(delta))
	    {
		chr=tmp;
		if(sc->compare(chr,best)<0)
		    best=chr;
		accepted=true;
	    }

	    // exponentially decay the temperature
	    if(++evals%m_schedule == 0)
		m_temp*=m_alpha;

	    if(accepted)
		break;
	}
    }

    // restore the best found
    chr=best;
}
