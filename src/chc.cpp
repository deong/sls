/*!
 * \file chc.cpp
 *
 * Eshelman's CHC algorithm
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <iostream>
#include "chc.h"
#include "chromosome.h"
#include "encoding.h"
#include "population.h"
#include "crossover.h"
#include "mutation.h"
#include "comparator.h"
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
chc<Chromosome,Encoding>::chc() :
    m_sel_scheme(0),
    m_rep_scheme(0),
    m_cross_op(0),
    m_distmeas(0),
    m_difference_threshold(0),
    m_divergence_rate(0.35)
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
chc<Chromosome,Encoding>::~chc()
{
    delete m_sel_scheme;
    delete m_rep_scheme;
    delete m_cross_op;
    delete m_comp;
    delete m_distmeas;
}

/*!
 * \brief initialize the algorithm
 */
template <template <typename> class Chromosome, typename Encoding>
void chc<Chromosome,Encoding>::initialize()
{
    // initialize the base class members
    ea<Chromosome,Encoding>::initialize();

    // initialize the comparator
    comparator_factory<Chromosome,Encoding> cf;
    m_comp = cf.construct();
	
    // sort the population by fitness
    this->m_population.sort(m_comp);
    
    // configure the selection and replacement schemes
    m_sel_scheme = new random_selection<Chromosome,Encoding>;
    m_rep_scheme = new truncation_replacement<Chromosome,Encoding>;
	m_rep_scheme->initialize();
	
    // figure out which crossover operator to use
    m_cross_op = new hux_crossover<Chromosome,Encoding>;

    // configure the neighborhood operator to use to compute distance
    m_distmeas = neighborhood_factory<Chromosome,Encoding>::construct();
    
    // set up the difference threshold to l/4
    m_difference_threshold = this->m_population[0].length()/4;
}

/*!
 * \brief run the CHC search algorithm
 */
template <template <typename> class Chromosome, typename Encoding>
void chc<Chromosome,Encoding>::run()
{
    while(!this->terminate())
    {
        this->run_one_generation();
        this->generation_completed(this->m_population);
    }
    this->compute_metrics();
    this->report_metrics(cout);
}

/*!
 * \brief run a single generation of the algorithm
 */
template <template <typename> class Chromosome, typename Encoding>
void chc<Chromosome,Encoding>::run_one_generation()
{
    mtrandom mt;
    population<Chromosome,Encoding> offspring;
    
    m_sel_scheme->set_population(&this->m_population);

    for(unsigned int i=0; i<this->m_population.size()/2; i++)
    {
        // select two parents
        Chromosome<Encoding> p1 = m_sel_scheme->select_parent();
        Chromosome<Encoding> p2 = m_sel_scheme->select_parent();

        // perform crossover unless incestuous
        if(m_distmeas->distance_between(p1, p2) > m_difference_threshold)
        {
            Chromosome<Encoding> c1 = p1;
            Chromosome<Encoding> c2 = p2;

            m_cross_op->crossover(p1, p2, c1, c2);

            // evaluate the offspring
			if(this->m_repair)
				this->m_repair->repair(c1,this->m_fitfunc);
			c1.evaluate(this->m_fitfunc);
			this->chromosome_evaluated(c1);
			
			if(this->m_repair)
				this->m_repair->repair(c2,this->m_fitfunc);
			c2.evaluate(this->m_fitfunc);
			this->chromosome_evaluated(c2);
			
			offspring.add(c1);
			offspring.add(c2);
        }
    }

    // construct the next population -- if it is the same as the
    // previous population, then decrement the difference threshold
    population<Chromosome,Encoding> next;
    m_rep_scheme->merge_populations(this->m_population, offspring, next);
    if(next == this->m_population)
    {
        m_difference_threshold--;
    }
    else
    {
        this->m_population = next;
    }
    
    if(m_difference_threshold < 0)
    {
        diverge();
        m_difference_threshold = (int)(m_divergence_rate * (1.0 - m_divergence_rate) *
                                       this->m_population[0].length());
    }
}

/*!
 * \brief perform a "cataclysmic mutation" to restart the search
 */
template <template <typename> class Chromosome, typename Encoding>
void chc<Chromosome,Encoding>::diverge()
{
    mtrandom mt;
    unsigned int n = this->m_population[0].length();
    int flips = (int)(m_divergence_rate * n);
    for(unsigned int i=1; i<this->m_population.size(); i++)
    {
        this->m_population[i] = this->m_population[0];

        vector<int> order = mt.permutation(n);
        for(int j=0; j<flips; j++)
        {
            this->m_population[i][order[j]] = (this->m_population[i][order[j]] == 0) ? 1 : 0;
        }
		if(this->m_repair)
			this->m_repair->repair(this->m_population[i],this->m_fitfunc);
		this->m_population[i].evaluate(this->m_fitfunc);
		this->chromosome_evaluated(this->m_population[i]);
    }
}
