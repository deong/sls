/*!
 * \file genitor.cpp
 *
 * Whitley's steady state ga which first advocated the use of
 * ranking for parent selection
 *
 * Deon Garrett
 * deong@acm.org
 */

#include <iostream>
#include "genitor.h"
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
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
genitor<Chromosome,Encoding>::genitor() :
    m_sel_scheme(0),
    m_rep_scheme(0),
    m_cross_op(0),
    m_cross_rate(1.0),
    m_mut_op(0)
{
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
genitor<Chromosome,Encoding>::~genitor()
{
    delete m_sel_scheme;
    delete m_rep_scheme;
    delete m_cross_op;
    delete m_mut_op;
    delete m_comp;
}

/*!
 * \brief initialize the algorithm components
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void genitor<Chromosome,Encoding>::initialize()
{
    // initialize the base class members
    ea<Chromosome,Encoding>::initialize();

    // initialize the comparator
    comparator_factory<Chromosome,Encoding> cf;
    m_comp = cf.construct();
    
    // sort the population by fitness
    this->m_population.sort(m_comp);
    
    // configure the selection and replacement schemes
    m_sel_scheme = new ranking_selection<Chromosome,Encoding>;
    m_sel_scheme->initialize();
    
    m_rep_scheme = new replace_worst<Chromosome,Encoding>;
    m_rep_scheme->initialize();
    
    // figure out which crossover operator to use
    m_cross_op = crossover_operator_factory<Chromosome,Encoding>::construct();

    // set up the crossover rate
    configuration::double_parameter(keywords::CROSSOVER_RATE, m_cross_rate, false);
    
    // figure out which mutation operator to use
    m_mut_op = mutation_operator_factory<Chromosome,Encoding>::construct();
}

/*!
 * \brief run the genitor algorithm
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void genitor<Chromosome,Encoding>::run()
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
 * \brief run one generation of the genitor algorithm
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void genitor<Chromosome,Encoding>::run_one_generation()
{
    mtrandom mt;

    m_sel_scheme->set_population(&this->m_population);
    
    // select two parents
    Chromosome<Encoding> p1 = m_sel_scheme->select_parent();
    Chromosome<Encoding> p2 = m_sel_scheme->select_parent();

    // perform crossover with some probability
    Chromosome<Encoding> c1 = p1;
    Chromosome<Encoding> c2 = p2;
    if(mt.random() < m_cross_rate)
    {
        m_cross_op->crossover(p1, p2, c1, c2);
    }
    Chromosome<Encoding> offspring;
    if(mt.random() < 0.5)
    {
        offspring = c1;
    }
    else
    {
        offspring = c2;
    }
    
    // mutate the offspring
    m_mut_op->mutate(offspring);

    // evaluate the offspring
    offspring.evaluate(this->m_fitfunc);
    chromosome_evaluated(offspring);
    
    // put the offspring into the population in place of the
    // worst individual
    m_rep_scheme->merge_individual(offspring, this->m_population, this->m_population);
}
