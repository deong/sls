/*!
 * \file replacement.cpp
 *
 * In evolutionary algorithms, the population at time t+1 is a function
 * of the population at time t and the population of offspring produced
 * at time t.  There are, however, different measures of time.  In a
 * generational GA, the population and offspring are complete populations
 * which are combined somehow into the new population.  In a steady state
 * algorithm, offspring are produced and inserted into the new population
 * one at a time.  These classes provide both behaviours.
 *
 * Deon Garrett
 * deong@acm.org
 */

#include <iostream>
#include "replacement.h"
#include "chromosome.h"
#include "encoding.h"
#include "population.h"
#include "comparator.h"
#include "configuration.h"
#include "keywords.h"
#include "utilities.h"

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
replacement_scheme<Chromosome,Encoding>::replacement_scheme()
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
replacement_scheme<Chromosome,Encoding>::~replacement_scheme()
{
    delete m_comp;
}

/*!
 * \brief initialize the replacement scheme
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
void replacement_scheme<Chromosome,Encoding>::initialize()
{
    comparator_factory<Chromosome,Encoding> cf;
    m_comp = cf.construct();
}

/*!
 * \brief construct next population from parent and offspring populations
 *
 * empty method defined so subclasses can choose whether to operate on
 * a population basis or an individual chromosome basis
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
void replacement_scheme<Chromosome,Encoding>::merge_populations(population<Chromosome,Encoding>& parents,
                                                                population<Chromosome,Encoding>& offspring,
                                                                population<Chromosome,Encoding>& next) const
{
}

/*!
 * \brief merge a single individual in to form a new population
 *
 * empty method defined so subclasses can choose whether to operate on
 * a population basis or an individual chromosome basis
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
void replacement_scheme<Chromosome,Encoding>::merge_individual(Chromosome<Encoding>& sol,
                                                               population<Chromosome,Encoding>& parents,
                                                               population<Chromosome,Encoding>& next) const
{
}

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
generational_replacement<Chromosome,Encoding>::generational_replacement()
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
generational_replacement<Chromosome,Encoding>::~generational_replacement()
{
}

/*!
 * \brief overwrite old parent population with offspring
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
void generational_replacement<Chromosome,Encoding>::merge_populations(population<Chromosome,Encoding>& parents,
                                                                      population<Chromosome,Encoding>& offspring,
                                                                      population<Chromosome,Encoding>& next) const
{
    next = offspring;
}

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
elitist_replacement<Chromosome,Encoding>::elitist_replacement()
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
elitist_replacement<Chromosome,Encoding>::~elitist_replacement()
{
}

/*!
 * \brief initialize the number of elites to keep each generation
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
void elitist_replacement<Chromosome,Encoding>::initialize()
{
    replacement_scheme<Chromosome,Encoding>::initialize();
    m_elites = 1;
    configuration::unsigned_integer_parameter(keywords::ELITES, m_elites, false);
}

/*!
 * \brief overwrite parents with offspring but keep best N parents
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
void elitist_replacement<Chromosome,Encoding>::merge_populations(population<Chromosome,Encoding>& parents,
                                                                 population<Chromosome,Encoding>& offspring,
                                                                 population<Chromosome,Encoding>& next) const
{
    next = offspring;
    
    // find the best individual in the parent population
    int best_idx = 0;
    Chromosome<Encoding> best = parents[0];
    for(unsigned int i=1; i<parents.size(); i++)
    {
        if(this->m_comp->compare(parents[i],best) == -1)
        {
            best = parents[i];
            best_idx = i;
        }
    }

    // find the worst individual in the offspring population
    int worst_idx = 0;
    Chromosome<Encoding> worst = offspring[0];
    for(unsigned int i=1; i<offspring.size(); i++)
    {
        if(this->m_comp->compare(worst,offspring[i]) == -1)
        {
            worst = offspring[i];
            worst_idx = i;
        }
    }

    // put the best parent into the next population in place
    // of the worst offspring
    next[worst_idx] = best;
}

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
truncation_replacement<Chromosome,Encoding>::truncation_replacement()
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
truncation_replacement<Chromosome,Encoding>::~truncation_replacement()
{
}

/*!
 * \brief take best N of parents+offspring
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
void truncation_replacement<Chromosome,Encoding>::merge_populations(population<Chromosome,Encoding>& parents,
                                                                    population<Chromosome,Encoding>& offspring,
                                                                    population<Chromosome,Encoding>& next) const
{
    /*
     * precondition: parent population is sorted 
     */
    
    // sort the offspring
    offspring.sort(this->m_comp);

    // store the next generation temporarily 
    population<Chromosome,Encoding> best(parents.size());
    
    unsigned int p = 0;         // parent index
    unsigned int o = 0;         // offspring index
    for(unsigned int i=0; i<parents.size(); i++)
    {
        if(o >= offspring.size() || this->m_comp->compare(parents[p],offspring[o]) == -1)
        {
            best.add(parents[p++]);
        }
        else
        {
            best.add(offspring[o++]);
        }
    }
    next = best;
}

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
steady_state_replacement<Chromosome,Encoding>::steady_state_replacement()
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
steady_state_replacement<Chromosome,Encoding>::~steady_state_replacement()
{
}

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
replace_worst<Chromosome,Encoding>::replace_worst()
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
replace_worst<Chromosome,Encoding>::~replace_worst()
{
}

/*!
 * \brief replace worst individual in the population
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
void replace_worst<Chromosome,Encoding>::merge_individual(Chromosome<Encoding>& sol,
                                                          population<Chromosome,Encoding>& pop,
                                                          population<Chromosome,Encoding>& next) const
{
    unsigned int i=0;
    while(i<pop.size() && this->m_comp->compare(pop[i],sol) == -1)
    {
        next[i] = pop[i];
        i++;
    }
    if(i<pop.size()-1)
    {
        next[i++] = sol;
    }
    while(i<pop.size()-1)
    {
        next[i] = pop[i-1];
        i++;
    }
}

/*!
 * \brief create and initialize a replacement scheme
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
replacement_scheme<Chromosome,Encoding>* replacement_scheme_factory<Chromosome,Encoding>::construct()
{
    string rsname;
    configuration::string_parameter(keywords::REPLACEMENT_SCHEME, rsname, true);

    if(rsname == keywords::GENERATIONAL_REPLACEMENT)
    {
        generational_replacement<Chromosome,Encoding>* gr = new  generational_replacement<Chromosome,Encoding>;
        gr->initialize();
        return gr;
    }
    else if(rsname == keywords::ELITIST_REPLACEMENT)
    {
        elitist_replacement<Chromosome,Encoding>* er = new elitist_replacement<Chromosome,Encoding>;
        er->initialize();
        return er;
    }
    else if(rsname == keywords::TRUNCATION_REPLACEMENT)
    {
        truncation_replacement<Chromosome,Encoding>* tr = new truncation_replacement<Chromosome,Encoding>;
        tr->initialize();
        return tr;
    }
    else if(rsname == keywords::REPLACE_WORST_REPLACEMENT)
    {
        replace_worst<Chromosome,Encoding>* rwr = new replace_worst<Chromosome,Encoding>;
        rwr->initialize();
        return rwr;
    }
    else
    {
        error("invalid replacement_scheme specified: " + rsname);
        return 0;
    }
}
