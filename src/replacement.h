/*!
 * \file replacement.h
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

#ifndef _REPLACEMENT_H_
#define _REPLACEMENT_H_

#include "chromosome.h"
#include "population.h"
#include "comparator.h"

using namespace std;

/*!
 * \class replacement_scheme
 *
 * \author deong
 * \date 05/11/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class replacement_scheme
{
private:
    // disable copying of functional classes
    replacement_scheme(const replacement_scheme& that);
    replacement_scheme& operator=(const replacement_scheme& that);

protected:
    comparator<Chromosome,Encoding>* m_comp;
    
public:
    replacement_scheme();
    virtual ~replacement_scheme();

    virtual void initialize();
    virtual void merge_populations(population<Chromosome,Encoding>& parents,
                                   population<Chromosome,Encoding>& offspring,
                                   population<Chromosome,Encoding>& next) const;
    virtual void merge_individual(Chromosome<Encoding>& sol,
                                  population<Chromosome,Encoding>& pop,
                                  population<Chromosome,Encoding>& next) const;
};

/*!
 * \class generational_replacement
 *
 * replace entire populations with their offspring
 *
 * \author deong
 * \date 05/11/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class generational_replacement : public replacement_scheme<Chromosome,Encoding>
{
private:
    // disable copying of functional class
    generational_replacement(const generational_replacement& that);
    generational_replacement& operator=(const generational_replacement& that);

public:
    generational_replacement();
    virtual ~generational_replacement();

    virtual void merge_populations(population<Chromosome,Encoding>& parents,
                                   population<Chromosome,Encoding>& offspring,
                                   population<Chromosome,Encoding>& next) const;
};

/*!
 * \class elitist_replacement
 *
 * replace entire population with offspring, but keep a small number of the best
 * parents in favor of the worst offspring.
 *
 * \author deong
 * \date 05/11/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class elitist_replacement : public generational_replacement<Chromosome,Encoding>
{
private:
    // disable copying of functional class
    elitist_replacement(const elitist_replacement& that);
    elitist_replacement& operator=(const elitist_replacement& that);

protected:
    unsigned int m_elites;
    
public:
    elitist_replacement();
    virtual ~elitist_replacement();

    virtual void initialize();
    virtual void merge_populations(population<Chromosome,Encoding>& parents,
                                   population<Chromosome,Encoding>& offspring,
                                   population<Chromosome,Encoding>& next) const;
};

/*!
 * \class truncation_replacement
 *
 * very strongly elitist, takes the best N of the union of parents and
 * offspring
 *
 * \author deong
 * \date 05/11/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class truncation_replacement : public replacement_scheme<Chromosome,Encoding>
{
private:
    // disable copying
    truncation_replacement(const truncation_replacement& that);
    truncation_replacement& operator=(const truncation_replacement& that);

public:
    truncation_replacement();
    virtual ~truncation_replacement();

    virtual void merge_populations(population<Chromosome,Encoding>& parents,
                                   population<Chromosome,Encoding>& offspring,
                                   population<Chromosome,Encoding>& next) const;
};

/*!
 * \class steady_state_replacement
 *
 * replace one chromosome at a time
 *
 * \author deong
 * \date 05/11/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class steady_state_replacement : public replacement_scheme<Chromosome,Encoding>
{
private:
    // disable copying of functional class
    steady_state_replacement(const steady_state_replacement& that);
    steady_state_replacement& operator=(const steady_state_replacement& that);

public:
    steady_state_replacement();
    virtual ~steady_state_replacement();

    virtual void merge_individual(Chromosome<Encoding>& sol,
                                  population<Chromosome,Encoding>& pop,
                                  population<Chromosome,Encoding>& next) const = 0;
};

/*!
 * \class replace_worst
 *
 * steady state method in which the worst population member is replaced
 *
 * \author deong
 * \date 05/11/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class replace_worst : public steady_state_replacement<Chromosome,Encoding>
{
private:
    // disable copying functional class
    replace_worst(const replace_worst& that);
    replace_worst& operator=(const replace_worst& that);

public:
    replace_worst();
    virtual ~replace_worst();

    virtual void merge_individual(Chromosome<Encoding>& sol,
                                  population<Chromosome,Encoding>& pop,
                                  population<Chromosome,Encoding>& next) const;
};

/*!
 * \class replacement_scheme_factory
 *
 * \author deong
 * \date 05/11/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class replacement_scheme_factory
{
public:
    static replacement_scheme<Chromosome,Encoding>* construct();
};

#include "replacement.cpp"

#endif
