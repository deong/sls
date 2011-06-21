/*!
 * \file selection.h
 *
 * evolutionary algorithms select individuals for reproduction based
 * upon fitness.  these classes provide commonly used mechanisms for
 * performing this selection.
 *
 * Deon Garrett
 * deong@acm.org
 */

#ifndef _SELECTION_H_
#define _SELECTION_H_

#include "chromosome.h"
#include "population.h"

using namespace std;

/*!
 * \class selection_scheme
 *
 * \author deong
 * \date 05/11/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class selection_scheme
{
private:
    // disable copying of functional classes
    selection_scheme(const selection_scheme& that);
    selection_scheme& operator=(const selection_scheme& that);
    
protected:
    population<Chromosome,Encoding>* m_population;

public:
    selection_scheme();
    virtual ~selection_scheme();
    virtual void set_population(const population<Chromosome,Encoding>* pop);
    virtual Chromosome<Encoding>& select_parent() const = 0;
};

/*!
 * \class tournament_selection
 *
 * \author deong
 * \date 05/11/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class tournament_selection : public selection_scheme<Chromosome,Encoding>
{
private:
    // disable copying of functional classes
    tournament_selection(const tournament_selection& that);
    tournament_selection& operator=(const tournament_selection& that);

    // comparison method for chromosomes
    comparator<Chromosome,Encoding>* m_comp;

    // determine whether to delete the comparator in the destructor
    bool m_delete_comp;
    
public:
    tournament_selection();
    tournament_selection(comparator<Chromosome,Encoding>* comp);
    virtual ~tournament_selection();

    virtual void initialize();
    virtual Chromosome<Encoding>& select_parent() const;
};

/*!
 * \class ranking_selection
 *
 * \author deong
 * \date 05/11/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class ranking_selection : public selection_scheme<Chromosome,Encoding>
{
private:
    // disable copying of functional classes
    ranking_selection(const ranking_selection& that);
    ranking_selection& operator=(const ranking_selection& that);

protected:
    double m_bias;
    
public:
    ranking_selection();
    virtual ~ranking_selection();

    virtual void initialize();
    virtual Chromosome<Encoding>& select_parent() const;
};

/*!
 * \class random_selection
 *
 * \author deong
 * \date 05/11/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class random_selection : public selection_scheme<Chromosome,Encoding>
{
private:
    // disable copying of functional classes
    random_selection(const random_selection& that);
    random_selection& operator=(const random_selection& that);

public:
    random_selection();
    virtual ~random_selection();

    virtual Chromosome<Encoding>& select_parent() const;
};

/*!
 * \class selection_scheme_factory
 *
 * \author deong
 * \date 05/11/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class selection_scheme_factory
{
public:
    static selection_scheme<Chromosome,Encoding>* construct();
};

#include "selection.cpp"

#endif
    
