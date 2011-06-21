/*!
 * \file population.cpp
 *
 * evolutionary algorithms utilize a population of candidate solutions.
 * this class maintains this population for use by the evolutionary
 * algorithms.
 *
 * Deon Garrett
 * deong@acm.org
 */

#include <deque>
#include <algorithm>
#include <functional>
#include <iostream>
#include <fstream>
#include <cstring>
#include "chromosome.h"
#include "encoding.h"
#include "population.h"
#include "comparator.h"
#include "problems.h"

using namespace std;

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
population<Chromosome,Encoding>::population() :
    m_count(0)
{
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
population<Chromosome,Encoding>::population(unsigned int n) :
    m_individuals(n),
    m_count(0)
{
}

/*!
 * \brief copy constructor
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
population<Chromosome,Encoding>::population(const population& that)
{
    m_individuals = that.m_individuals;
    m_count = that.m_count;
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
population<Chromosome,Encoding>::~population()
{
}

/*!
 * \brief assignment operator
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
population<Chromosome,Encoding>& population<Chromosome,Encoding>::operator=(const population& that)
{
    m_individuals = that.m_individuals;
    m_count = that.m_count;
    return *this;
}

/*!
 * \brief equality operator
 *
 * \bug only works if both populations contain the same chromosomes in
 * the same order
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
bool population<Chromosome,Encoding>::operator==(const population& that) const
{
    if(m_count != that.m_count)
    {
        return false;
    }
    for(unsigned int i=0; i<m_count; i++)
    {
        if(m_individuals[i] != that.m_individuals[i])
        {
            return false;
        }
    }
    return true;
}

/*!
 * \brief inequality operator
 *
 * \bug incorrect if the two populations have the same chromosomes
 * in different orders
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
bool population<Chromosome,Encoding>::operator!=(const population& that) const
{
    if(m_count != that.m_count)
    {
        return true;
    }
    for(unsigned int i=0; i<m_count; i++)
    {
        if(m_individuals[i] != that.m_individuals[i])
        {
            return true;
        }
    }
    return false;
}
    
/*!
 * \brief return the number of chromosomes in the population
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
inline unsigned int population<Chromosome,Encoding>::size() const
{
    return m_count;
}

/*!
 * \brief construct the union of two populations (with duplicates)
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
population<Chromosome,Encoding>
population<Chromosome,Encoding>::operator+(const population<Chromosome,Encoding>& that) const
{
    population<Chromosome,Encoding> msunion = *this;
    for(unsigned int i=0; i<that.m_count; i++)
    {
        msunion.add(that.m_individuals[i]);
    }
    return msunion;
}

/*!
 * \brief construct set difference between two populations (with duplicates)
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
population<Chromosome,Encoding>
population<Chromosome,Encoding>::operator-(const population<Chromosome,Encoding>& that) const
{
    population<Chromosome,Encoding> diff = *this;
    for(unsigned int i=0; i<that.m_count; i++)
    {
        diff.remove(that.m_individuals[i]);
    }
    return diff;
}

/*!
 * \brief add in elements of a given population to the current population
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
population<Chromosome,Encoding>
population<Chromosome,Encoding>::operator+=(const population<Chromosome,Encoding>& that)
{
    *this = *this + that;
    return *this;
}

/*!
 * \brief subtract elements of a given population from the current population
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
population<Chromosome,Encoding>
population<Chromosome,Encoding>::operator-=(const population<Chromosome,Encoding>& that)
{
    *this = *this - that;
    return *this;
}

/*!
 * \brief retrieve chromosomes by index
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
inline Chromosome<Encoding>& population<Chromosome,Encoding>::operator[](unsigned int i)
{
    return m_individuals[i];
}

/*!
 * \brief retrieve chromosomes by index
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
inline const Chromosome<Encoding>& population<Chromosome,Encoding>::operator[](unsigned int i) const
{
    return m_individuals[i];
}

/*!
 * \brief resize the population
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void population<Chromosome,Encoding>::resize(unsigned int sz)
{
    m_individuals.resize(sz);
}

/*!
 * \brief add an individual to the population
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
bool population<Chromosome,Encoding>::add(const Chromosome<Encoding>& ind)
{
    if(m_count < m_individuals.size())
    {
        m_individuals[m_count++] = ind;
    }
    else
    {
        m_individuals.push_back(ind);
        m_count++;
    }
    return true;
}

/*!
 * \brief remove an individual from the population
 *
 * if more than one copy is present, the first one is removed
 * 
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void population<Chromosome,Encoding>::remove(const Chromosome<Encoding>& ind)
{
    for(unsigned int i=0; i<m_count; i++)
    {
        if(m_individuals[i] == ind)
        {
            m_individuals.erase(m_individuals.begin()+i);
            m_count--;
            return;
        }
    }
}

/*!
 * \brief remove an individual by index
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void population<Chromosome,Encoding>::remove_at(unsigned int index)
{
    m_individuals.erase(m_individuals.begin()+index);
    m_count--;
}

/*!
 * \brief clear all individuals from the population
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void population<Chromosome,Encoding>::clear()
{
    m_individuals.clear();
    m_count = 0;
}

/*!
 * \brief return an iterator to the beginning of the population
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
typename population<Chromosome,Encoding>::iterator population<Chromosome,Encoding>::begin()
{
    return m_individuals.begin();
}

/*!
 * \brief return an iterator to the end of the population
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
typename population<Chromosome,Encoding>::iterator population<Chromosome,Encoding>::end()
{
    return m_individuals.end();
}

/*!
 * \brief return a const_iterator to the beginning of the population
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
typename population<Chromosome,Encoding>::const_iterator population<Chromosome,Encoding>::begin() const
{
    return m_individuals.begin();
}

/*!
 * \brief return a const_iterator to the end of the population
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
typename population<Chromosome,Encoding>::const_iterator population<Chromosome,Encoding>::end() const
{
    return m_individuals.end();
}

/*!
 * \brief sort the population using the specified comparator
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void population<Chromosome,Encoding>::sort(comparator<Chromosome,Encoding>* comp)
{
    std::sort(m_individuals.begin(), m_individuals.begin()+m_count,
              sorting_comparator<Chromosome,Encoding>(comp));
}

/*!
 * \brief print the individuals in the population
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
ostream& operator<<(ostream& s, const population<Chromosome,Encoding>& p)
{
    s << "Population" << endl;
    for(unsigned int i=0; i<p.size(); i++)
    {
        s << p[i] << endl;
    }
    return s;
}

/*!
 * \brief read the population from a file
 *
 * \author deong@acm.org
 * \date 10/13/2008
 */
template <template <typename> class Chromosome, typename Encoding>
istream& operator>>(istream& istr, population<Chromosome,Encoding>& p)
{
    typename Encoding::ProblemFactoryType pf;
    typename Encoding::ProblemType* prob;
    prob=pf.construct();
    
    while(true)
    {
        Chromosome<Encoding> sol(prob);
        istr >> sol;
        if(istr.eof())
            break;
        p.add(sol);
    }
    return istr;
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
sorting_comparator<Chromosome,Encoding>::sorting_comparator(comparator<Chromosome,Encoding>* comp) :
    m_comp(comp)
{
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
sorting_comparator<Chromosome,Encoding>::~sorting_comparator()
{
}

/*!
 * \brief compare chromosomes using the specified comparator
 *
 * \author deong
 * \date 05/10/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
bool sorting_comparator<Chromosome,Encoding>::operator()(const Chromosome<Encoding>& c1,
                                                         const Chromosome<Encoding>& c2) const
{
    return m_comp->compare(c1,c2) == -1;
}
