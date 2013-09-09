/*!
 * \file move.cpp
 *
 * attempt to abstract away the notion of a local search move
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <utility>
#include <deque>
#include <algorithm>
#include "move.h"
#include "chromosome.h"
#include "encoding.h"

using namespace std;

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
move<Chromosome,Encoding>::move()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
move<Chromosome,Encoding>::~move()
{
}

/*!
 * \brief compare two moves for equality
 *
 * if method=ANY, then moves are equal if they put any single allele value
 * in the same state.  if method=ALL, then all alleles must be assigned the
 * same values for the moves to be equivalent.
 */
template <template <typename> class Chromosome, typename Encoding>
bool move<Chromosome,Encoding>::operator==(const move<Chromosome,Encoding>& that) const
{
    if(m_components.size()!=that.m_components.size())
    {
	return false;
    }

    for(unsigned int i=0; i<m_components.size(); i++)
    {
	if(m_components[i].first!=that.m_components[i].first ||
	   m_components[i].second!=that.m_components[i].second)
	{
	    return false;
	}
    }
    return true;
}

/*!
 * \brief return a specified move component
 */
template <template <typename> class Chromosome, typename Encoding>
pair<unsigned int,typename Encoding::Genotype>& move<Chromosome,Encoding>::operator[](unsigned int i)
{
    return m_components[i];
}

/*!
 * \brief return a specified move component
 */
template <template <typename> class Chromosome, typename Encoding>
const pair<unsigned int,typename Encoding::Genotype>& move<Chromosome,Encoding>::operator[](unsigned int i) const
{
    return m_components[i];
}

/*!
 * \brief return the number of components in the move
 */
template <template <typename> class Chromosome, typename Encoding>
unsigned int move<Chromosome,Encoding>::length() const
{
    return m_components.size();
}

/*!
 * \brief return an iterator to the front of the move
 */
template <template <typename> class Chromosome, typename Encoding>
typename move<Chromosome,Encoding>::iterator move<Chromosome,Encoding>::begin()
{
    return m_components.begin();
}

/*!
 * \brief return an iterator to the end of the move
 */
template <template <typename> class Chromosome, typename Encoding>
typename move<Chromosome,Encoding>::iterator move<Chromosome,Encoding>::end()
{
    return m_components.end();
}

/*!
 * \brief return a const_iterator to the front of the move
 */
template <template <typename> class Chromosome, typename Encoding>
typename move<Chromosome,Encoding>::const_iterator move<Chromosome,Encoding>::begin() const
{
    return m_components.begin();
}

/*!
 * \brief return a const_iterator to the end of the list
 */
template <template <typename> class Chromosome, typename Encoding>
typename move<Chromosome,Encoding>::const_iterator move<Chromosome,Encoding>::end() const
{
    return m_components.end();
}

/*!
 * \brief add a new move component
 *
 * adds the component in sorted order by index, overwriting any previous
 * assignment of the given position
 */
template <template <typename> class Chromosome, typename Encoding>
void move<Chromosome,Encoding>::add_component(unsigned int pos, typename Encoding::Genotype val)
{
    unsigned int i=0;

    // find the position in the list where the new item should go
    while(i<m_components.size() && m_components[i].first<pos)
	i++;

    // if at end of list, append it
    if(i==m_components.size())
    {
	m_components.push_back(make_pair<unsigned int,typename Encoding::Genotype>(pos,val));
	return;
    }

    // if there is already an item there, overwrite it
    if(m_components[i].first==pos)
    {
	m_components[i].second=val;
	return;
    }

    // otherwise, insert the new item at position i
    m_components.insert(m_components.begin()+i,
			make_pair<unsigned int,typename Encoding::Genotype>(pos, val));
}

/*!
 * \brief remove a component from the move
 */
template <template <typename> class Chromosome, typename Encoding>
void move<Chromosome,Encoding>::remove_component(unsigned int pos, typename Encoding::Genotype val)
{
    typename deque<pair<unsigned int, typename Encoding::Genotype> >::iterator it;
    
    while((it=find(m_components.begin(),m_components.end(),
		   make_pair<unsigned int, typename Encoding::Genotype>(pos,val)))!=m_components.end())
    {
	m_components.erase(it);
    }
}

/*!
 * \brief clear out all move components
 */
template <template <typename> class Chromosome, typename Encoding>
void move<Chromosome,Encoding>::reset()
{
    m_components.clear();
}

/*!
 * \brief apply the move to a chromosome
 */
template <template <typename> class Chromosome, typename Encoding>
void move<Chromosome,Encoding>::apply(Chromosome<Encoding>& chr) const
{
    for(unsigned int i=0; i<m_components.size(); i++)
    {
	pair<unsigned int, typename Encoding::Genotype> com=m_components[i];
	chr[com.first] = com.second;
    }
}

/*!
 * \brief print a human-readable representation of the move
 *
 * useful for debugging purposes
 */
template <template <typename> class Chromosome, typename Encoding>
ostream& operator<<(ostream& ostr, const move<Chromosome,Encoding>& m)
{
    ostr << "move: ";
    for(unsigned int i=0; i<m._comp.size(); i++)
    {
	ostr << "<" << m._comp[i].first << "," << m._comp[i].second << "> ";
    }
    return ostr;
}
