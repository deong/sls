/*!
 * \file lsmove.h
 *
 * attempt to abstract away the notion of a move in a local search
 * algorithm
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _LSMOVE_H_
#define _LSMOVE_H_

#include <deque>
#include <utility>
#include <ostream>
#include "encoding.h"
#include "problems.h"

using namespace std;

template <template <typename> class Chromosome, typename Encoding> class lsmove;
template <template <typename> class Chromosome, typename Encoding>
ostream& operator<<(ostream& ostr, const lsmove<Chromosome,Encoding>& m);

/*!
 * \class lsmove
 */
template <template <typename> class Chromosome, typename Encoding>
class lsmove
{
public:
	typedef typename deque<pair<unsigned int,typename Encoding::Genotype> >::iterator iterator;
	typedef typename deque<pair<unsigned int,typename Encoding::Genotype> >::const_iterator const_iterator;

protected:
	deque<pair<unsigned int,typename Encoding::Genotype> > m_components;

public:
	lsmove();
	~lsmove();
	bool operator==(const lsmove<Chromosome,Encoding>& that) const;
	bool operator!=(const lsmove<Chromosome,Encoding>& that) const;
	unsigned int length() const;
	pair<unsigned int,typename Encoding::Genotype>& operator[](unsigned int i);
	const pair<unsigned int,typename Encoding::Genotype>& operator[](unsigned int i) const;
	typename lsmove<Chromosome,Encoding>::iterator begin();
	typename lsmove<Chromosome,Encoding>::iterator end();
	typename lsmove<Chromosome,Encoding>::const_iterator begin() const;
	typename lsmove<Chromosome,Encoding>::const_iterator end() const;
	void add_component(unsigned int pos, typename Encoding::Genotype val);
	void remove_component(unsigned int pos, typename Encoding::Genotype val);
	void apply(Chromosome<Encoding>& chr) const;
	void reset();

public:
	friend ostream& operator<<<>(ostream& ostr, const lsmove<Chromosome,Encoding>& m);
};

#include "lsmove.cpp"

#endif
