/*!
 * \file move.h
 *
 * attempt to abstract away the notion of a move in a local search
 * algorithm
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _MOVE_H_
#define _MOVE_H_

#include <deque>
#include <utility>
#include "encoding.h"
#include "problems.h"

using namespace std;

template <template <typename> class Chromosome, typename Encoding> class move;
template <template <typename> class Chromosome, typename Encoding>
ostream& operator<<(ostream& ostr, const move<Chromosome,Encoding>& m);

/*!
 * \class move
 */
template <template <typename> class Chromosome, typename Encoding>
class move
{
public:
	typedef typename deque<pair<unsigned int,typename Encoding::Genotype> >::iterator iterator;
	typedef typename deque<pair<unsigned int,typename Encoding::Genotype> >::const_iterator const_iterator;

protected:
	deque<pair<unsigned int,typename Encoding::Genotype> > m_components;

public:
	move();
	~move();
	bool operator==(const move<Chromosome,Encoding>& that) const;
	bool operator!=(const move<Chromosome,Encoding>& that) const;
	unsigned int length() const;
	pair<unsigned int,typename Encoding::Genotype>& operator[](unsigned int i);
	const pair<unsigned int,typename Encoding::Genotype>& operator[](unsigned int i) const;
	typename move<Chromosome,Encoding>::iterator begin();
	typename move<Chromosome,Encoding>::iterator end();
	typename move<Chromosome,Encoding>::const_iterator begin() const;
	typename move<Chromosome,Encoding>::const_iterator end() const;
	void add_component(unsigned int pos, typename Encoding::Genotype val);
	void remove_component(unsigned int pos, typename Encoding::Genotype val);
	void apply(Chromosome<Encoding>& chr) const;
	void reset();

public:
	friend ostream& operator<<<>(ostream& ostr, const move<Chromosome,Encoding>& m);
};

#include "move.cpp"

#endif
