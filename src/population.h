/*!
 * \file population.h
 *
 * evolutionary algorithms utilize a population of candidate chromosomes.
 * this class maintains this population for use by the evolutionary
 * algorithms.
 *
 * Deon Garrett
 * deong@acm.org
 */

#ifndef _POPULATION_H_
#define _POPULATION_H_

#include <deque>
#include "chromosome.h"
#include "encoding.h"
#include "comparator.h"

using namespace std;

// forward declarations
template <template <typename> class Chromosome, typename Encoding> class population;
template <template <typename> class Chromosome, typename Encoding>
ostream& operator<<(ostream& s, const population<Chromosome,Encoding>& p);
template <template <typename> class Chromosome, typename Encoding>
istream& operator>>(istream& s, population<Chromosome,Encoding>& p);

/*!
 * \class population
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class population
{
public:
    typedef typename deque<Chromosome<Encoding> >::iterator iterator;
    typedef typename deque<Chromosome<Encoding> >::const_iterator const_iterator;
    
protected:
    // store the individuals in a deque for efficiency
    deque<Chromosome<Encoding> > m_individuals;

    // keep track of the number of individuals currently allocated
    unsigned int m_count;

public:
    // ctors and dtor
    population();
    population(unsigned int n);
    population(const population& that);
    virtual ~population();

    // assignment operator
    population& operator=(const population& that);

    // equality operator
    bool operator==(const population& that) const;
    bool operator!=(const population& that) const;

    // population size
    inline unsigned int size() const;
    
    // define + and - to do multiset union and difference
    population operator+(const population& that) const;
    population operator-(const population& that) const;

    // similarly, allow union-equal and difference-equal operators
    population operator+=(const population& that);
    population operator-=(const population& that);
    
    // access a particular individual
    inline Chromosome<Encoding>& operator[](unsigned int i);
    inline const Chromosome<Encoding>& operator[](unsigned int i) const;

    // resize the population (efficiency hack only)
    void resize(unsigned int sz);
    
    // add and remove chromosomes from the population
    virtual bool add(const Chromosome<Encoding>& ind);
    virtual void remove(const Chromosome<Encoding>& ind);
    virtual void remove_at(unsigned int index);

    // iterator functions
    iterator begin();
    iterator end();
    const_iterator begin() const;
    const_iterator end() const;
    
    // empty the population
    void clear();
    
    // sort the population by the specified comparison operator
    void sort(comparator<Chromosome,Encoding>* comp);
    
    // initialize a population by reading it from a file
    void from_file(const string& filename);
    
    friend ostream& operator<<<>(ostream& ostr, const population<Chromosome,Encoding>& p);
    friend istream& operator>><>(istream& istr, population<Chromosome,Encoding>& p);
};

/*!
 * \class bridge comparator used for sorting of population
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class sorting_comparator
{
private:
    comparator<Chromosome,Encoding>* m_comp;

public:
    sorting_comparator(comparator<Chromosome,Encoding>* comp);
    ~sorting_comparator();
    
    inline bool operator()(const Chromosome<Encoding>& c1, const Chromosome<Encoding>& c2) const;
};

#include "population.cpp"

#endif
