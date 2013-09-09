/*!
 * \file chromosome.h
 *
 * for stochastic local search, we need the notion of a candidate
 * solution to the problem.  the solution may be encoded in some form
 * (such as bit strings in genetic algorithms). 
 *
 * Deon Garrett
 * deong@acm.org
 */

#ifndef _CHROMOSOME_H_
#define _CHROMOSOME_H_

#include <iostream>
#include <vector>
#include <functional>
#include "encoding.h"
#include "problems.h"

using namespace std;

// forward declarations required for stream insertion operator
template <typename Encoding> class chromosome;
template <typename Encoding> ostream& operator<<(ostream&, const chromosome<Encoding>&);
template <typename Encoding> istream& operator>>(istream&, chromosome<Encoding>&);

/*!
 * \class basic_chromosome
 *
 * provides aspects common to all chromosome types
 */
template <typename Encoding>
class basic_chromosome
{
public:
    typedef typename Encoding::iterator iterator;
    typedef typename Encoding::const_iterator const_iterator;
    
public:
    vector<typename Encoding::FitnessType> fitness;
    
public:
    // ctors and dtor
    basic_chromosome();
    basic_chromosome(const typename Encoding::ProblemType* p);
    basic_chromosome(const basic_chromosome& that);
    virtual ~basic_chromosome();

    // assignment and comparison operators
    basic_chromosome& operator=(const basic_chromosome& that);
    bool operator!=(const basic_chromosome& that) const;
    bool operator==(const basic_chromosome& that) const; 
    bool operator<(const basic_chromosome& that) const;
    
    typename basic_chromosome<Encoding>::iterator begin();
    typename basic_chromosome<Encoding>::const_iterator begin() const;
    typename basic_chromosome<Encoding>::iterator end();
    typename basic_chromosome<Encoding>::const_iterator end() const;
    
    // direct access to the encoding
    typename Encoding::Genotype& operator[](int index);
    const typename Encoding::Genotype& operator[](int index) const;
    
    // access to important member variables
    unsigned int length() const;

    // get the genotype
    vector<typename Encoding::Genotype>& genotype();
    const vector<typename Encoding::Genotype>& genotype() const;
    
    // compute distances to other chromosomes
    double fitness_distance(const basic_chromosome& that) const;

    // determine if the chromosome is feasible (true by default)
    virtual bool feasible() const;
    
    // randomize the encoding
    void randomize();
    
    // evaluating the chromosome
    virtual void evaluate(const typename Encoding::ProblemType* prob);

protected:
    Encoding* m_parameters;
    bool m_feasible;
};

/*!
 * \class chromosome
 */
template <typename Encoding>
class chromosome : public basic_chromosome<Encoding>
{
public:
    // ctors and dtor
    chromosome();
    chromosome(const typename Encoding::ProblemType* p);
    chromosome(const chromosome& that);
    virtual ~chromosome();

    // assignment and comparison operators
    chromosome& operator=(const chromosome& that);
    bool operator!=(const chromosome& that) const;
    bool operator==(const chromosome& that) const;
    
    friend ostream& operator<<<>(ostream& ostr, const chromosome<Encoding>& sol);
    friend istream& operator>><>(istream& istr, chromosome<Encoding>& sol);
};

/*!
 * \class chromosome<integer_encoding>
 */
template <>
class chromosome<integer_encoding> : public basic_chromosome<integer_encoding>
{
public:
    // ctors and dtor
    chromosome();
    chromosome(const integer_encoding::ProblemType* p);
    chromosome(const chromosome& that);
    virtual ~chromosome();

    // assignment and comparison operators
    chromosome& operator=(const chromosome& that);
    bool operator!=(const chromosome& that) const;
    bool operator==(const chromosome& that) const;
    
    const vector<int>& legal_values(unsigned int index) const;
    virtual void evaluate(const integer_encoding::ProblemType* prob);

    friend ostream& operator<<<>(ostream& ostr, const chromosome<integer_encoding>& sol);
    friend istream& operator>><>(istream& istr, chromosome<integer_encoding>& sol);
};

/*!
 * \class chromosome<gap_encoding>
 */
template <>
class chromosome<gap_encoding> : public basic_chromosome<gap_encoding>
{
protected:
    vector<int> m_cap_used;
    vector<vector<int> > m_agt_map;
    
public:
    // ctors and dtor
    chromosome();
    chromosome(const gap_encoding::ProblemType* p);
    chromosome(const chromosome& that);
    virtual ~chromosome();

    // assignment and comparison operators
    chromosome& operator=(const chromosome& that);
    bool operator!=(const chromosome& that) const;
    bool operator==(const chromosome& that) const;

    void compute_infeasibility(const gap_encoding::ProblemType* p);

    unsigned int agents() const;
    unsigned int tasks() const;
    virtual void penalize(const gap_encoding::ProblemType* prob);
    virtual void repair(const gap_encoding::ProblemType* prob);
    virtual void evaluate(const gap_encoding::ProblemType* prob);

    friend ostream& operator<<<>(ostream& ostr, const chromosome<gap_encoding>& sol);
    friend istream& operator>><>(istream& istr, chromosome<gap_encoding>& sol);
};

/*!
 * \class chromosome<gsap_encoding>
 */
template <>
class chromosome<gsap_encoding> : public basic_chromosome<gsap_encoding>
{
protected:
    vector<unsigned int> m_cap_used;
    unsigned int n_unassigned;
	unsigned int m_excess;
	
public:
    // ctors and dtor
    chromosome();
    chromosome(const gsap_encoding::ProblemType* p);
    chromosome(const chromosome& that);
    virtual ~chromosome();

    // assignment and comparison operators
    chromosome& operator=(const chromosome& that);
    bool operator!=(const chromosome& that) const;
    bool operator==(const chromosome& that) const;

    void compute_infeasibility(const gsap_encoding::ProblemType* p);

    unsigned int agents() const;
    unsigned int tasks() const;
    virtual void penalize(const gsap_encoding::ProblemType* prob);
    virtual void repair(const gsap_encoding::ProblemType* prob);
    virtual void evaluate(const gsap_encoding::ProblemType* prob);

    friend ostream& operator<<<>(ostream& ostr, const chromosome<gsap_encoding>& sol);
    friend istream& operator>><>(istream& istr, chromosome<gsap_encoding>& sol);
};

#include "chromosome.cpp"

#endif

