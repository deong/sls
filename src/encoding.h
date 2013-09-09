/*!
 * \file encoding.h
 *
 * for stochastic local search, we need the notion of a candidate
 * solution to the problem.  the solution may be encoded in some form
 * (such as bit strings in genetic algorithms).  these classes provide
 * the representations for several different types of encodings
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _ENCODING_H_
#define _ENCODING_H_

#include <iostream>
#include <vector>
#include "problems.h"

using namespace std;

/*!
 * \class encoding
 */
template <typename G, typename P>
class encoding
{
public:
    typedef G Genotype;
    typedef P Phenotype;

protected:
    vector<G> m_genotype;
    vector<P> m_phenotype;

public:
    encoding();
    encoding(const problem* p);
    virtual ~encoding();

    // assignment and comparison operators
    encoding& operator=(const encoding& that);
    bool operator==(const encoding& that);

    // direct access to the encoding information
    inline G& operator[](unsigned int index);
    inline const G& operator[](unsigned int index) const;
    
    // direct access to the information stored here
    vector<G>& genotype();
    const vector<G>& genotype() const;
    vector<P>& phenotype();
    const vector<P>& phenotype() const;

    // retrieve the encoding length
    inline unsigned int length() const;

    virtual void randomize() = 0;
    virtual void decode() = 0;
};

/*!
 * \class bit_vector_encoding
 *
 * bit_vector_encoding is used for any encoding that uses
 * a vector of bits for the genotype
 */
template <typename P>
class bit_vector_encoding : public encoding<int,P>
{
public:
    typedef vector<int>::iterator iterator;
    typedef vector<int>::const_iterator const_iterator;
    
public:
    bit_vector_encoding();
    bit_vector_encoding(const problem* p);
    virtual ~bit_vector_encoding();
    typename bit_vector_encoding<P>::iterator begin();
    typename bit_vector_encoding<P>::const_iterator begin() const;
    typename bit_vector_encoding<P>::iterator end();
    typename bit_vector_encoding<P>::const_iterator end() const;
    virtual void randomize();
    virtual void decode() = 0;
};

/*!
 * \class boolean_encoding
 *
 * boolean_encoding is suitable for problems in which each parameter
 * can be represented as a boolean value (knapsack problems for
 * example)
 */
class boolean_encoding : public bit_vector_encoding<int>
{
public:
    typedef bit_string_problem ProblemType;
    typedef bit_string_problem_factory ProblemFactoryType;
    typedef int FitnessType;
    
public:
    boolean_encoding();
    boolean_encoding(const bit_string_problem* p);
    virtual ~boolean_encoding();
    static void initialize_parameters(const bit_string_problem* prob);
    static void clear_parameters();
    inline virtual void encode(unsigned int p);
    inline virtual void decode();
    friend ostream& operator<<(ostream& s, const boolean_encoding& e);
    friend istream& operator>>(istream& s, boolean_encoding& e);
};

/*!
 * \class numeric_parameters
 *
 * helper class for numeric encodings
 */
class numeric_parameters
{
public:
    static void initialize_parameters(const numeric_problem* p);
    static void cleanup();
};

/*!
 * \class numeric_encoding
 *
 * numeric_encoding is a mixin which adds parameter range bounds for
 * real-valued encodings
 */
class numeric_encoding
{
protected:
    static vector<pair<double,double> > m_range;
    
public:
    numeric_encoding();
    numeric_encoding(const numeric_problem* p);
    virtual ~numeric_encoding();
    static const pair<double,double>& parameter_range(int pnum);
    friend class numeric_parameters;
};

/*!
 * \class binary_parameters
 *
 * helper class for binary encoding
 */
class binary_parameters
{
public:
    static void initialize_parameters(const numeric_problem* p);
    static void cleanup();
};

/*!
 * \class binary_encoding
 *
 * binary_encoding encodes a real-valued parameter vector into a
 * binary string
 */
class binary_encoding : public bit_vector_encoding<double>, public numeric_encoding
{
public:
    typedef numeric_problem ProblemType;
    typedef numeric_problem_factory ProblemFactoryType;
    typedef double FitnessType;
    
protected:
    // how many bits store each parameter
    static vector<int> m_bpp;

    // for simplicity, keep track of the total length of the bit string
    static unsigned int m_len;
    
    // is the string gray coded
    static vector<bool> m_gray;

public:
    binary_encoding();
    binary_encoding(const numeric_problem* p);
    virtual ~binary_encoding();
    static void initialize_parameters(const numeric_problem* p);
    static void clear_parameters();
    inline virtual void encode(const vector<double>& params);
    inline virtual void decode();
    friend ostream& operator<<(ostream& s, const binary_encoding& e);
    friend istream& operator>>(istream& s, binary_encoding& e);
    friend class binary_parameters;
};

/*!
 * \class real_encoding
 *
 * directly stores real numbers as part of the encoding
 */
class real_encoding : public encoding<double,double>, public numeric_encoding
{
public:
    typedef numeric_problem ProblemType;
    typedef numeric_problem_factory ProblemFactoryType;
    typedef double FitnessType;
    typedef vector<double>::iterator iterator;
    typedef vector<double>::const_iterator const_iterator;
    
public:
    real_encoding();
    real_encoding(const numeric_problem* p);
    virtual ~real_encoding();
    static void initialize_parameters(const numeric_problem* p);
    static void clear_parameters();
    real_encoding::iterator begin();
    real_encoding::const_iterator begin() const;
    real_encoding::iterator end();
    real_encoding::const_iterator end() const;
    virtual void randomize();
    inline virtual void decode();
    friend ostream& operator<<(ostream& s, const real_encoding& e);
    friend istream& operator>>(istream& s, real_encoding& e);
};

/*!
 * \class permutation_encoding
 *
 * permutation_encoding is suitable for combinatorial optimization
 * problems whose chromosomes can be represented as permutations (tsp
 * for example)
 */
class permutation_encoding : public encoding<int,int>
{
public:
    typedef permutation_problem ProblemType;
    typedef permutation_problem_factory ProblemFactoryType;
    typedef long FitnessType;
    typedef vector<int>::iterator iterator;
    typedef vector<int>::const_iterator const_iterator;
    
public:
    permutation_encoding();
    permutation_encoding(const permutation_problem* p);
    virtual ~permutation_encoding();
    permutation_encoding::iterator begin();
    permutation_encoding::const_iterator begin() const;
    permutation_encoding::iterator end();
    permutation_encoding::const_iterator end() const;
    static void initialize_parameters(const permutation_problem* p);
    static void clear_parameters();
    virtual void randomize();
    inline virtual void decode();
    friend ostream& operator<<(ostream& s, const permutation_encoding& e);
    friend istream& operator>>(istream& s, permutation_encoding& e);
};

/*!
 * \class integer_encoding
 *
 * integer_encoding is suitable for problems characterized by vectors
 * of integers.  it supports constrained integer encodings through a
 * method denoting the acceptable values at each locus
 */
class integer_encoding : public encoding<int,int>
{
public:
    typedef integer_problem ProblemType;
    typedef integer_problem_factory ProblemFactoryType;
    typedef int FitnessType;
    typedef vector<int>::iterator iterator;
    typedef vector<int>::const_iterator const_iterator;
    
protected:
    static vector<vector<int> > m_legal_values;
    
public:
    integer_encoding();
    integer_encoding(const integer_problem* p);
    virtual ~integer_encoding();
    integer_encoding::iterator begin();
    integer_encoding::const_iterator begin() const;
    integer_encoding::iterator end();
    integer_encoding::const_iterator end() const;
    static void initialize_parameters(const integer_problem* p);
    static void clear_parameters();
    const vector<int>& legal_values(unsigned int index) const;
    virtual void randomize();
    inline virtual void decode();
    friend ostream& operator<<(ostream& s, const integer_encoding& e);
    friend istream& operator>>(istream& s, integer_encoding& e);
};

/*!
 * \class gap_encoding
 *
 * specialized encoding for generalized assignment problem
 */
class gap_encoding : public encoding<int,int>
{
public:
    typedef gap_problem ProblemType;
    typedef gap_problem_factory ProblemFactoryType;
    typedef int FitnessType;
    typedef vector<int>::iterator iterator;
    typedef vector<int>::const_iterator const_iterator;
    
public:
    static vector<double> m_alpha;
    
protected:
    static unsigned int m_agents;
    static unsigned int m_tasks;
    
public:
    gap_encoding();
    gap_encoding(const gap_problem* p);
    virtual ~gap_encoding();
    unsigned int agents() const;
    unsigned int tasks() const;
    gap_encoding::iterator begin();
    gap_encoding::const_iterator begin() const;
    gap_encoding::iterator end();
    gap_encoding::const_iterator end() const;
    virtual void randomize();
    inline virtual void decode();
    static void initialize_parameters(const gap_problem* p);
    static void clear_parameters();
    friend ostream& operator<<(ostream& s, const gap_encoding& e);
    friend istream& operator>>(istream& s, gap_encoding& e);
};

/**
 * \class gsap_encoding
 * \brief encoding for the generalized sailor assignment problem
 */
class gsap_encoding : public encoding<int,int>
{
public:
    typedef gsap_problem ProblemType;
    typedef gsap_problem_factory ProblemFactoryType;
    typedef int FitnessType;
    typedef vector<int>::iterator iterator;
    typedef vector<int>::const_iterator const_iterator;

public:
    static int m_unass_pen;
    static int m_cap_pen;
    
protected:
    static vector<gsap_problem::task_element> elements;
    static vector<vector<int> > agents_for_task;
    
public:
    gsap_encoding();
    gsap_encoding(const ProblemType* p);
    virtual ~gsap_encoding();
    gsap_encoding::iterator begin();
    gsap_encoding::const_iterator begin() const;
    gsap_encoding::iterator end();
    gsap_encoding::const_iterator end() const;
    virtual void randomize();
    inline virtual void decode();
    static void initialize_parameters(const gsap_problem* p);
    static void clear_parameters();
    friend ostream& operator<<(ostream& s, const gsap_encoding& e);
    friend istream& operator>>(istream& s, gsap_encoding& e);
};

#include "encoding.cpp"

#endif
