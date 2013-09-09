/*!
 * \file terminators.h
 *
 * in stochastic local search, there is typically no one method by
 * which we can decide when to terminate the search.  these classes
 * provide several alternative mechanism by which the search process
 * can be halted.
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _TERMINATORS_H_
#define _TERMINATORS_H_

#include <list>
#include <string>
#include "chromosome.h"
#include "population.h"
#include "factory.h"

using namespace std;

/*!
 * \class terminator
 *
 * basic class representing a criteria for halting a search algorithm
 */
template <template <typename> class Chromosome, typename Encoding>
class terminator
{
private:
    // disable copying of terminators
    terminator(const terminator& that);
    terminator& operator=(const terminator& that);

private:
    string m_prefix;
    
public:
    terminator();
    virtual ~terminator();

    virtual void initialize(const string& prefix);
    virtual void reset()=0;
    virtual void chromosome_evaluated(const Chromosome<Encoding>& sol);
    virtual void generation_completed();
    virtual void generation_completed(const population<Chromosome,Encoding>& pop);
    virtual bool terminate() const = 0;
};

/*!
 * \class evaluation_limit
 */
template <template <typename> class Chromosome, typename Encoding>
class evaluation_limit : public terminator<Chromosome,Encoding>
{
private:
    unsigned int m_evals;
    unsigned int m_max_evals;

    // disable copying of terminators
    evaluation_limit(const evaluation_limit& that);
    evaluation_limit& operator=(const evaluation_limit& that);

public:
    evaluation_limit();
    virtual ~evaluation_limit();
    virtual void initialize(const string& prefix);
    virtual void reset();
    virtual void chromosome_evaluated(const Chromosome<Encoding>& sol);
    virtual bool terminate() const;
};

/*!
 * \class generation_limit
 */
template <template <typename> class Chromosome, typename Encoding>
class generation_limit : public terminator<Chromosome,Encoding>
{
private:
    unsigned int m_gens;
    unsigned int m_max_gens;

    // disable copying of terminators
    generation_limit(const generation_limit& that);
    generation_limit& operator=(const generation_limit& that);

public:
    generation_limit();
    virtual ~generation_limit();
    virtual void initialize(const string& prefix);
    virtual void reset();
    virtual void generation_completed();
    virtual void generation_completed(const population<Chromosome,Encoding>& pop);
    virtual bool terminate() const;
};

/*!
 * \class null_terminator
 */
template <template <typename> class Chromosome, typename Encoding>
class null_terminator : public terminator<Chromosome,Encoding>
{
public:
    null_terminator();
    virtual ~null_terminator();
    virtual void reset();
    virtual bool terminate() const;
};


/*!
 * \class terminator_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class terminator_factory : public factory
{
public:
    list<terminator<Chromosome,Encoding>*> construct();
};

#include "terminators.cpp"

#endif
