/*! 
 * \file comparator.h
 *
 * defines classes which determine how chromosomes are compared
 * to one another by fitness.
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _COMPARATOR_H_
#define _COMPARATOR_H_

#include "chromosome.h"
#include "factory.h"

/*!
 * \class comparator
 */
template <template <typename> class Chromosome, typename Encoding>
class comparator
{
public:
    comparator();
    virtual ~comparator();
    virtual void initialize();
    inline bool operator()(const Chromosome<Encoding>& c1, const Chromosome<Encoding>& c2) const;
    virtual int compare(const Chromosome<Encoding>& c1, const Chromosome<Encoding>& c2) const = 0;
};

/*!
 * \class fitness_comparator
 *
 * compare chromosomes by value of first objective function.  used to fake
 * single objective optimization.
 */
template <template <typename> class Chromosome, typename Encoding>
class fitness_comparator : public comparator<Chromosome,Encoding>
{
public:
    fitness_comparator();
    virtual ~fitness_comparator();

    inline virtual int compare(const Chromosome<Encoding>& c1, const Chromosome<Encoding>& c2) const;
};

/*!
 * \class pareto_dominance_comparator
 */
template <template <typename> class Chromosome, typename Encoding>
class pareto_dominance_comparator : public comparator<Chromosome,Encoding>
{
public:
    pareto_dominance_comparator();
    virtual ~pareto_dominance_comparator();

    inline virtual int compare(const Chromosome<Encoding>& c1, const Chromosome<Encoding>& c2) const;
};

/*!
 * \class scalarizing_comparator
 */
template <template <typename> class Chromosome, typename Encoding>
class scalarizing_comparator : public comparator<Chromosome,Encoding>
{
public:
    vector<double> weights;

public:
    scalarizing_comparator();
    virtual ~scalarizing_comparator();

    void randomize_weights(unsigned int nobj);
    virtual void initialize(const string& prefix);
    virtual void initialize();
    inline virtual int compare(const Chromosome<Encoding>& c1, const Chromosome<Encoding>& c2) const;
    inline double difference(const Chromosome<Encoding>& c1, const Chromosome<Encoding>& c2) const;
};

/*!
 * \class weak_dominance_comparator
 */
template <template <typename> class Chromosome, typename Encoding>
class weak_dominance_comparator : public comparator<Chromosome,Encoding>
{
public:
    weak_dominance_comparator();
    virtual ~weak_dominance_comparator();

    inline virtual int compare(const Chromosome<Encoding>& c1, const Chromosome<Encoding>& c2) const;
};

/*!
 * \class strong_dominance_comparator
 */
template <template <typename> class Chromosome, typename Encoding>
class strong_dominance_comparator : public comparator<Chromosome,Encoding>
{
public:
    strong_dominance_comparator();
    virtual ~strong_dominance_comparator();

    inline virtual int compare(const Chromosome<Encoding>& c1, const Chromosome<Encoding>& c2) const;
};

/*!
 * \class epsilon_dominance_comparator
 */
template <template <typename> class Chromosome, typename Encoding>
class epsilon_dominance_comparator : public comparator<Chromosome,Encoding>
{
protected:
    vector<double> m_epsilon;
    
public:
    epsilon_dominance_comparator();
    virtual ~epsilon_dominance_comparator();
    virtual void initialize(const string& prefix);
    virtual void initialize();
    inline virtual int compare(const Chromosome<Encoding>& c1, const Chromosome<Encoding>& c2) const;
};

/*!
 * \class single_objective_comparator
 *
 * more general case of fitness_comparator where client can specify on which
 * objective function to base the comparison
 */
template <template <typename> class Chromosome, typename Encoding>
class single_objective_comparator : public comparator<Chromosome,Encoding>
{
protected:
    unsigned int m_obj;

public:
    single_objective_comparator(unsigned int objnum);
    virtual ~single_objective_comparator();

    inline virtual int compare(const Chromosome<Encoding>& c1, const Chromosome<Encoding>& c2) const;
};

/*!
 * \class comparator_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class comparator_factory : public factory
{
public:
    comparator<Chromosome,Encoding>* construct();
};

#include "comparator.cpp"

#endif
