/*!
 * \file localsearch.h
 *
 * local search algorithms
 *
 * Deon Garrett
 * deong@acm.org
 */

#ifndef _LOCALSEARCH_H_
#define _LOCALSEARCH_H_

#include "sls.h"
#include "chromosome.h"
#include "comparator.h"
#include "neighborhood.h"
#include "problems.h"
#include "metrics.h"
#include "terminators.h"
#include "factory.h"
#include "repair.h"

/*!
 * \class local_search
 *
 * \author deong
 * \date 05/09/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class local_search
{
private:
    // disable copying
    local_search(const local_search& that);
    local_search& operator=(const local_search& that);

protected:
    neighborhood<Chromosome,Encoding>* m_nf;
    list<metric<Chromosome,Encoding>*> m_metrics;
    list<terminator<Chromosome,Encoding>*> m_terminators;
    repair_operator<Chromosome,Encoding>* m_repair;
    bool debug_generations;
    string m_prefix;
    
public:
    local_search();
    virtual ~local_search();
    void chromosome_evaluated(const Chromosome<Encoding>& chr);
    void generation_completed();
    virtual void set_prefix(const string& prefix);
    virtual bool terminate() const;
    virtual void initialize();
    virtual void reset();
    neighborhood<Chromosome,Encoding>* get_neighborhood();
    virtual void improve(Chromosome<Encoding>& chr,
                         comparator<Chromosome,Encoding>* comp,
                         const typename Encoding::ProblemType* prob)=0;
};

/*!
 * \class next_descent
 *
 * \author deong
 * \date 05/09/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class next_descent : public local_search<Chromosome,Encoding>
{
private:
    // disable copying
    next_descent(const next_descent& that);
    next_descent& operator=(const next_descent& that);

public:
    next_descent();
    virtual ~next_descent();

    virtual void improve(Chromosome<Encoding>& chr,
                         comparator<Chromosome,Encoding>* comp,
                         const typename Encoding::ProblemType* prob);
};

/*!
 * \class steepest_descent
 *
 * \author deong
 * \date 05/09/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class steepest_descent : public local_search<Chromosome,Encoding>
{
private:
    // disable copying
    steepest_descent(const steepest_descent& that);
    steepest_descent& operator=(const steepest_descent& that);

public:
    steepest_descent();
    virtual ~steepest_descent();

    virtual void improve(Chromosome<Encoding>& chr,
                         comparator<Chromosome,Encoding>* comp,
                         const typename Encoding::ProblemType* prob);
};

/*!
 * \class local_search_factory
 *
 * \author deong
 * \date 05/09/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class local_search_factory : public factory
{
public:
    local_search<Chromosome,Encoding>* construct();
};

#include "localsearch.cpp"

#endif
