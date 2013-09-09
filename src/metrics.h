/*!
 * \file metrics.h
 *
 * the performance of a search algorithm can be characterised by
 * any of several performance metrics.  these classes define some
 * useful metrics for stochastic local search.
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _METRICS_H_
#define _METRICS_H_

#include <iostream>
#include <deque>
#include <list>
#include <vector>
#include <string>
#include "chromosome.h"
#include "population.h"
#include "comparator.h"
#include "factory.h"

using namespace std;

/*!
 * \class metric
 */
template <template <typename> class Chromosome, typename Encoding>
class metric
{
private:
    // prevent copying of performance metrics
    metric(const metric& that);
    metric& operator=(const metric& that);

public:
    metric();
    virtual ~metric();

    virtual void chromosome_evaluated(const Chromosome<Encoding>& sol);
    virtual void generation_completed();
    virtual void generation_completed(const population<Chromosome,Encoding>& sol);
    virtual void reset()=0;
    virtual void compute();
    virtual void report(ostream& ostr) const = 0;
};

/*!
 * \class evaluation_counter
 */
template <template <typename> class Chromosome, typename Encoding>
class evaluation_counter : public metric<Chromosome,Encoding>
{
private:
    unsigned int m_num_evals;
    
    // prevent copying of performance metrics
    evaluation_counter(const evaluation_counter& that);
    evaluation_counter& operator=(const evaluation_counter& that);

public:
    evaluation_counter();
    virtual ~evaluation_counter();
    virtual void reset();
    virtual void chromosome_evaluated(const Chromosome<Encoding>& sol);
    virtual void report(ostream& ostr) const;
};

/*!
 * \class best_solution
 */
template <template <typename> class Chromosome, typename Encoding>
class best_solution : public metric<Chromosome,Encoding>
{
private:
    Chromosome<Encoding> m_best;
    unsigned int  m_num_evals;
    unsigned int  m_evals_to_best;
	bool m_report_all_best;
    comparator<Chromosome,Encoding>* m_comp;
    
    // prevent copying of performance metrics
    best_solution(const best_solution& that);
    best_solution& operator=(const best_solution& that);

public:
    best_solution();
    virtual ~best_solution();
    virtual void initialize();
    virtual void initialize(const string& prefix);
    virtual void reset();
    virtual void chromosome_evaluated(const Chromosome<Encoding>& sol);
    virtual void report(ostream& ostr) const;
};

/*!
 * \class generation_counter
 */
template <template <typename> class Chromosome, typename Encoding>
class generation_counter : public metric<Chromosome,Encoding>
{
private:
    unsigned int m_num_gens;

    // prevent copying of performance metrics
    generation_counter(const generation_counter& that);
    generation_counter& operator=(const generation_counter& that);

public:
    generation_counter();
    virtual ~generation_counter();
    virtual void reset();
    virtual void generation_completed(const population<Chromosome,Encoding>& pop);
    virtual void generation_completed();
    virtual void report(ostream& ostr) const;
};

/*!
 * \class hypervolume
 */
template <template <typename> class Chromosome, typename Encoding>
class hypervolume : public metric<Chromosome,Encoding>
{
protected:
    vector<typename Encoding::FitnessType> m_ref_point;
    double m_hypervolume;
    population<Chromosome,Encoding> m_pop;
    
private:
    // disable copying of performance metrics
    hypervolume(const hypervolume& that);
    hypervolume& operator=(const hypervolume& that);

    bool weakly_dominated_by_archive(const vector<typename Encoding::FitnessType>& p,
                                     deque<vector<typename Encoding::FitnessType> >& ps) const;
    bool weakly_dominated(const vector<typename Encoding::FitnessType>& p1,
                          const vector<typename Encoding::FitnessType>& p2) const;
    
public:
    hypervolume();
    virtual ~hypervolume();
    virtual void initialize(const string& prefix);
    virtual void initialize();
    virtual void reset();
    virtual void generation_completed(const population<Chromosome,Encoding>& pop);
    virtual void compute();
    virtual void report(ostream& ostr) const;
};

/**
 * \class population_entropy
 * \brief measures the average entropy in the population over all genes
 */
template <template <typename> class Chromosome, typename Encoding>
class population_entropy : public metric<Chromosome,Encoding>
{
public:
    population_entropy();
    virtual ~population_entropy();
	
	virtual void reset();
	virtual void generation_completed(const population<Chromosome,Encoding>& pop);
	virtual void report(ostream& ostr) const;
	
private:
	// disable copying
	population_entropy(const population_entropy& that);
	population_entropy& operator=(const population_entropy& that);
	
private:
	double m_avg_entropy;
};

/*!
 * \class metric_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class metric_factory : public factory
{
public:
    list<metric<Chromosome,Encoding>*> construct();
};

#include "metrics.cpp"

#endif
