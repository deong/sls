/*!
 * \file landscape.h
 *
 * simple landscape analysis of the mQAP
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _LANDSCAPE_H_
#define _LANDSCAPE_H_

#include "sls.h"
#include "localsearch.h"
#include "pfront.h"
#include "comparator.h"
#include "metrics.h"
#include "terminators.h"

/*!
 * \class landscape
 *
 * base class for landscape analysis tools
 */
template <template <typename> class Chromosome, typename Encoding>
class landscape : public sls<Chromosome,Encoding>
{
public:
    enum landscape_tool
    {
        PERTURB_WEIGHT_VECTOR,
        FITNESS_DISTANCE_CORRELATION,
        RANDOM_WALK,
        RANDOM_WALK_BETWEEN_OPTIMA,
        PARETO_PLATEAUS
    };

protected:
    local_search<Chromosome,Encoding>* m_hc;
    unsigned int m_lsiter;

public:
    landscape();
    virtual ~landscape();

    virtual void initialize();
    virtual void run() = 0;
};

/*!
 * \class fitness_distance_correlation
 */
template <template <typename> class Chromosome, typename Encoding>
class fitness_distance_correlation : public landscape<Chromosome,Encoding>
{
public:
    virtual void run();
};

/*!
 * \class random_walk
 */
template <template <typename> class Chromosome, typename Encoding>
class random_walk : public landscape<Chromosome,Encoding>
{
public:
    virtual void run();
};

/*!
 * \class random_walk_between_optima
 */
template <template <typename> class Chromosome, typename Encoding>
class random_walk_between_optima : public landscape<Chromosome,Encoding>
{
public:
    virtual void run();
};

/*!
 * \class ruggedness
 */
template <template <typename> class Chromosome, typename Encoding>
class ruggedness : public landscape<Chromosome,Encoding>
{
public:
    virtual void run();
};

/*!
 * \class ruggedness_between_optima
 */
template <template <typename> class Chromosome, typename Encoding>
class ruggedness_between_optima : public landscape<Chromosome,Encoding>
{
public:
    virtual void run();
};

/*!
 * \class perturbation_search
 */
template <template <typename> class Chromosome, typename Encoding>
class perturbation_search : public landscape<Chromosome,Encoding>
{
public:
    virtual void run();
};

/*!
 * \class infeasibility_region
 */
template <template <typename> class Chromosome, typename Encoding>
class infeasibility_region : public landscape<Chromosome,Encoding>
{
public:
    virtual void run();
};

/*!
 * \class hypervolume_analysis
 */
template <template <typename> class Chromosome, typename Encoding>
class hypervolume_analysis : public landscape<Chromosome,Encoding>
{
public:
    virtual void initialize();
    virtual void run();
};

/*!
 * \class pareto_plateaus
 * \brief computes pareto plateau structure
 */
template <template <typename> class Chromosome, typename Encoding>
class pareto_plateaus : public landscape<Chromosome,Encoding>
{
public:
    virtual void initialize();
    virtual void run();
protected:
    void nondominated_sort(population<Chromosome,Encoding>& pop,
                           vector<pareto_front<Chromosome,Encoding> >& fronts);
};

/*!
 * \class landscape_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class landscape_factory
{
public:
    static landscape<Chromosome,Encoding>* construct();
};

#include "landscape.cpp"

#endif

    
