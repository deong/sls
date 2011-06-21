/*!
 * \file landscape.h
 *
 * simple landscape analysis of the mQAP
 *
 * Deon Garrett
 * deong@acm.org
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
 *
 * \author deong
 * \date 05/09/2007
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
 *
 * \author deong
 * \date 05/09/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class fitness_distance_correlation : public landscape<Chromosome,Encoding>
{
public:
    virtual void run();
};

/*!
 * \class random_walk
 *
 * \author deong
 * \date 05/09/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class random_walk : public landscape<Chromosome,Encoding>
{
public:
    virtual void run();
};

/*!
 * \class random_walk_between_optima
 *
 * \author deong
 * \date 05/09/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class random_walk_between_optima : public landscape<Chromosome,Encoding>
{
public:
    virtual void run();
};

/*!
 * \class ruggedness
 *
 * \author deong
 * \date 05/09/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class ruggedness : public landscape<Chromosome,Encoding>
{
public:
    virtual void run();
};

/*!
 * \class ruggedness_between_optima
 *
 * \author deong
 * \date 05/09/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class ruggedness_between_optima : public landscape<Chromosome,Encoding>
{
public:
    virtual void run();
};

/*!
 * \class perturbation_search
 *
 * \author deong
 * \date 05/09/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class perturbation_search : public landscape<Chromosome,Encoding>
{
public:
    virtual void run();
};

/*!
 * \class infeasibility_region
 *
 * \author deong
 * \date 05/09/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class infeasibility_region : public landscape<Chromosome,Encoding>
{
public:
    virtual void run();
};

/*!
 * \class hypervolume_analysis
 *
 * \author deong
 * \date 06/06/2007
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
 *
 * \author deong@acm.org
 * \date 10/14/2008
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
 *
 * \author deong
 * \date 05/09/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class landscape_factory
{
public:
    static landscape<Chromosome,Encoding>* construct();
};

#include "landscape.cpp"

#endif

    
