/*!
 * \file convergence.cpp
 *
 * classes which are used to determine whether an evolutionary
 * algorithm has converged
 *
 * Deon Garrett
 * deong@acm.org
 */

#include <string>
#include <list>
#include "convergence.h"
#include "chromosome.h"
#include "population.h"

using namespace std;

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
convergence<Chromosome,Encoding>::convergence()
{
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
convergence<Chromosome,Encoding>::~convergence()
{
}

/*!
 * \brief callback for when chromosome is evaluated
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void convergence<Chromosome,Encoding>::chromosome_evaluated(const Chromosome<Encoding>& sol)
{
}

/*!
 * \brief callback for when generation is completed
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void convergence<Chromosome,Encoding>::generation_completed(const population<Chromosome,Encoding>& pop)
{
}

/*!
 * \brief create and return an initialized convergence operator
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
list<convergence<Chromosome,Encoding>*> convergence_factory<Chromosome,Encoding>::construct()
{
    list<convergence<Chromosome,Encoding>*> conv;
    return conv;
}
