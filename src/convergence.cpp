/*!
 * \file convergence.cpp
 *
 * classes which are used to determine whether an evolutionary
 * algorithm has converged
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <string>
#include <list>
#include "convergence.h"
#include "chromosome.h"
#include "population.h"

using namespace std;

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
convergence<Chromosome,Encoding>::convergence()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
convergence<Chromosome,Encoding>::~convergence()
{
}

/*!
 * \brief callback for when chromosome is evaluated
 */
template <template <typename> class Chromosome, typename Encoding>
void convergence<Chromosome,Encoding>::chromosome_evaluated(const Chromosome<Encoding>& sol)
{
}

/*!
 * \brief callback for when generation is completed
 */
template <template <typename> class Chromosome, typename Encoding>
void convergence<Chromosome,Encoding>::generation_completed(const population<Chromosome,Encoding>& pop)
{
}

/*!
 * \brief create and return an initialized convergence operator
 */
template <template <typename> class Chromosome, typename Encoding>
list<convergence<Chromosome,Encoding>*> convergence_factory<Chromosome,Encoding>::construct()
{
	list<convergence<Chromosome,Encoding>*> conv;
	return conv;
}
