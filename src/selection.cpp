/*!
 * \file selection.cpp
 *
 * evolutionary algorithms select individuals for reproduction based
 * upon fitness.  these classes provide commonly used mechanisms for
 * performing this selection.
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <string>
#include "selection.h"
#include "chromosome.h"
#include "encoding.h"
#include "population.h"
#include "comparator.h"
#include "mtrandom.h"
#include "kvparse/kvparse.h"
#include "keywords.h"
#include "utilities.h"

using namespace std;

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
selection_scheme<Chromosome,Encoding>::selection_scheme()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
selection_scheme<Chromosome,Encoding>::~selection_scheme()
{
}

/*!
 * \brief set the mating pool from the current population
 */
template <template <typename> class Chromosome, typename Encoding>
void selection_scheme<Chromosome,Encoding>::set_population(const population<Chromosome,Encoding>* pop)
{
	m_population = const_cast<population<Chromosome,Encoding>*>(pop);
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
tournament_selection<Chromosome,Encoding>::tournament_selection() :
	m_comp(0),
	m_delete_comp(false)
{
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
tournament_selection<Chromosome,Encoding>::tournament_selection(comparator<Chromosome,Encoding>* comp)
{
	m_comp = comp;
	m_delete_comp = false;
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
tournament_selection<Chromosome,Encoding>::~tournament_selection()
{
	if(m_delete_comp) {
		delete m_comp;
	}
}

/*!
 * \brief initialize the selection comparator object
 */
template <template <typename> class Chromosome, typename Encoding>
void tournament_selection<Chromosome,Encoding>::initialize()
{
	comparator_factory<Chromosome,Encoding> cf;
	m_comp = cf.construct();
	m_delete_comp = true;
}

/*!
 * \brief select a parent from the gene pool
 */
template <template <typename> class Chromosome, typename Encoding>
Chromosome<Encoding>& tournament_selection<Chromosome,Encoding>::select_parent() const
{
	mtrandom mt;
	int n = this->m_population->size();
	Chromosome<Encoding>& candidate1 = (*(this->m_population))[mt.random(n)];
	Chromosome<Encoding>& candidate2 = (*(this->m_population))[mt.random(n)];

	if(m_comp->compare(candidate1,candidate2) == -1) {
		return candidate1;
	} else {
		return candidate2;
	}
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
ranking_selection<Chromosome,Encoding>::ranking_selection() :
	m_bias(1.5)
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
ranking_selection<Chromosome,Encoding>::~ranking_selection()
{
}

/*!
 * \brief initialize the linear bias
 */
template <template <typename> class Chromosome, typename Encoding>
void ranking_selection<Chromosome,Encoding>::initialize()
{
	kvparse::parameter_value(keywords::RANKING_BIAS, m_bias, false);
}

/*!
 * \brief select a parent from the gene pool
 */
template <template <typename> class Chromosome, typename Encoding>
Chromosome<Encoding>& ranking_selection<Chromosome,Encoding>::select_parent() const
{
	mtrandom mt;

	int index = static_cast<int>(this->m_population->size() *
	                             (m_bias - sqrt(m_bias * m_bias - 4.0 * (m_bias - 1) * mt.random()))
	                             / 2.0 / (m_bias - 1));
	return (*(this->m_population))[index];
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
random_selection<Chromosome,Encoding>::random_selection()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
random_selection<Chromosome,Encoding>::~random_selection()
{
}

/*!
 * \brief select a parent at random from the gene pool
 */
template <template <typename> class Chromosome, typename Encoding>
Chromosome<Encoding>& random_selection<Chromosome,Encoding>::select_parent() const
{
	mtrandom mt;
	return (*(this->m_population))[mt.random(this->m_population->size())];
}

/*!
 * \brief create and return an initialized selection scheme
 */
template <template <typename> class Chromosome, typename Encoding>
selection_scheme<Chromosome,Encoding>* selection_scheme_factory<Chromosome,Encoding>::construct()
{
	string ssname;
	kvparse::parameter_value(keywords::SELECTION_SCHEME, ssname, true);

	if(ssname == keywords::TOURNAMENT_SELECTION) {
		tournament_selection<Chromosome,Encoding>* ts =
		    new tournament_selection<Chromosome,Encoding>;
		ts->initialize();
		return ts;
	} else if(ssname == keywords::RANKING_SELECTION) {
		ranking_selection<Chromosome,Encoding>* rs =
		    new ranking_selection<Chromosome,Encoding>;
		rs->initialize();
		return rs;
	} else if(ssname == keywords::RANDOM_SELECTION) {
		return new random_selection<Chromosome,Encoding>;
	} else {
		error("invalid selection_scheme specified: " + ssname);
		return 0;
	}
}
