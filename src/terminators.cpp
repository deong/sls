/*!
 * \file terminators.cpp
 *
 * in stochastic local search, there is typically no one method by
 * which we can decide when to terminate the search.  these classes
 * provide several alternative mechanism by which the search process
 * can be halted.
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <list>
#include <string>
#include "terminators.h"
#include "chromosome.h"
#include "encoding.h"
#include "population.h"
#include "configuration.h"
#include "keywords.h"
#include "utilities.h"

using namespace std;

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
terminator<Chromosome,Encoding>::terminator()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
terminator<Chromosome,Encoding>::~terminator()
{
}

/*!
 * \brief initialize the terminator components
 */
template <template <typename> class Chromosome, typename Encoding>
void terminator<Chromosome,Encoding>::initialize(const string& prefix)
{
	m_prefix=prefix;
}

/*!
 * \brief notify of chromosome evaluation
 */
template <template <typename> class Chromosome, typename Encoding>
void terminator<Chromosome,Encoding>::chromosome_evaluated(const Chromosome<Encoding>& sol)
{
}

/*!
 * \brief notify of generation completion
 */
template <template <typename> class Chromosome, typename Encoding>
void terminator<Chromosome,Encoding>::generation_completed()
{
}

/*!
 * \brief notify of generation completion
 */
template <template <typename> class Chromosome, typename Encoding>
void terminator<Chromosome,Encoding>::generation_completed(const population<Chromosome,Encoding>& pop)
{
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
evaluation_limit<Chromosome,Encoding>::evaluation_limit() :
	m_evals(0),
	m_max_evals(0)
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
evaluation_limit<Chromosome,Encoding>::~evaluation_limit()
{
}

/*!
 * \brief initialize the maximum number of allowed evaluations
 */
template <template <typename> class Chromosome, typename Encoding>
void evaluation_limit<Chromosome,Encoding>::initialize(const string& prefix)
{
	terminator<Chromosome,Encoding>::initialize(prefix);
	configuration::unsigned_integer_parameter(prefix+keywords::MAX_EVALUATIONS, m_max_evals, true);
}

/*!
 * \brief reset the eval counter
 */
template <template <typename> class Chromosome, typename Encoding>
void evaluation_limit<Chromosome,Encoding>::reset()
{
	m_evals=0;
}

/*!
 * \brief increment the evaluation count
 */
template <template <typename> class Chromosome, typename Encoding>
void evaluation_limit<Chromosome,Encoding>::chromosome_evaluated(const Chromosome<Encoding>& sol)
{
	m_evals++;
}

/*!
 * \brief return true if the evaluation count exceeds the maximum allowed
 */
template <template <typename> class Chromosome, typename Encoding>
bool evaluation_limit<Chromosome,Encoding>::terminate() const
{
	return m_evals >= m_max_evals;
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
generation_limit<Chromosome,Encoding>::generation_limit() :
	m_gens(0),
	m_max_gens(0)
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
generation_limit<Chromosome,Encoding>::~generation_limit()
{
}

/*!
 * \brief initialize the maximum allowable number of generations
 */
template <template <typename> class Chromosome, typename Encoding>
void generation_limit<Chromosome,Encoding>::initialize(const string& prefix)
{
	configuration::unsigned_integer_parameter(prefix+keywords::MAX_GENERATIONS,
	        m_max_gens, true);
}

/*!
 * \brief reset the generation count
 */
template <template <typename> class Chromosome, typename Encoding>
void generation_limit<Chromosome,Encoding>::reset()
{
	m_gens=0;
}

/*!
 * \brief increment the generation count
 */
template <template <typename> class Chromosome, typename Encoding>
void generation_limit<Chromosome,Encoding>::generation_completed()
{
	m_gens++;
}

/*!
 * \brief increment the generation count
 */
template <template <typename> class Chromosome, typename Encoding>
void generation_limit<Chromosome,Encoding>::generation_completed(const population<Chromosome,Encoding>& pop)
{
	m_gens++;
}

/*!
 * \brief return true if the generation count exceeds the maximum allowed
 */
template <template <typename> class Chromosome, typename Encoding>
bool generation_limit<Chromosome,Encoding>::terminate() const
{
	return m_gens >= m_max_gens;
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
null_terminator<Chromosome,Encoding>::null_terminator()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
null_terminator<Chromosome,Encoding>::~null_terminator()
{
}

/*!
 * \brief reset of null terminator does nothing
 */
template <template <typename> class Chromosome, typename Encoding>
void null_terminator<Chromosome,Encoding>::reset()
{
}

/*!
 * \brief always return false
 *
 * useful for algorithms with implicitly defined termination conditions
 */
template <template <typename> class Chromosome, typename Encoding>
bool null_terminator<Chromosome,Encoding>::terminate() const
{
	return false;
}

/*!
 * \brief return a list of initialized terminators
 */
template <template <typename> class Chromosome, typename Encoding>
list<terminator<Chromosome,Encoding>*> terminator_factory<Chromosome,Encoding>::construct()
{
	list<terminator<Chromosome,Encoding>*> termlist;

	list<string> tnames;
	configuration::list_parameter(this->m_prefix+keywords::TERMINATOR, tnames, true);

	for(list<string>::iterator it=tnames.begin(); it!=tnames.end(); it++) {
		string tname = (*it);
		if(tname == keywords::EVALUATION_LIMIT) {
			evaluation_limit<Chromosome,Encoding>* term = new evaluation_limit<Chromosome,Encoding>;
			term->initialize(this->m_prefix);
			termlist.push_back(term);
		} else if(tname == keywords::GENERATION_LIMIT) {
			generation_limit<Chromosome,Encoding>* term = new generation_limit<Chromosome,Encoding>;
			term->initialize(this->m_prefix);
			termlist.push_back(term);
		} else if(tname == keywords::NULL_TERMINATOR) {
			null_terminator<Chromosome,Encoding>* term = new null_terminator<Chromosome,Encoding>;
			termlist.push_back(term);
		} else {
			error("illegal terminator specified: " + tname);
		}
	}
	return termlist;
}
