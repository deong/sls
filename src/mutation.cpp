/*!
 * \file mutation.cpp
 *
 * evolutionary algorithm mutation operators
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <string>
#include <algorithm>
#include <cassert>
#include "mutation.h"
#include "chromosome.h"
#include "encoding.h"
#include "mtrandom.h"
#include "configuration.h"
#include "keywords.h"
#include "utilities.h"

using namespace std;

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
mutation_operator<Chromosome,Encoding>::mutation_operator()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
mutation_operator<Chromosome,Encoding>::~mutation_operator()
{
}

/*!
 * \brief initialize mutation operator parameters
 *
 * empty method provided here so that subclasses need only implement
 * the method if they require non-standard initialization
 */
template <template <typename> class Chromosome, typename Encoding>
void mutation_operator<Chromosome,Encoding>::initialize()
{
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
bitwise_mutation_impl<Chromosome,Encoding>::bitwise_mutation_impl()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
bitwise_mutation_impl<Chromosome,Encoding>::~bitwise_mutation_impl()
{
}

/*!
 * \brief initialize the mutation rate
 */
template <template <typename> class Chromosome, typename Encoding>
void bitwise_mutation_impl<Chromosome,Encoding>::initialize()
{
	configuration::double_parameter(keywords::MUTATION_RATE, m_rate, true);
}

/*!
 * \brief flip individual bits with a given probability
 */
template <template <typename> class Chromosome, typename Encoding>
void bitwise_mutation_impl<Chromosome,Encoding>::mutate(Chromosome<Encoding>& sol) const
{
	mtrandom mt;

	for(unsigned int i=0; i<sol.length(); i++) {
		if(mt.random() < m_rate) {
			sol[i] = (sol[i] == 0) ? 1 : 0;
		}
	}
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
swap_mutation_impl<Chromosome,Encoding>::swap_mutation_impl() :
	m_rate(0.0)
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
swap_mutation_impl<Chromosome,Encoding>::~swap_mutation_impl()
{
}

/*!
 * \brief initialize the mutation rate
 */
template <template <typename> class Chromosome, typename Encoding>
void swap_mutation_impl<Chromosome,Encoding>::initialize()
{
	configuration::double_parameter(keywords::MUTATION_RATE, m_rate, true);
}

/*!
 * \brief swap alleles with a given probability
 */
template <template <typename> class Chromosome, typename Encoding>
void swap_mutation_impl<Chromosome,Encoding>::mutate(Chromosome<Encoding>& sol) const
{
	mtrandom mt;
	for(unsigned int i=0; i<sol.length(); i++) {
		if(mt.random() < m_rate) {
			int other = mt.random(sol.length());
			swap(sol[i], sol[other]);
		}
	}
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
gaussian_mutation_impl<Chromosome,Encoding>::gaussian_mutation_impl() :
	m_rate(0),
	m_mu(0.0),
	m_sigma(1.0)
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
gaussian_mutation_impl<Chromosome,Encoding>::~gaussian_mutation_impl()
{
}

/*!
 * \brief initialize the distribution parameters
 */
template <template <typename> class Chromosome, typename Encoding>
void gaussian_mutation_impl<Chromosome,Encoding>::initialize()
{
	configuration::double_parameter(keywords::MUTATION_RATE, m_rate, true);
	configuration::double_parameter(keywords::GAUSSIAN_MUTATION_MU, m_mu, true);
	configuration::double_parameter(keywords::GAUSSIAN_MUTATION_SIGMA, m_sigma, true);
}

/*!
 * \brief add a random gaussian to the parameters with a specified probability
 */
template <template <typename> class Chromosome, typename Encoding>
void gaussian_mutation_impl<Chromosome,Encoding>::mutate(Chromosome<Encoding>& sol) const
{
	mtrandom mt;
	for(unsigned int i=0; i<sol.length(); i++) {
		pair<double,double> range = Encoding::parameter_range(i);
		if(mt.random() < m_rate) {
			sol[i] = sol[i] + mt.gaussian(m_mu, m_sigma);
			sol[i] = max(min(sol[i],range.second),range.first);
		}
	}
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
polynomial_mutation_impl<Chromosome,Encoding>::polynomial_mutation_impl() :
	m_rate(0),
	m_eta(0)
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
polynomial_mutation_impl<Chromosome,Encoding>::~polynomial_mutation_impl()
{
}

/*!
 * \brief initialize the distribution parameters
 */
template <template <typename> class Chromosome, typename Encoding>
void polynomial_mutation_impl<Chromosome,Encoding>::initialize()
{
	configuration::double_parameter(keywords::MUTATION_RATE, m_rate, true);
	configuration::double_parameter(keywords::POLYNOMIAL_ETA, m_eta, true);
}

/*!
 * \brief add a polynomially generated random variate to the parameters
 */
template <template <typename> class Chromosome, typename Encoding>
void polynomial_mutation_impl<Chromosome,Encoding>::mutate(Chromosome<Encoding>& sol) const
{
	mtrandom mt;
	for(unsigned int j=0; j<sol.length(); j++) {
		if(mt.random() < m_rate) {
			pair<double,double> range = Encoding::parameter_range(j);
			double y = sol[j];
			double yl = range.first;
			double yu = range.second;
			double delta1 = (y-yl)/(yu-yl);
			double delta2 = (yu-y)/(yu-yl);
			double rnd = mt.random();
			double mut_pow = 1.0/(m_eta+1.0);
			double deltaq;
			if(rnd <= 0.5) {
				double xy = 1.0-delta1;
				double val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(m_eta+1.0)));
				deltaq =  pow(val,mut_pow) - 1.0;
			} else {
				double xy = 1.0-delta2;
				double val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(m_eta+1.0)));
				deltaq = 1.0 - (pow(val,mut_pow));
			}
			y = y + deltaq*(yu-yl);
			sol[j] = max(min(y,yu),yl);
		}
	}
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
shift_mutation_impl<Chromosome,Encoding>::shift_mutation_impl() :
	m_rate(0)
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
shift_mutation_impl<Chromosome,Encoding>::~shift_mutation_impl()
{
}

/*!
 * \brief initialize the mutation rate
 */
template <template <typename> class Chromosome, typename Encoding>
void shift_mutation_impl<Chromosome,Encoding>::initialize()
{
	configuration::double_parameter(keywords::MUTATION_RATE, m_rate, true);
}

/*!
 * \brief randomly reassign selected tasks
 */
template <template <typename> class Chromosome, typename Encoding>
void shift_mutation_impl<Chromosome,Encoding>::mutate(Chromosome<Encoding>& sol) const
{
	mtrandom mt;
	for(unsigned int i=0; i<sol.length(); i++) {
		if(mt.random() < m_rate) {
			sol[i] = mt.random(sol.agents());
		}
	}
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
sss_mutation_impl<Chromosome,Encoding>::sss_mutation_impl() :
	m_rate(0)
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
sss_mutation_impl<Chromosome,Encoding>::~sss_mutation_impl()
{
}

/*!
 * \brief initialize the mutation rate
 */
template <template <typename> class Chromosome, typename Encoding>
void sss_mutation_impl<Chromosome,Encoding>::initialize()
{
	configuration::double_parameter(keywords::MUTATION_RATE, m_rate, true);
}

/*!
 * \brief randomly reassign selected tasks or swap tasks
 */
template <template <typename> class Chromosome, typename Encoding>
void sss_mutation_impl<Chromosome,Encoding>::mutate(Chromosome<Encoding>& sol) const
{
	mtrandom mt;
	for(unsigned int i=0; i<sol.length(); i++) {
		if(mt.random() < 0.5) {
			if(mt.random() < m_rate) {
				sol[i] = mt.random(sol.agents());
			}
		} else {
			int other = mt.random(sol.length());
			swap(sol[i], sol[other]);
		}
	}
}

/*!
 * \brief construct mutation operators for binary encodings
 */
template <template <typename> class Chromosome, typename Encoding>
mutation_operator<Chromosome,Encoding>* bit_mutation_operator_factory<Chromosome,Encoding>::construct()
{
	string opname;
	configuration::string_parameter(keywords::MUTATION_OPERATOR, opname, true);

	if(opname == keywords::BITWISE_MUTATION) {
		bitwise_mutation<Chromosome,Encoding>* m = new bitwise_mutation<Chromosome,Encoding>;
		m->initialize();
		return m;
	} else {
		error("illegal mutation_operator specified: " + opname);
		return 0;
	}
}

/*!
 * \brief construct mutation operators for permutation encodings
 */
template <template <typename> class Chromosome, typename Encoding>
mutation_operator<Chromosome,Encoding>* permutation_mutation_operator_factory<Chromosome,Encoding>::construct()
{
	string opname;
	configuration::string_parameter(keywords::MUTATION_OPERATOR, opname, true);

	if(opname == keywords::SWAP_MUTATION) {
		swap_mutation<Chromosome,Encoding>* m = new swap_mutation<Chromosome,Encoding>;
		m->initialize();
		return m;
	} else {
		error("illegal mutation operator specified: " + opname);
		return 0;
	}
}

/*!
 * \brief construct mutation operators for real encodings
 */
template <template <typename> class Chromosome, typename Encoding>
mutation_operator<Chromosome,Encoding>* real_mutation_operator_factory<Chromosome,Encoding>::construct()
{
	string opname;
	configuration::string_parameter(keywords::MUTATION_OPERATOR, opname, true);

	if(opname == keywords::GAUSSIAN_MUTATION) {
		gaussian_mutation<Chromosome,Encoding>* m = new gaussian_mutation<Chromosome,Encoding>;
		m->initialize();
		return m;
	} else if(opname == keywords::POLYNOMIAL_MUTATION) {
		polynomial_mutation<Chromosome,Encoding>* m = new polynomial_mutation<Chromosome,Encoding>;
		m->initialize();
		return m;
	} else {
		error("illegal mutation operator specified: " + opname);
		return 0;
	}
}

/*!
 * \brief construct mutation operators for integer encodings
 */
template <template <typename> class Chromosome, typename Encoding>
mutation_operator<Chromosome,Encoding>* integer_mutation_operator_factory<Chromosome,Encoding>::construct()
{
	error("no mutation operator defined for integer encoding.");
	return 0;
}

/*!
 * \brief construct mutation operators for gap encodings
 */
template <template <typename> class Chromosome, typename Encoding>
mutation_operator<Chromosome,Encoding>* gap_mutation_operator_factory<Chromosome,Encoding>::construct()
{
	string opname;
	configuration::string_parameter(keywords::MUTATION_OPERATOR, opname, true);

	if(opname == keywords::SHIFT_MUTATION) {
		shift_mutation<Chromosome,Encoding>* m = new shift_mutation<Chromosome,Encoding>;
		m->initialize();
		return m;
	} else if(opname == keywords::SWAP_MUTATION) {
		swap_mutation<Chromosome,Encoding>* m = new swap_mutation<Chromosome,Encoding>;
		m->initialize();
		return m;
	} else if(opname==keywords::SSS_MUTATION) {
		sss_mutation<Chromosome,Encoding>* m=new sss_mutation<Chromosome,Encoding>;
		m->initialize();
		return m;
	} else {
		error("illegal mutation operator specified: " + opname);
		return 0;
	}
}

/*!
 * \brief construct mutation operators for gsap encodings
 */
template <template <typename> class Chromosome, typename Encoding>
mutation_operator<Chromosome,Encoding>* gsap_mutation_operator_factory<Chromosome,Encoding>::construct()
{
	string opname;
	configuration::string_parameter(keywords::MUTATION_OPERATOR, opname, true);

	if(opname == keywords::SHIFT_MUTATION) {
		shift_mutation<Chromosome,Encoding>* m = new shift_mutation<Chromosome,Encoding>;
		m->initialize();
		return m;
	} else if(opname == keywords::SWAP_MUTATION) {
		swap_mutation<Chromosome,Encoding>* m = new swap_mutation<Chromosome,Encoding>;
		m->initialize();
		return m;
	} else if(opname==keywords::SSS_MUTATION) {
		sss_mutation<Chromosome,Encoding>* m=new sss_mutation<Chromosome,Encoding>;
		m->initialize();
		return m;
	} else {
		error("illegal mutation operator specified: " + opname);
		return 0;
	}
}
