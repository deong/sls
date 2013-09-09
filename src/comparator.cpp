/*!
 * \file comparator.cpp
 *
 * defines classes which determine how chromosomes are compared to
 * one another by fitness.
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <vector>
#include "comparator.h"
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
comparator<Chromosome,Encoding>::comparator()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
comparator<Chromosome,Encoding>::~comparator()
{
}

/*!
 * \brief functor to allow dispatch to virtual compare method
 */
template <template <typename> class Chromosome, typename Encoding>
inline bool comparator<Chromosome,Encoding>::operator()(const Chromosome<Encoding>& c1,
        const Chromosome<Encoding>& c2) const
{
	return compare(c1,c2) == -1;
}

/*!
 * \brief initialize the comparator
 *
 * empty definition provided so that derived classes need only provide
 * the method if they require initialization
 */
template <template <typename> class Chromosome, typename Encoding>
void comparator<Chromosome,Encoding>::initialize()
{
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
fitness_comparator<Chromosome,Encoding>::fitness_comparator()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
fitness_comparator<Chromosome,Encoding>::~fitness_comparator()
{
}

/*!
 * \brief compare based on first objective function
 */
template <template <typename> class Chromosome, typename Encoding>
inline int fitness_comparator<Chromosome,Encoding>::compare(const Chromosome<Encoding>& c1,
        const Chromosome<Encoding>& c2) const
{
	if(c1.fitness[0] < c2.fitness[0]) {
		return -1;
	} else if(c1.fitness[0] > c2.fitness[0]) {
		return 1;
	} else {
		return 0;
	}
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
pareto_dominance_comparator<Chromosome,Encoding>::pareto_dominance_comparator()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
pareto_dominance_comparator<Chromosome,Encoding>::~pareto_dominance_comparator()
{
}

/*!
 * \brief compare based on pareto dominance
 */
template <template <typename> class Chromosome, typename Encoding>
inline int pareto_dominance_comparator<Chromosome,Encoding>::compare(const Chromosome<Encoding>& c1,
        const Chromosome<Encoding>& c2) const
{
	int res = 0;
	for(unsigned int i=0; i<c1.fitness.size(); i++) {
		if(c1.fitness[i] > c2.fitness[i]) {
			if(res == -1) {
				return 0;
			} else {
				res = 1;
			}
		} else if(c1.fitness[i] < c2.fitness[i]) {
			if(res == 1) {
				return 0;
			} else {
				res = -1;
			}
		}
	}
	return res;
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
scalarizing_comparator<Chromosome,Encoding>::scalarizing_comparator()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
scalarizing_comparator<Chromosome,Encoding>::~scalarizing_comparator()
{
}

/*!
 * \brief read in the weight vector used to scalarize objectives
 */
template <template <typename> class Chromosome, typename Encoding>
void scalarizing_comparator<Chromosome,Encoding>::initialize(const string& prefix)
{
	string keyword=prefix+keywords::WEIGHT_VECTOR;
	if(configuration::keyword_exists(keyword)) {
		configuration::vector_parameter<double>(keyword,weights,true);
	}
}

/*!
 * \brief initialize with default (empty) prefix
 */
template <template <typename> class Chromosome, typename Encoding>
void scalarizing_comparator<Chromosome,Encoding>::initialize()
{
	this->initialize("");
}

/*!
 * \brief create a random weight vector
 */
template <template <typename> class Chromosome, typename Encoding>
void scalarizing_comparator<Chromosome,Encoding>::randomize_weights(unsigned int nobj)
{
	weights.clear();
	weights.resize(nobj);

	mtrandom mt;
	double remainder = 1.0;
	for(unsigned int i=0; i<nobj-1; i++) {
		weights[i] = mt.random(0.0, remainder);
		remainder -= weights[i];
	}
	weights[nobj-1] = remainder;
}

/*!
 * \brief compare based on current scalarization of objective function values
 */
template <template <typename> class Chromosome, typename Encoding>
inline int scalarizing_comparator<Chromosome,Encoding>::compare(const Chromosome<Encoding>& c1,
        const Chromosome<Encoding>& c2) const
{
	double val1 = 0;
	double val2 = 0;

	for(unsigned int i=0; i<c1.fitness.size(); i++) {
		val1 += weights[i] * c1.fitness[i];
		val2 += weights[i] * c2.fitness[i];
	}

	if(val1 < val2) {
		return -1;
	} else if (val1 > val2) {
		return 1;
	} else {
		return 0;
	}
}

/*!
 * \brief compute the distance in objective space between two individuals
 */
template <template <typename> class Chromosome, typename Encoding>
inline double scalarizing_comparator<Chromosome,Encoding>::difference(const Chromosome<Encoding>& c1,
        const Chromosome<Encoding>& c2) const
{
	double val1=0;
	double val2=0;

	for(unsigned int i=0; i<c1.fitness.size(); i++) {
		val1+=weights[i]*c1.fitness[i];
		val2+=weights[i]*c2.fitness[i];
	}

	return val1-val2;
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
weak_dominance_comparator<Chromosome,Encoding>::weak_dominance_comparator()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
weak_dominance_comparator<Chromosome,Encoding>::~weak_dominance_comparator()
{
}

/*!
 * \brief compare based on weak dominance
 */
template <template <typename> class Chromosome, typename Encoding>
inline int weak_dominance_comparator<Chromosome,Encoding>::compare(const Chromosome<Encoding>& c1,
        const Chromosome<Encoding>& c2) const
{
	// look at the first element of the fitness vector to determine
	// which state we're in (see the documentation for a description
	// of the finite state machine that implements weak dominance)
	// note that weak dominance is ambiguous in that if two vectors
	// are identical then each weakly dominates the other.  In that
	// case, I return the -1 (indicating the first weakly dominates
	// the second)
	int state = 0;
	if(c1.fitness[0] <= c2.fitness[0]) {
		state = -1;
	} else {
		state = 1;
	}

	for(unsigned int i=1; i<c1.fitness.size(); i++) {
		if(state == 1) {
			if(c1.fitness[i] < c2.fitness[i]) {
				return 0;
			}
		} else {
			if(c1.fitness[i] > c2.fitness[i]) {
				return 0;
			}
		}
	}
	return state;
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
strong_dominance_comparator<Chromosome,Encoding>::strong_dominance_comparator()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
strong_dominance_comparator<Chromosome,Encoding>::~strong_dominance_comparator()
{
}

/*!
 * \brief compare based on strong dominance
 */
template <template <typename> class Chromosome, typename Encoding>
inline int strong_dominance_comparator<Chromosome,Encoding>::compare(const Chromosome<Encoding>& c1,
        const Chromosome<Encoding>& c2) const
{
	// look at the first element of the fitness vector to determine
	// which state we're in (see the documentation for a description
	// of the finite state machine that implements strong dominance)
	int state = 0;
	if(c1.fitness[0] == c2.fitness[0]) {
		return 0;
	} else if(c1.fitness[0] < c2.fitness[0]) {
		state = -1;
	} else {
		state = 1;
	}

	for(unsigned int i=1; i<c1.fitness.size(); i++) {
		if(c1.fitness[i] == c2.fitness[i]) {
			return 0;
		} else if(c1.fitness[i] < c2.fitness[i]) {
			if(state == 1) {
				return 0;
			}
		} else if(c1.fitness[i] > c2.fitness[i]) {
			if(state == -1) {
				return 0;
			}
		}
	}
	return state;
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
epsilon_dominance_comparator<Chromosome,Encoding>::epsilon_dominance_comparator()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
epsilon_dominance_comparator<Chromosome,Encoding>::~epsilon_dominance_comparator()
{
}

/*!
 * \brief read in the epsilon vector
 *
 * this defines the granularity of the imposed mesh
 */
template <template <typename> class Chromosome, typename Encoding>
void epsilon_dominance_comparator<Chromosome,Encoding>::initialize(const string& prefix)
{
	configuration::vector_parameter<double>(prefix+keywords::EPSILON, m_epsilon, true);
}

/*!
 * \brief initialize with default (empty) prefix
 */
template <template <typename> class Chromosome, typename Encoding>
void epsilon_dominance_comparator<Chromosome,Encoding>::initialize()
{
	this->initialize("");
}

/*!
 * \brief compare chromosomes based on epsilon dominance
 *
 * A <<eps B iff A is in a better "box" than B or they are in the same box,
 * but A is closer to the good corner.
 */
template <template <typename> class Chromosome, typename Encoding>
inline int epsilon_dominance_comparator<Chromosome,Encoding>::compare(const Chromosome<Encoding>& c1,
        const Chromosome<Encoding>& c2) const
{
	// epsilon dominance is equivalent to weak dominance except instead
	// of c1[i] <= c2[i] for all i, we need (1+epsilon)*c1[i] <= c2[i] for
	// all i
	int state = 0;
	if((1+m_epsilon[0])*c1.fitness[0] <= c2.fitness[0]) {
		state = -1;
	} else {
		state = 1;
	}

	for(unsigned int i=1; i<c1.fitness.size(); i++) {
		if(state == 1) {
			if((1+m_epsilon[i])*c1.fitness[i] < c2.fitness[i]) {
				return 0;
			}
		} else {
			if((1+m_epsilon[i])*c1.fitness[i] > c2.fitness[i]) {
				return 0;
			}
		}
	}
	return state;
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
single_objective_comparator<Chromosome,Encoding>::single_objective_comparator(unsigned int objnum)
{
	m_obj = objnum;
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
single_objective_comparator<Chromosome,Encoding>::~single_objective_comparator()
{
}

/*!
 * \brief compare chromosomes based on selected objective function value
 */
template <template <typename> class Chromosome, typename Encoding>
inline int single_objective_comparator<Chromosome,Encoding>::compare(const Chromosome<Encoding>& c1,
        const Chromosome<Encoding>& c2) const
{
	if(c1.fitness[m_obj] < c2.fitness[m_obj]) {
		return -1;
	} else if(c1.fitness[m_obj] > c2.fitness[m_obj]) {
		return 1;
	} else {
		return 0;
	}
}

/*!
 * \brief create and return an initialized comparator
 */
template <template <typename> class Chromosome, typename Encoding>
comparator<Chromosome,Encoding>* comparator_factory<Chromosome,Encoding>::construct()
{
	string comp;
	configuration::string_parameter(this->m_prefix+keywords::COMPARATOR, comp, true);
	if(comp == keywords::FITNESS_COMPARATOR) {
		return new fitness_comparator<Chromosome,Encoding>;
	} else if(comp == keywords::PARETO_DOMINANCE_COMPARATOR) {
		return new pareto_dominance_comparator<Chromosome,Encoding>;
	} else if(comp == keywords::WEAK_DOMINANCE_COMPARATOR) {
		return new weak_dominance_comparator<Chromosome,Encoding>;
	} else if(comp == keywords::STRONG_DOMINANCE_COMPARATOR) {
		return new strong_dominance_comparator<Chromosome,Encoding>;
	} else if(comp == keywords::EPSILON_DOMINANCE_COMPARATOR) {
		epsilon_dominance_comparator<Chromosome,Encoding>* c =
		    new epsilon_dominance_comparator<Chromosome,Encoding>;
		c->initialize(this->m_prefix);
		return c;
	} else if(comp == keywords::SCALARIZING_COMPARATOR) {
		scalarizing_comparator<Chromosome,Encoding>* c =
		    new scalarizing_comparator<Chromosome,Encoding>;
		c->initialize(this->m_prefix);
		return c;
	} else {
		cerr << "invalid comparator: " << comp << " specified" << endl;
		exit(1);
		return 0;
	}
}
