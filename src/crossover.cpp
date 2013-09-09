/*!
 * \file crossover.cpp
 *
 * evolutionary algorithm crossover operators
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <string>
#include <algorithm>
#include "crossover.h"
#include "encoding.h"
#include "mtrandom.h"
#include "chromosome.h"
#include "configuration.h"
#include "keywords.h"
#include "utilities.h"
#include "mtrandom.h"

using namespace std;

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
crossover_operator<Chromosome,Encoding>::crossover_operator()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
crossover_operator<Chromosome,Encoding>::~crossover_operator()
{
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
uniform_crossover_impl<Chromosome,Encoding>::uniform_crossover_impl()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
uniform_crossover_impl<Chromosome,Encoding>::~uniform_crossover_impl()
{
}

/*!
 * \brief perform uniform crossover
 */
template <template <typename> class Chromosome, typename Encoding>
void uniform_crossover_impl<Chromosome,Encoding>::crossover(const Chromosome<Encoding>& p1,
        const Chromosome<Encoding>& p2,
        Chromosome<Encoding>& c1,
        Chromosome<Encoding>& c2) const
{
	mtrandom mt;
	for(unsigned int i=0; i<p1.length(); i++) {
		if(mt.random() < 0.5) {
			c1[i] = p1[i];
			c2[i] = p2[i];
		} else {
			c1[i] = p2[i];
			c2[i] = p1[i];
		}
	}
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
one_point_crossover_impl<Chromosome,Encoding>::one_point_crossover_impl()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
one_point_crossover_impl<Chromosome,Encoding>::~one_point_crossover_impl()
{
}

/*!
 * \brief perform one-point crossover
 */
template <template <typename> class Chromosome, typename Encoding>
void one_point_crossover_impl<Chromosome,Encoding>::crossover(const Chromosome<Encoding>& p1,
        const Chromosome<Encoding>& p2,
        Chromosome<Encoding>& c1,
        Chromosome<Encoding>& c2) const
{
	mtrandom mt;
	int xpt = mt.random(p1.length());
	for(int i=0; i<xpt; i++) {
		c1[i] = p1[i];
		c2[i] = p2[i];
	}
	for(int i=xpt; i<int(p1.length()); i++) {
		c1[i] = p2[i];
		c2[i] = p1[i];
	}
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
two_point_crossover_impl<Chromosome,Encoding>::two_point_crossover_impl()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
two_point_crossover_impl<Chromosome,Encoding>::~two_point_crossover_impl()
{
}

/*!
 * \brief perform two-point crossover
 */
template <template <typename> class Chromosome, typename Encoding>
void two_point_crossover_impl<Chromosome,Encoding>::crossover(const Chromosome<Encoding>& p1,
        const Chromosome<Encoding>& p2,
        Chromosome<Encoding>& c1,
        Chromosome<Encoding>& c2) const
{
	mtrandom mt;
	int xpt1;
	int xpt2;

	xpt1 = mt.random(p1.length()) + 1;
	do {
		xpt2 = mt.random(p1.length()-1);
	} while (xpt1 == xpt2);

	if(xpt1 > xpt2) {
		int tmp = xpt1;
		xpt1 = xpt2;
		xpt2 = tmp;
	}

	for(int i=0; i<xpt1; i++) {
		c1[i] = p1[i];
		c2[i] = p2[i];
	}
	for(int i=xpt1; i<xpt2; i++) {
		c1[i] = p2[i];
		c2[i] = p1[i];
	}
	for(int i=xpt2; i<int(p1.length()); i++) {
		c1[i] = p1[i];
		c2[i] = p2[i];
	}
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
hux_crossover_impl<Chromosome,Encoding>::hux_crossover_impl()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
hux_crossover_impl<Chromosome,Encoding>::~hux_crossover_impl()
{
}

/*!
 * \brief perform HUX crossover
 */
template <template <typename> class Chromosome, typename Encoding>
void hux_crossover_impl<Chromosome,Encoding>::crossover(const Chromosome<Encoding>& p1,
        const Chromosome<Encoding>& p2,
        Chromosome<Encoding>& c1,
        Chromosome<Encoding>& c2) const
{
	mtrandom mt;

	vector<int> diff;
	for(unsigned int i=0; i<p1.length(); i++) {
		if(p1[i] != p2[i]) {
			diff.push_back(i);
		}
	}

	unsigned int d = diff.size()/2;
	mt.shuffle(diff);

	c1 = p1;
	c2 = p2;
	for(unsigned int i=0; i<d; i++) {
		c1[diff[i]] = p2[diff[i]];
		c2[diff[i]] = p1[diff[i]];
	}
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
sbx_crossover_impl<Chromosome,Encoding>::sbx_crossover_impl() :
	m_eta(0)
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
sbx_crossover_impl<Chromosome,Encoding>::~sbx_crossover_impl()
{
}

/*!
 * \brief initialize the eta parameter of the SBX crossover operator
 */
template <template <typename> class Chromosome, typename Encoding>
void sbx_crossover_impl<Chromosome,Encoding>::initialize()
{
	configuration::double_parameter(keywords::SBX_ETA, m_eta, true);
}

/*!
 * \brief perform SBX crossover
 */
template <template <typename> class Chromosome, typename Encoding>
void sbx_crossover_impl<Chromosome,Encoding>::crossover(const Chromosome<Encoding>& p1,
        const Chromosome<Encoding>& p2,
        Chromosome<Encoding>& c1,
        Chromosome<Encoding>& c2) const
{
	mtrandom mt;
	double y1, y2, yl, yu;
	double ci1, ci2;
	double alpha, beta, betaq;

	for (unsigned int i=0; i<p1.length(); i++) {
		if (fabs(p1[i]-p2[i]) > 1e-10) {
			if (p1[i] < p2[i]) {
				y1 = p1[i];
				y2 = p2[i];
			} else {
				y1 = p2[i];
				y2 = p1[i];
			}

			pair<double,double> range = Encoding::parameter_range(i);
			yl = range.first;
			yu = range.second;

			double rn = mt.random();

			beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
			alpha = 2.0 - pow(beta,-(m_eta+1.0));
			if (rn <= (1.0/alpha)) {
				betaq = pow ((rn*alpha),(1.0/(m_eta+1.0)));
			} else {
				betaq = pow ((1.0/(2.0 - rn*alpha)),(1.0/(m_eta+1.0)));
			}
			ci1 = 0.5*((y1+y2)-betaq*(y2-y1));
			beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
			alpha = 2.0 - pow(beta,-(m_eta+1.0));
			if (rn <= (1.0/alpha)) {
				betaq = pow ((rn*alpha),(1.0/(m_eta+1.0)));
			} else {
				betaq = pow ((1.0/(2.0 - rn*alpha)),(1.0/(m_eta+1.0)));
			}
			ci2 = 0.5*((y1+y2)+betaq*(y2-y1));

			// make sure the variables are in bounds
			ci1 = max(ci1,yl);
			ci2 = max(ci2,yl);
			ci1 = min(ci1,yu);
			ci2 = min(ci2,yu);

			if (mt.random() < 0.5) {
				c1[i] = ci2;
				c2[i] = ci1;
			} else {
				c1[i] = ci1;
				c2[i] = ci2;
			}
		} else {
			c1[i] = p1[i];
			c2[i] = p2[i];
		}
	}
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
cycle_crossover_impl<Chromosome,Encoding>::cycle_crossover_impl()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
cycle_crossover_impl<Chromosome,Encoding>::~cycle_crossover_impl()
{
}

/*!
 * \brief perform cycle crossover
 */
template <template <typename> class Chromosome, typename Encoding>
void cycle_crossover_impl<Chromosome,Encoding>::crossover(const Chromosome<Encoding>& p1,
        const Chromosome<Encoding>& p2,
        Chromosome<Encoding>& c1,
        Chromosome<Encoding>& c2) const
{
	mtrandom mt;

	// find a cycle starting from a random point
	int n = (int)p1.length();
	int cycle_len = 0;
	vector<int> perm = mt.permutation(n);
	int index;

	bool mix = true;            // mix cycles from different parents
	while(cycle_len < n) {
		// pick the first point in the cycle
		for(index=0; index<n; index++) {
			if(perm[index] != -1) {
				break;
			}
		}
		index = perm[index];

		// store the points in the cycle
		int start = p1[index];
		vector<int> cycle;
		do {
			cycle.push_back(index);
			perm[index] = -1;   // this point has been included in a cycle

			// set index equal to position of p2[index] in p1
			for(int i=0; i<(int)p2.length(); i++) {
				if(p1[i] == p2[index]) {
					index = i;
					break;
				}
			}
			cycle_len++;
		} while(p1[index] != start);

		// now that we have a cycle, update the offspring
		for(unsigned int i=0; i<cycle.size(); i++) {
			if(mix) {
				c1[cycle[i]] = p2[cycle[i]];
				c2[cycle[i]] = p1[cycle[i]];
			} else {
				c1[cycle[i]] = p1[cycle[i]];
				c2[cycle[i]] = p2[cycle[i]];
			}
		}
		mix = !mix;
	}
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
order_crossover_impl<Chromosome,Encoding>::order_crossover_impl()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
order_crossover_impl<Chromosome,Encoding>::~order_crossover_impl()
{
}

/*!
 * \brief perform order crossover
 */
template <template <typename> class Chromosome, typename Encoding>
void order_crossover_impl<Chromosome,Encoding>::crossover(const Chromosome<Encoding>& p1,
        const Chromosome<Encoding>& p2,
        Chromosome<Encoding>& c1,
        Chromosome<Encoding>& c2) const
{
	mtrandom mt;
	int n = (int)p1.length();

	int xpt1 = mt.random(1,n-1);
	int xpt2 = mt.random(1,n-1);
	while(xpt1 == xpt2) {
		xpt2 = mt.random(1,n-1);
	}

	if(xpt1 > xpt2) {
		swap(xpt1, xpt2);
	}

	// store the previously used alleles in a vector
	vector<int> used_in_c1(n);
	vector<int> used_in_c2(n);

	// copy the middle sections directly into the offspring
	for(int i=xpt1; i<xpt2; i++) {
		c1[i] = p1[i];
		c2[i] = p2[i];
		used_in_c1.push_back(p1[i]);
		used_in_c2.push_back(p2[i]);
	}

	// now starting at the second crossover point, use the order
	// genes in the second parent to determine the order in the
	// first offspring
	int index = xpt2;
	vector<int>::iterator start, end, iter;
	for(int i=xpt2; i<n; i++) {
		iter = find(used_in_c1.begin(), used_in_c1.end(), p2[i]);
		if(iter == end) {
			c1[index++] = p2[i];
			used_in_c1.push_back(p2[i]);
			if(index == n) {
				index = 0;
			}
		}
	}
	for(int i=0; i<xpt1; i++) {
		iter = find(used_in_c1.begin(), used_in_c1.end(), p2[i]);
		if(iter == end) {
			c1[index++] = p2[i];
			used_in_c1.push_back(p2[i]);
			if(index == n) {
				index = 0;
			}
		}
	}

	// now repeat the process to create the second offspring
	index = xpt2;
	for(int i=xpt2; i<n; i++) {
		iter = find(used_in_c2.begin(), used_in_c2.end(), p1[i]);
		if(iter == end) {
			c2[index++] = p1[i];
			used_in_c2.push_back(p1[i]);
			if(index == n) {
				index = 0;
			}
		}
	}
	for(int i=0; i<xpt2; i++) {
		iter = find(used_in_c2.begin(), used_in_c2.end(), p1[i]);
		if(iter == end) {
			c2[index++] = p1[i];
			used_in_c2.push_back(p1[i]);
			if(index == n) {
				index = 0;
			}
		}
	}
}

/*!
 * \brief create an initialized crossover operator for bit-vector encodings
 */
template <template <typename> class Chromosome, typename Encoding>
crossover_operator<Chromosome,Encoding>* bit_vector_crossover_operator_factory<Chromosome,Encoding>::construct()
{
	string coname;
	configuration::string_parameter(keywords::CROSSOVER_OPERATOR, coname, true);
	if(coname == keywords::UNIFORM_CROSSOVER) {
		return new uniform_crossover<Chromosome,Encoding>;
	} else if(coname == keywords::ONE_POINT_CROSSOVER) {
		return new one_point_crossover<Chromosome,Encoding>;
	} else if(coname == keywords::TWO_POINT_CROSSOVER) {
		return new two_point_crossover<Chromosome,Encoding>;
	} else if(coname == keywords::HUX_CROSSOVER) {
		return new hux_crossover<Chromosome,Encoding>;
	} else {
		error("illegal crossover_operator specified: " + coname);
		return 0;
	}
}

/*!
 * \brief create an initialized crossover operator for real encodings
 */
template <template <typename> class Chromosome, typename Encoding>
crossover_operator<Chromosome,Encoding>* real_crossover_operator_factory<Chromosome,Encoding>::construct()
{
	string coname;
	configuration::string_parameter(keywords::CROSSOVER_OPERATOR, coname, true);
	if(coname == keywords::UNIFORM_CROSSOVER) {
		return new uniform_crossover<Chromosome,Encoding>;
	} else if(coname == keywords::ONE_POINT_CROSSOVER) {
		return new one_point_crossover<Chromosome,Encoding>;
	} else if(coname == keywords::TWO_POINT_CROSSOVER) {
		return new two_point_crossover<Chromosome,Encoding>;
	} else if(coname == keywords::SBX_CROSSOVER) {
		sbx_crossover<Chromosome,Encoding>* sbx = new sbx_crossover<Chromosome,Encoding>;
		sbx->initialize();
		return sbx;
	} else {
		error("illegal crossover_operator specified: " + coname);
		return 0;
	}
}

/*!
 * \brief create an initialized crossover operator for permutation encodings
 */
template <template <typename> class Chromosome, typename Encoding>
crossover_operator<Chromosome,Encoding>* permutation_crossover_operator_factory<Chromosome,Encoding>::construct()
{
	string coname;
	configuration::string_parameter(keywords::CROSSOVER_OPERATOR, coname, true);
	if(coname == keywords::ORDER_CROSSOVER) {
		return new order_crossover<Chromosome,Encoding>;
	} else if(coname == keywords::CYCLE_CROSSOVER) {
		return new cycle_crossover<Chromosome,Encoding>;
	} else {
		error("illegal crossover_operator specified: " + coname);
		return 0;
	}
}

/*!
 * \brief create an initialized crossover operator for integer encodings
 */
template <template <typename> class Chromosome, typename Encoding>
crossover_operator<Chromosome,Encoding>* integer_crossover_operator_factory<Chromosome,Encoding>::construct()
{
	string coname;
	configuration::string_parameter(keywords::CROSSOVER_OPERATOR, coname, true);
	if(coname == keywords::UNIFORM_CROSSOVER) {
		return new uniform_crossover<Chromosome,Encoding>;
	} else if(coname == keywords::ONE_POINT_CROSSOVER) {
		return new one_point_crossover<Chromosome,Encoding>;
	} else if(coname == keywords::TWO_POINT_CROSSOVER) {
		return new two_point_crossover<Chromosome,Encoding>;
	} else {
		error("illegal crossover_operator specified: " + coname);
		return 0;
	}
}

/*!
 * \brief create an initialized crossover operator for gap encodings
 */
template <template <typename> class Chromosome, typename Encoding>
crossover_operator<Chromosome,Encoding>* gap_crossover_operator_factory<Chromosome,Encoding>::construct()
{
	string coname;
	configuration::string_parameter(keywords::CROSSOVER_OPERATOR, coname, true);
	if(coname == keywords::UNIFORM_CROSSOVER) {
		return new uniform_crossover<Chromosome,Encoding>;
	} else if(coname == keywords::ONE_POINT_CROSSOVER) {
		return new one_point_crossover<Chromosome,Encoding>;
	} else if(coname == keywords::TWO_POINT_CROSSOVER) {
		return new two_point_crossover<Chromosome,Encoding>;
	} else {
		error("illegal crossover_operator specified: " + coname);
		return 0;
	}
}

/*!
 * \brief create an initialized crossover operator for gsap encodings
 */
template <template <typename> class Chromosome, typename Encoding>
crossover_operator<Chromosome,Encoding>* gsap_crossover_operator_factory<Chromosome,Encoding>::construct()
{
	string coname;
	configuration::string_parameter(keywords::CROSSOVER_OPERATOR, coname, true);
	if(coname == keywords::UNIFORM_CROSSOVER) {
		return new uniform_crossover<Chromosome,Encoding>;
	} else if(coname == keywords::ONE_POINT_CROSSOVER) {
		return new one_point_crossover<Chromosome,Encoding>;
	} else if(coname == keywords::TWO_POINT_CROSSOVER) {
		return new two_point_crossover<Chromosome,Encoding>;
	} else {
		error("illegal crossover_operator specified: " + coname);
		return 0;
	}
}

