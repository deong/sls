/*!
 * \file neighborhood.cpp
 *
 * Defines neighborhood operators for the various encodings
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <iostream>
#include <algorithm>
#include <cstdlib>
#include "neighborhood.h"
#include "lsmove.h"
#include "chromosome.h"
#include "encoding.h"
#include "kvparse/kvparse.h"
#include "keywords.h"
#include "mtrandom.h"

using namespace std;

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
neighborhood<Chromosome,Encoding>::neighborhood()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
neighborhood<Chromosome,Encoding>::~neighborhood()
{
}

/*!
 * \brief initialize the neighborhood with a seed value
 */
template <template <typename> class Chromosome, typename Encoding>
void neighborhood<Chromosome,Encoding>::initialize(const Chromosome<Encoding>& sol)
{
	m_base = sol;
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
hamming_neighborhood_impl<Chromosome,Encoding>::hamming_neighborhood_impl()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
hamming_neighborhood_impl<Chromosome,Encoding>::~hamming_neighborhood_impl()
{
}

/*!
 * \brief compute the hamming distance between two chromosomes
 */
template <template <typename> class Chromosome, typename Encoding>
typename Encoding::Genotype hamming_neighborhood_impl<Chromosome,Encoding>::distance_between(const Chromosome<Encoding>& c1,
        const Chromosome<Encoding>& c2) const
{
	typename Encoding::Genotype d = 0;
	for(unsigned int i=0; i<c1.length(); i++) {
		if(c1[i] != c2[i]) {
			d++;
		}
	}
	return d;
}

/*!
 * \brief compute distance between chromosomes and number of infeasible solutions on path
 */
template <template <typename> class Chromosome, typename Encoding>
void hamming_neighborhood_impl<Chromosome,Encoding>::feasible_distance_between(const Chromosome<Encoding>& c1,
        const Chromosome<Encoding>& c2,
        typename Encoding::Genotype& dist,
        typename Encoding::Genotype& infeas,
        const typename Encoding::ProblemType* prob) const
{
	dist = 0;
	infeas = 0;
	Chromosome<Encoding> temp(c1);
	for(unsigned int i=0; i<c1.length(); i++) {
		if(temp[i] != c2[i]) {
			temp[i] = c2[i];
			temp.evaluate(prob);

			dist++;
			if(!temp.feasible()) {
				infeas++;
			}
		}
	}
}

/*!
 * \brief randomize the order in which to flip bits
 */
template <template <typename> class Chromosome, typename Encoding>
void hamming_neighborhood_impl<Chromosome,Encoding>::initialize(const Chromosome<Encoding>& sol)
{
	neighborhood<Chromosome,Encoding>::initialize(sol);

	mtrandom mt;
	order = mt.permutation(sol.length());
	index = 0;
}

/*!
 * \brief check if more neighbors exist
 */
template <template <typename> class Chromosome, typename Encoding>
bool hamming_neighborhood_impl<Chromosome,Encoding>::has_more_neighbors() const
{
	return index < this->m_base.length();
}

/*!
 * \brief return the next available move
 */
template <template <typename> class Chromosome, typename Encoding>
lsmove<Chromosome,Encoding> hamming_neighborhood_impl<Chromosome,Encoding>::next_neighbor()
{
	this->m_current= this->m_base;
	lsmove<Chromosome,Encoding> mv;
	mv.add_component(order[index],!this->m_current[order[index]]);
	mv.apply(this->m_current);
	index++;
	return mv;
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
swap_neighborhood_impl<Chromosome,Encoding>::swap_neighborhood_impl()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
swap_neighborhood_impl<Chromosome,Encoding>::~swap_neighborhood_impl()
{
}

/*!
 * \brief compute swap distance between chromosomes
 */
template <template <typename> class Chromosome, typename Encoding>
typename Encoding::Genotype swap_neighborhood_impl<Chromosome,Encoding>::distance_between(const Chromosome<Encoding>& c1,
        const Chromosome<Encoding>& c2) const
{
	typename Encoding::Genotype dist = 0;
	Chromosome<Encoding> temp(c1);
	for(int i=0; i<int(temp.length()); i++) {
		typename Chromosome<Encoding>::iterator it = find(temp.begin(), temp.end(), c2[i]);
		int pos = int(it - temp.begin());
		if(pos != i) {
			swap(temp[i],temp[pos]);
			dist++;
		}
	}
	return dist;
}

/*!
 * \brief compute swap distance and number of infeasible solutions along path
 */
template <template <typename> class Chromosome, typename Encoding>
void swap_neighborhood_impl<Chromosome,Encoding>::feasible_distance_between(const Chromosome<Encoding>& c1,
        const Chromosome<Encoding>& c2,
        typename Encoding::Genotype& dist,
        typename Encoding::Genotype& infeas,
        const typename Encoding::ProblemType* prob) const
{
	dist = 0;
	infeas = 0;
	Chromosome<Encoding> temp(c1);
	for(int i=0; i<int(temp.length()); i++) {
		typename Chromosome<Encoding>::iterator it = find(temp.begin(), temp.end(), c2[i]);
		int pos = int(it - temp.begin());
		if(pos != i) {
			swap(temp[i],temp[pos]);
			temp.evaluate(prob);

			dist++;
			if(!temp.feasible()) {
				infeas++;
			}
		}
	}
}

/*!
 * \brief randomize the order of allele visitation
 */
template <template <typename> class Chromosome, typename Encoding>
void swap_neighborhood_impl<Chromosome,Encoding>::initialize(const Chromosome<Encoding>& sol)
{
	neighborhood<Chromosome,Encoding>::initialize(sol);

	mtrandom mt;
	order = mt.permutation(sol.length());
	i = 0;
	j = 1;
}

/*!
 * \brief check if more neighbors exist
 */
template <template <typename> class Chromosome, typename Encoding>
bool swap_neighborhood_impl<Chromosome,Encoding>::has_more_neighbors() const
{
	return i != this->m_base.length() - 1;
}

/*!
 * \brief return the next available swap
 */
template <template <typename> class Chromosome, typename Encoding>
lsmove<Chromosome,Encoding> swap_neighborhood_impl<Chromosome,Encoding>::next_neighbor()
{
	this->m_current= this->m_base;
	lsmove<Chromosome,Encoding> m;
	m.add_component(order[i],this->m_current[order[j]]);
	m.add_component(order[j++],this->m_current[order[i]]);
	if(j == this->m_base.length()) {
		i++;
		j = i+1;
	}
	return m;
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
shift_neighborhood_impl<Chromosome,Encoding>::shift_neighborhood_impl()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
shift_neighborhood_impl<Chromosome,Encoding>::~shift_neighborhood_impl()
{
}

/*!
 * \brief compute shift distance between chromosomes
 */
template <template <typename> class Chromosome, typename Encoding>
typename Encoding::Genotype shift_neighborhood_impl<Chromosome,Encoding>::distance_between(const Chromosome<Encoding>& c1,
        const Chromosome<Encoding>& c2) const
{
	typename Encoding::Genotype dist = 0;
	for(unsigned int i=0; i<c1.length(); i++) {
		if(c1[i] != c2[i]) {
			dist++;
		}
	}
	return dist;
}

/*!
 * \brief compute shift distance with infeasibility of intermediates
 */
template <template <typename> class Chromosome, typename Encoding>
void shift_neighborhood_impl<Chromosome,Encoding>::feasible_distance_between(const Chromosome<Encoding>& c1,
        const Chromosome<Encoding>& c2,
        typename Encoding::Genotype& dist,
        typename Encoding::Genotype& infeas,
        const typename Encoding::ProblemType* prob) const
{
	dist = 0;
	infeas = 0;
	Chromosome<Encoding> temp(c1);
	for(unsigned int i=0; i<c1.length(); i++) {
		if(temp[i] != c2[i]) {
			temp[i] = c2[i];
			temp.evaluate(prob);

			dist++;
			if(!temp.feasible()) {
				infeas++;
			}
		}
	}
}

/*!
 * \brief randomize agent/task visitation orders
 */
template <template <typename> class Chromosome, typename Encoding>
void shift_neighborhood_impl<Chromosome,Encoding>::initialize(const Chromosome<Encoding>& sol)
{
	neighborhood<Chromosome,Encoding>::initialize(sol);

	mtrandom mt;
	task_index = 0;
	agent_index = 0;
	task_order = mt.permutation(sol.tasks());
	agent_order = mt.permutation(sol.agents());
}

/*!
 * \brief check if more neighbors exist
 */
template <template <typename> class Chromosome, typename Encoding>
bool shift_neighborhood_impl<Chromosome,Encoding>::has_more_neighbors() const
{
	return task_index != this->m_base.tasks();
}

/*!
 * \brief return next available move
 */
template <template <typename> class Chromosome, typename Encoding>
lsmove<Chromosome,Encoding> shift_neighborhood_impl<Chromosome,Encoding>::next_neighbor()
{
	this->m_current= this->m_base;
	lsmove<Chromosome,Encoding> m;
	m.add_component(task_order[task_index],agent_order[agent_index++]);
	if(agent_index == this->m_base.agents()) {
		mtrandom mt;
		agent_index = 0;
		agent_order = mt.permutation(this->m_base.agents());
		task_index++;
	}
	return m;
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
sss_neighborhood_impl<Chromosome,Encoding>::sss_neighborhood_impl()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
sss_neighborhood_impl<Chromosome,Encoding>::~sss_neighborhood_impl()
{
}

/*!
 * \brief compute shift/swap distance between chromosomes
 */
template <template <typename> class Chromosome, typename Encoding>
typename Encoding::Genotype sss_neighborhood_impl<Chromosome,Encoding>::distance_between(const Chromosome<Encoding>& c1,
        const Chromosome<Encoding>& c2) const
{
	typename Encoding::Genotype dist = 0;
	Chromosome<Encoding> temp(c1);
	for(unsigned int i=0; i<temp.length(); i++) {
		// check if the two chromosomes differ in the given location
		if(temp[i] != c2[i]) {
			// if so, then we need to figure out whether an appropriate swap exists as follows:
			// we must change c1[i] to c2[i] to move c1 one step closer to c2.  However,
			// if another element of c1 is equal to c2[i], then we can swap c1[i] and c1[other]
			// and kill two birds with one move
			typename Chromosome<Encoding>::iterator it = find(temp.begin()+i, temp.end(), c2[i]);
			int pos = 0;
			while(it != temp.end()) {
				pos = static_cast<int>(it - temp.begin());
				if(c2[pos] == temp[i]) {
					// if there exists an element of c1 equal to c2[pos], then we can
					// do a swap of c1[i] and c1[pos]
					break;
				} else {
					// otherwise, find the next element of c1 equal to c2[i]
					it = find(it+1, temp.end(), c2[i]);
				}
			}

			if(it != temp.end()) {
				// if we broke out of the loop before it == temp.end(), then a good
				// swap was found
				swap(temp[i], temp[pos]);
			} else {
				// otherwise, no swap does any good, so we just perform a shift
				temp[i] = c2[i];
			}

			// increment the number of moves regardless of what type of move was made
			dist++;
		}
	}
	return dist;
}

/*!
 * \brief compute shift/swap distance with infeasibility of intermediates
 */
template <template <typename> class Chromosome, typename Encoding>
void sss_neighborhood_impl<Chromosome,Encoding>::feasible_distance_between(const Chromosome<Encoding>& c1,
        const Chromosome<Encoding>& c2,
        typename Encoding::Genotype& dist,
        typename Encoding::Genotype& infeas,
        const typename Encoding::ProblemType* prob) const
{
	dist = 0;
	infeas = 0;
	Chromosome<Encoding> temp(c1);
	for(unsigned int i=0; i<temp.length(); i++) {
		// check if the two chromosomes differ in the given location
		if(temp[i] != c2[i]) {
			// if so, then we need to figure out whether an appropriate swap exists as follows:
			// we must change c1[i] to c2[i] to move c1 one step closer to c2.  However,
			// if another element of c1 is equal to c2[i], then we can swap c1[i] and c1[other]
			// and kill two birds with one move
			typename Chromosome<Encoding>::iterator it = find(temp.begin()+i, temp.end(), c2[i]);
			int pos = 0;
			while(it != temp.end()) {
				pos = static_cast<int>(it - temp.begin());
				if(c2[pos] == temp[i]) {
					// if there exists an element of c1 equal to c2[pos], then we can
					// do a swap of c1[i] and c1[pos]
					break;
				} else {
					// otherwise, find the next element of c1 equal to c2[i]
					it = find(it+1, temp.end(), c2[i]);
				}
			}

			if(it != temp.end()) {
				// if we broke out of the loop before it == temp.end(), then a good
				// swap was found
				swap(temp[i], temp[pos]);
			} else {
				// otherwise, no swap does any good, so we just perform a shift
				temp[i] = c2[i];
			}

			// check if the move is feasible
			temp.evaluate(prob);

			dist++;
			if(!temp.feasible()) {
				infeas++;
			}
		}
	}
}

/*!
 * \brief initialize each constituent neighborhood
 */
template <template <typename> class Chromosome, typename Encoding>
void sss_neighborhood_impl<Chromosome,Encoding>::initialize(const Chromosome<Encoding>& sol)
{
	shift_part.initialize(sol);
	swap_part.initialize(sol);
}

/*!
 * \brief return false only if both neighborhoods have been exhausted
 */
template <template <typename> class Chromosome, typename Encoding>
bool sss_neighborhood_impl<Chromosome,Encoding>::has_more_neighbors() const
{
	return (swap_part.has_more_neighbors() || shift_part.has_more_neighbors());
}

/*!
 * \brief return next shift or swap lsmove at random
 */
template <template <typename> class Chromosome, typename Encoding>
lsmove<Chromosome,Encoding> sss_neighborhood_impl<Chromosome,Encoding>::next_neighbor()
{
	mtrandom mt;

	if(mt.random() < 0.5) {
		if(shift_part.has_more_neighbors()) {
			return shift_part.next_neighbor();
		} else {
			return swap_part.next_neighbor();
		}
	} else {
		if(swap_part.has_more_neighbors()) {
			return swap_part.next_neighbor();
		} else {
			return shift_part.next_neighbor();
		}
	}
}

/*!
 * \brief create a binary neighborhood operator
 */
template <template <typename> class Chromosome, typename Encoding>
neighborhood<Chromosome,Encoding>* bit_neighborhood_factory<Chromosome,Encoding>::construct()
{
	string op;
	kvparse::parameter_value(keywords::NEIGHBORHOOD, op, true);
	if(op == keywords::HAMMING_NEIGHBORHOOD) {
		hamming_neighborhood<Chromosome,Encoding>* n = new hamming_neighborhood<Chromosome,Encoding>;
		return n;
	} else if(op == keywords::SWAP_NEIGHBORHOOD) {
		swap_neighborhood<Chromosome,Encoding>* n = new swap_neighborhood<Chromosome,Encoding>;
		return n;
	} else {
		cerr << "illegal neighborhood operator: " << op << " specified" << endl;
		exit(1);
		return 0;
	}
}

/*!
 * \brief create a permutation neighborhood operator
 */
template <template <typename> class Chromosome, typename Encoding>
neighborhood<Chromosome,Encoding>* permutation_neighborhood_factory<Chromosome,Encoding>::construct()
{
	string op;
	kvparse::parameter_value(keywords::NEIGHBORHOOD, op, true);
	if(op == keywords::SWAP_NEIGHBORHOOD) {
		swap_neighborhood<Chromosome,Encoding>* n = new swap_neighborhood<Chromosome,Encoding>;
		return n;
	} else {
		cerr << "illegal neighborhood operator specified: " << op << endl;
		exit(1);
		return 0;
	}
}

/*!
 * \brief create a gap neighborhood operator
 */
template <template <typename> class Chromosome, typename Encoding>
neighborhood<Chromosome,Encoding>* gap_neighborhood_factory<Chromosome,Encoding>::construct()
{
	string op;
	kvparse::parameter_value(keywords::NEIGHBORHOOD, op, true);
	if(op == keywords::SWAP_NEIGHBORHOOD) {
		swap_neighborhood<Chromosome,Encoding>* n = new swap_neighborhood<Chromosome,Encoding>;
		return n;
	} else if(op == keywords::SHIFT_NEIGHBORHOOD) {
		shift_neighborhood<Chromosome,Encoding>* n = new shift_neighborhood<Chromosome,Encoding>;
		return n;
	} else if(op == keywords::SSS_NEIGHBORHOOD) {
		sss_neighborhood<Chromosome,Encoding>* n = new sss_neighborhood<Chromosome,Encoding>;
		return n;
	} else {
		cerr << "illegal neighborhood operator specified: " << op << endl;
		exit(1);
		return 0;
	}
}

/*!
 * \brief create a gsap neighborhood operator
 */
template <template <typename> class Chromosome, typename Encoding>
neighborhood<Chromosome,Encoding>* gsap_neighborhood_factory<Chromosome,Encoding>::construct()
{
	string op;
	kvparse::parameter_value(keywords::NEIGHBORHOOD, op, true);
	if(op == keywords::SWAP_NEIGHBORHOOD) {
		swap_neighborhood<Chromosome,Encoding>* n = new swap_neighborhood<Chromosome,Encoding>;
		return n;
	} else if(op == keywords::SHIFT_NEIGHBORHOOD) {
		shift_neighborhood<Chromosome,Encoding>* n = new shift_neighborhood<Chromosome,Encoding>;
		return n;
	} else if(op == keywords::SSS_NEIGHBORHOOD) {
		sss_neighborhood<Chromosome,Encoding>* n = new sss_neighborhood<Chromosome,Encoding>;
		return n;
	} else {
		cerr << "illegal neighborhood operator specified: " << op << endl;
		exit(1);
		return 0;
	}
}

/*!
 * \brief other neighborhoods not implemented
 */
template <template <typename> class Chromosome, typename Encoding>
neighborhood<Chromosome,Encoding>* neighborhood_factory<Chromosome,Encoding>::construct()
{
	cerr << "neighborhood not implemented" << endl;
	exit(1);
	return 0;
}

