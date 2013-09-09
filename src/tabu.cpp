/*!
 * \file tabu.cpp
 *
 * tabu search algorithms
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cfloat>
#include "sls.h"
#include "localsearch.h"
#include "move.h"
#include "neighborhood.h"
#include "chromosome.h"
#include "encoding.h"
#include "comparator.h"
#include "problems.h"
#include "configuration.h"
#include "keywords.h"

using namespace std;

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
tabu_list<Chromosome,Encoding>::tabu_list()
{
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
tabu_list<Chromosome,Encoding>::~tabu_list()
{
}

/**
 * \brief set the configuration prefix
 */
template <template <typename> class Chromosome, typename Encoding>
void tabu_list<Chromosome,Encoding>::set_prefix(const string& prefix)
{
	this->m_prefix=prefix;
}

/*!
 * \brief initialize the tabu list
 */
template <template <typename> class Chromosome, typename Encoding>
void tabu_list<Chromosome,Encoding>::initialize()
{
	_internal.clear();
	configuration::unsigned_integer_parameter(this->m_prefix+keywords::MIN_TABU_TENURE,_ttmin,true);
	configuration::unsigned_integer_parameter(this->m_prefix+keywords::MAX_TABU_TENURE,_ttmax,true);
}

/*!
 * \brief erase all items from the tabu list
 */
template <template <typename> class Chromosome, typename Encoding>
void tabu_list<Chromosome,Encoding>::clear()
{
	_internal.clear();
}

/*!
 * \brief update the tabu list with an accepted move
 */
template <template <typename> class Chromosome, typename Encoding>
void tabu_list<Chromosome,Encoding>::accept_move(const Chromosome<Encoding>& chr,
        const move<Chromosome,Encoding>& m,
        unsigned int iter)
{
	mtrandom mt;

	//! create a move that would undo the accepted move (prevent the chromosome
	//! from receiving the values it currently has at the affected positions)
	move<Chromosome,Encoding> revert;
	for(typename move<Chromosome,Encoding>::const_iterator mi=m.begin(); mi!=m.end(); mi++) {
		revert.add_component(mi->first,chr[mi->first]);
	}

	tlist_item i=make_pair(revert,iter+mt.random(static_cast<int>(_ttmin),
	                       static_cast<int>(_ttmax)));

	// put the new item on the back of the queue
	_internal.push_back(i);

	// if the queue is too large, remove an item from the front
	while(iter>_internal[0].second) {
		_internal.pop_front();
	}

}

/*!
 * \brief determine if a move is tabu at the given generation
 */
template <template <typename> class Chromosome, typename Encoding>
bool tabu_list<Chromosome,Encoding>::tabu(const move<Chromosome,Encoding>& m, unsigned int iter) const
{
	for(typename deque<tlist_item>::const_iterator i=_internal.begin(); i!=_internal.end(); i++) {
		//! check to see if any component of the move is part of a tabu move
		move<Chromosome,Encoding> cm=i->first;
		unsigned int tenure=i->second;

		//! for each component mi in m
		//!     if mi appears in curr.first and iter<=curr.second
		//!         return true
		//! return true
		for(typename move<Chromosome,Encoding>::const_iterator mi=m.begin(); mi!=m.end(); mi++) {
			if((find(cm.begin(),cm.end(),(*mi))!=cm.end()) &&
			        (iter<=tenure)) {
				return true;
			}
		}
	}
	return false;
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
tabu_search<Chromosome,Encoding>::tabu_search()
{
	_tlist=new tabu_list<Chromosome,Encoding>;
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
tabu_search<Chromosome,Encoding>::~tabu_search()
{
	delete _tlist;
}

/**
 * \brief set the configuration prefix
 */
template <template <typename> class Chromosome, typename Encoding>
void tabu_search<Chromosome,Encoding>::set_prefix(const string& prefix)
{
	local_search<Chromosome,Encoding>::set_prefix(prefix);
	_tlist->set_prefix(prefix);
}

/*!
 * \brief initialize the tabu search components
 */
template <template <typename> class Chromosome, typename Encoding>
void tabu_search<Chromosome,Encoding>::initialize()
{
	local_search<Chromosome,Encoding>::initialize();
	_tlist->initialize();
}

/*!
 * \brief reset the algorithm parameters
 */
template <template <typename> class Chromosome, typename Encoding>
void tabu_search<Chromosome,Encoding>::reset()
{
	local_search<Chromosome,Encoding>::reset();
	_tlist->clear();
	_tlist->initialize();
}

/*!
 * \brief optimize the given individual via tabu search
 */
template <template <typename> class Chromosome, typename Encoding>
void tabu_search<Chromosome,Encoding>::improve(Chromosome<Encoding>& chr,
        comparator<Chromosome,Encoding>* comp,
        const typename Encoding::ProblemType* prob)
{
	this->reset();

	unsigned int iter=0;

	//! keep track of current best known solution
	Chromosome<Encoding> best=chr;

	while(!this->terminate()) {
		//! keep track of best neighbor
		Chromosome<Encoding> best_neighbor=chr;
		move<Chromosome,Encoding> best_move;
		bool already_aspired=false;
		bool first_one=true;

		//! initialize the neighborhood
		this->m_nf->initialize(chr);

		//! for each move in neighborhood
		while(this->m_nf->has_more_neighbors()) {
			//! determine if move is tabu, aspired
			move<Chromosome,Encoding> m=this->m_nf->next_neighbor();
			bool istabu=_tlist->tabu(m,iter);
			bool aspired=false;
			Chromosome<Encoding> tmp=chr;
			m.apply(tmp);
			if(this->m_repair) {
				this->m_repair->repair(tmp,prob);
			}
			tmp.evaluate(prob);
			this->chromosome_evaluated(tmp);

			//! if better than the previous best, the move
			//! is aspired
			if(comp->compare(tmp,best)<0) {
				aspired=true;
			}

			//! if first aspired neighbor, or if better than
			//! any prior aspired neighbors, or if better than
			//! prior best neighbor and not tabu, make this the
			//! new best
			if((aspired && !already_aspired) ||
			        (aspired && already_aspired && (first_one || comp->compare(tmp,best_neighbor)<0)) ||
			        (!aspired && !already_aspired && (first_one || comp->compare(tmp,best_neighbor)<0) && !istabu)) {
				first_one=false;
				best_neighbor=tmp;
				best_move=m;
				if(aspired) {
					already_aspired=true;
				}
			}
		}

		//! if no best neighbor found, there was a problem
		if(first_one) {
			error("all moves are tabu!");
		}

		//! otherwise, accept the move and update the tabu list
		_tlist->accept_move(chr,best_move,iter);
		chr=best_neighbor;

		//! update best known solution
		if(comp->compare(chr,best)<0) {
			best=chr;
		}

		if(this->debug_generations) {
			cout << chr << endl;
		}

		this->generation_completed();
		iter++;
	}

	//! clear the tabu list
	_tlist->clear();
}
