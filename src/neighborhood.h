/*!
 * \file neighborhood.h
 *
 * Defines neighborhood operators for the various encodings
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _NEIGHBORHOOD_H_
#define _NEIGHBORHOOD_H_

#include <vector>
#include "encoding.h"
#include "move.h"

/*!
 * \class neighborhood
 *
 * provides iterator-like functionality over a local search neighborhood
 */
template <template <typename> class Chromosome, typename Encoding>
class neighborhood
{
protected:
	Chromosome<Encoding> m_base;
	Chromosome<Encoding> m_current;

private:
	// disable copying
	neighborhood(const neighborhood& that);
	neighborhood& operator=(const neighborhood& that);

public:
	neighborhood();
	virtual ~neighborhood();

	virtual typename Encoding::Genotype distance_between(const Chromosome<Encoding>& c1,
	        const Chromosome<Encoding>& c2) const = 0;
	virtual void feasible_distance_between(const Chromosome<Encoding>& c1,
	                                       const Chromosome<Encoding>& c2,
	                                       typename Encoding::Genotype& dist,
	                                       typename Encoding::Genotype& feas,
	                                       const typename Encoding::ProblemType* prob) const = 0;
	virtual void initialize(const Chromosome<Encoding>& sol);
	virtual bool has_more_neighbors() const = 0;
	virtual move<Chromosome,Encoding> next_neighbor() = 0;
};

/*!
 * \class hamming_neighborhood_impl
 *
 * all solutions reachable in a single bit-flip from a seed value
 */
template <template <typename> class Chromosome, typename Encoding>
class hamming_neighborhood_impl : public neighborhood<Chromosome,Encoding>
{
protected:
	vector<int> order;
	unsigned int index;

public:
	hamming_neighborhood_impl();
	virtual ~hamming_neighborhood_impl();

	virtual typename Encoding::Genotype distance_between(const Chromosome<Encoding>& c1,
	        const Chromosome<Encoding>& c2) const;
	virtual void feasible_distance_between(const Chromosome<Encoding>& c1,
	                                       const Chromosome<Encoding>& c2,
	                                       typename Encoding::Genotype& dist,
	                                       typename Encoding::Genotype& feas,
	                                       const typename Encoding::ProblemType* prob) const;
	virtual void initialize(const Chromosome<Encoding>& sol);
	virtual bool has_more_neighbors() const;
	virtual move<Chromosome,Encoding> next_neighbor();
};

/*!
 * \class hamming_neighborhood
 *
 * interface to hamming_neighborhood_impl
 */
template <template <typename> class Chromosome, typename Encoding> class hamming_neighborhood;
template <template <typename> class Chromosome>
class hamming_neighborhood<Chromosome,binary_encoding> :
	public hamming_neighborhood_impl<Chromosome,binary_encoding>
{
};

/*!
 * \class hamming_neighborhood
 *
 * interface to hamming_neighborhood_impl
 */
template <template <typename> class Chromosome>
class hamming_neighborhood<Chromosome,boolean_encoding> :
	public hamming_neighborhood_impl<Chromosome,boolean_encoding>
{
};

/*!
 * \class swap_neighborhood_impl
 *
 * all solutions reachable in a single swap from a seed value
 */
template <template <typename> class Chromosome, typename Encoding>
class swap_neighborhood_impl : public neighborhood<Chromosome,Encoding>
{
protected:
	vector<int> order;
	unsigned int i;
	unsigned int j;

private:
	// disable copying
	swap_neighborhood_impl(const swap_neighborhood_impl& that);
	swap_neighborhood_impl& operator=(const swap_neighborhood_impl& that);

public:
	swap_neighborhood_impl();
	virtual ~swap_neighborhood_impl();

	virtual typename Encoding::Genotype distance_between(const Chromosome<Encoding>& c1,
	        const Chromosome<Encoding>& c2) const;
	virtual void feasible_distance_between(const Chromosome<Encoding>& c1,
	                                       const Chromosome<Encoding>& c2,
	                                       typename Encoding::Genotype& dist,
	                                       typename Encoding::Genotype& feas,
	                                       const typename Encoding::ProblemType* prob) const;
	virtual void initialize(const Chromosome<Encoding>& sol);
	virtual bool has_more_neighbors() const;
	virtual move<Chromosome,Encoding> next_neighbor();
};

/*!
 * \class swap_neighborhood
 *
 * interface to swap_neighborhood_impl
 */
template <template <typename> class Chromosome, typename Encoding> class swap_neighborhood;
template <template <typename> class Chromosome>
class swap_neighborhood<Chromosome,permutation_encoding> :
	public swap_neighborhood_impl<Chromosome,permutation_encoding>
{
};

/*!
 * \class hamming_neighborhood
 *
 * interface to hamming_neighborhood_impl
 */
template <template <typename> class Chromosome>
class swap_neighborhood<Chromosome,boolean_encoding> :
	public swap_neighborhood_impl<Chromosome,boolean_encoding>
{
};

/*!
 * \class hamming_neighborhood
 *
 * interface to hamming_neighborhood_impl
 */
template <template <typename> class Chromosome>
class swap_neighborhood<Chromosome,binary_encoding> :
	public swap_neighborhood_impl<Chromosome,binary_encoding>
{
};

/*!
 * \class swap_neighborhood
 *
 * interface to swap_neighborhood_impl
 */
template <template <typename> class Chromosome>
class swap_neighborhood<Chromosome,gap_encoding> :
	public swap_neighborhood_impl<Chromosome,gap_encoding>
{
};

/*!
 * \class swap_neighborhood
 *
 * interface to swap_neighborhood_impl
 */
template <template <typename> class Chromosome>
class swap_neighborhood<Chromosome,gsap_encoding> :
	public swap_neighborhood_impl<Chromosome,gsap_encoding>
{
};

/*!
 * \class shift_neighborhood_impl
 *
 * all solutions reachable with a single gap shift operation
 */
template <template <typename> class Chromosome, typename Encoding>
class shift_neighborhood_impl : public neighborhood<Chromosome,Encoding>
{
protected:
	vector<int> task_order;
	vector<int> agent_order;
	unsigned int task_index;
	unsigned int agent_index;

public:
	shift_neighborhood_impl();
	virtual ~shift_neighborhood_impl();

	virtual typename Encoding::Genotype distance_between(const Chromosome<Encoding>& c1,
	        const Chromosome<Encoding>& c2) const;
	virtual void feasible_distance_between(const Chromosome<Encoding>& c1,
	                                       const Chromosome<Encoding>& c2,
	                                       typename Encoding::Genotype& dist,
	                                       typename Encoding::Genotype& feas,
	                                       const typename Encoding::ProblemType* prob) const;
	virtual void initialize(const Chromosome<Encoding>& sol);
	virtual bool has_more_neighbors() const;
	virtual move<Chromosome,Encoding> next_neighbor();
};

/*!
 * \class shift_neighborhood
 *
 * interface to shift_neighborhood_impl
 */
template <template <typename> class Chromosome, typename Encoding> class shift_neighborhood;
template <template <typename> class Chromosome>
class shift_neighborhood<Chromosome,gap_encoding> :
	public shift_neighborhood_impl<Chromosome,gap_encoding>
{
};
template <template <typename> class Chromosome>
class shift_neighborhood<Chromosome,gsap_encoding> :
	public shift_neighborhood_impl<Chromosome,gsap_encoding>
{
};

/*!
 * \class sss_neighborhood_impl
 *
 * implements "shift and subsequent swap" neighborhood
 */
template <template <typename> class Chromosome, typename Encoding>
class sss_neighborhood_impl : public neighborhood<Chromosome,Encoding>
{
protected:
	shift_neighborhood<Chromosome,Encoding> shift_part;
	swap_neighborhood<Chromosome,Encoding> swap_part;

public:
	sss_neighborhood_impl();
	virtual ~sss_neighborhood_impl();

	virtual typename Encoding::Genotype distance_between(const Chromosome<Encoding>& c1,
	        const Chromosome<Encoding>& c2) const;
	virtual void feasible_distance_between(const Chromosome<Encoding>& c1,
	                                       const Chromosome<Encoding>& c2,
	                                       typename Encoding::Genotype& dist,
	                                       typename Encoding::Genotype& feas,
	                                       const typename Encoding::ProblemType* prob) const;
	virtual void initialize(const Chromosome<Encoding>& sol);
	virtual bool has_more_neighbors() const;
	virtual move<Chromosome,Encoding> next_neighbor();
};

/*!
 * \class sss_neighborhood
 *
 * interface to sss_neighborhood_impl
 */
template <template <typename> class Chromosome, typename Encoding> class sss_neighborhood;
template <template <typename> class Chromosome>
class sss_neighborhood<Chromosome,gap_encoding> :
	public sss_neighborhood_impl<Chromosome,gap_encoding>
{
};
template <template <typename> class Chromosome>
class sss_neighborhood<Chromosome,gsap_encoding> :
	public sss_neighborhood_impl<Chromosome,gsap_encoding>
{
};

/*!
 * \class permutation_neighborhood_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class permutation_neighborhood_factory
{
public:
	static neighborhood<Chromosome,Encoding>* construct();
};

/*!
 * \class bit_neighborhood_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class bit_neighborhood_factory
{
public:
	static neighborhood<Chromosome,Encoding>* construct();
};

/*!
 * \class gap_neighborhood_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class gap_neighborhood_factory
{
public:
	static neighborhood<Chromosome,Encoding>* construct();
};

/*!
 * \class gsap_neighborhood_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class gsap_neighborhood_factory
{
public:
	static neighborhood<Chromosome,Encoding>* construct();
};

/*!
 * \class neighborhood_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class neighborhood_factory
{
public:
	static neighborhood<Chromosome,Encoding>* construct();
};

/*!
 * \class neighborhood_factory
 */
template <template <typename> class Chromosome>
class neighborhood_factory<Chromosome,binary_encoding> :
	public bit_neighborhood_factory<Chromosome,binary_encoding>
{
};

/*!
 * \class neighborhood_factory
 */
template <template <typename> class Chromosome>
class neighborhood_factory<Chromosome,boolean_encoding> :
	public bit_neighborhood_factory<Chromosome,boolean_encoding>
{
};

/*!
 * \class neighborhood_factory
 */
template <template <typename> class Chromosome>
class neighborhood_factory<Chromosome,permutation_encoding>
	: public permutation_neighborhood_factory<Chromosome,permutation_encoding>
{
};

/*!
 * \class neighborhood_factory
 */
template <template <typename> class Chromosome>
class neighborhood_factory<Chromosome,gap_encoding>
	: public gap_neighborhood_factory<Chromosome,gap_encoding>
{
};

/*!
 * \class neighborhood_factory
 */
template <template <typename> class Chromosome>
class neighborhood_factory<Chromosome,gsap_encoding>
	: public gap_neighborhood_factory<Chromosome,gsap_encoding>
{
};

#include "neighborhood.cpp"

#endif
