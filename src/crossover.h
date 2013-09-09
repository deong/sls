/*!
 * \file crossover.h
 *
 * evolutionary algorithm crossover operators
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _CROSSOVER_H_
#define _CROSSOVER_H_

#include "chromosome.h"
#include "encoding.h"

using namespace std;

/*!
 * \class crossover_operator
 *
 * base class for all crossover operators
 */
template <template <typename> class Chromosome, typename Encoding>
class crossover_operator
{
private:
    // disallow copying of functional classes
    crossover_operator(const crossover_operator& that);
    crossover_operator& operator=(const crossover_operator& that);

public:
    crossover_operator();
    virtual ~crossover_operator();

    virtual void crossover(const Chromosome<Encoding>& p1,
                           const Chromosome<Encoding>& p2,
                           Chromosome<Encoding>& c1,
                           Chromosome<Encoding>& c2) const = 0;
};

/*!
 * \class uniform_crossover_impl
 *
 * implementation class for uniform crossover
 */
template <template <typename> class Chromosome, typename Encoding>
class uniform_crossover_impl : public crossover_operator<Chromosome,Encoding> 
{
private:
    // disallow copying of functional classes
    uniform_crossover_impl(const uniform_crossover_impl& that);
    uniform_crossover_impl& operator=(const uniform_crossover_impl& that);

public:
    uniform_crossover_impl();
    virtual ~uniform_crossover_impl();

    virtual void crossover(const Chromosome<Encoding>& p1,
                           const Chromosome<Encoding>& p2,
                           Chromosome<Encoding>& c1,
                           Chromosome<Encoding>& c2) const;
};

/*!
 * \class uniform_crossover
 *
 * interface for uniform_crossover_impl
 */
template <template <typename> class Chromosome, typename Encoding> class uniform_crossover;
template <template <typename> class Chromosome>
class uniform_crossover<Chromosome,boolean_encoding> :
    public uniform_crossover_impl<Chromosome,boolean_encoding>
{
};

/*!
 * \class uniform_crossover
 *
 * interface for uniform_crossover_impl
 */
template <template <typename> class Chromosome>
class uniform_crossover<Chromosome,binary_encoding> :
    public uniform_crossover_impl<Chromosome,binary_encoding>
{
};

/*!
 * \class uniform_crossover
 *
 * interface for uniform_crossover_impl
 */
template <template <typename> class Chromosome>
class uniform_crossover<Chromosome,real_encoding> :
    public uniform_crossover_impl<Chromosome,real_encoding>
{
};

/*!
 * \class uniform_crossover
 *
 * interface for uniform_crossover_impl
 */
template <template <typename> class Chromosome>
class uniform_crossover<Chromosome,integer_encoding> :
    public uniform_crossover_impl<Chromosome,integer_encoding>
{
};

/*!
 * \class uniform_crossover
 *
 * interface for uniform_crossover_impl
 */
template <template <typename> class Chromosome>
class uniform_crossover<Chromosome,gap_encoding> :
    public uniform_crossover_impl<Chromosome,gap_encoding>
{
};

/*!
 * \class uniform_crossover
 *
 * interface for uniform_crossover_impl
 */
template <template <typename> class Chromosome>
class uniform_crossover<Chromosome,gsap_encoding> :
    public uniform_crossover_impl<Chromosome,gsap_encoding>
{
};

/*!
 * \class one_point_crossover_impl
 *
 * implementation class for one-point crossover
 */
template <template <typename> class Chromosome, typename Encoding>
class one_point_crossover_impl : public crossover_operator<Chromosome,Encoding>
{
private:
    // disable copying of functional classes
    one_point_crossover_impl(const one_point_crossover_impl& that);
    one_point_crossover_impl& operator=(const one_point_crossover_impl& that);

public:
    one_point_crossover_impl();
    virtual ~one_point_crossover_impl();

    virtual void crossover(const Chromosome<Encoding>& p1,
                           const Chromosome<Encoding>& p2,
                           Chromosome<Encoding>& c1,
                           Chromosome<Encoding>& c2) const;
};

/*!
 * \class one_point_crossover
 *
 * interface to one_point_crossover_impl
 */
template <template <typename> class Chromosome, typename Encoding> class one_point_crossover;
template <template <typename> class Chromosome>
class one_point_crossover<Chromosome,boolean_encoding> :
    public one_point_crossover_impl<Chromosome,boolean_encoding>
{
};

/*!
 * \class one_point_crossover
 *
 * interface to one_point_crossover_impl
 */
template <template <typename> class Chromosome>
class one_point_crossover<Chromosome,binary_encoding> :
    public one_point_crossover_impl<Chromosome,binary_encoding>
{
};

/*!
 * \class one_point_crossover
 *
 * interface to one_point_crossover_impl
 */
template <template <typename> class Chromosome>
class one_point_crossover<Chromosome,real_encoding> :
    public one_point_crossover_impl<Chromosome,real_encoding>
{
};

/*!
 * \class one_point_crossover
 *
 * interface to one_point_crossover_impl
 */
template <template <typename> class Chromosome>
class one_point_crossover<Chromosome,integer_encoding> :
    public one_point_crossover_impl<Chromosome,integer_encoding>
{
};

/*!
 * \class one_point_crossover
 *
 * interface to one_point_crossover_impl
 */
template <template <typename> class Chromosome>
class one_point_crossover<Chromosome,gap_encoding> :
    public one_point_crossover_impl<Chromosome,gap_encoding>
{
};

/*!
 * \class one_point_crossover
 *
 * interface to one_point_crossover_impl
 */
template <template <typename> class Chromosome>
class one_point_crossover<Chromosome,gsap_encoding> :
    public one_point_crossover_impl<Chromosome,gsap_encoding>
{
};

/*!
 * \class two_point_crossover_impl
 *
 * implementation class for two-point crossover
 */
template <template <typename> class Chromosome, typename Encoding>
class two_point_crossover_impl : public crossover_operator<Chromosome,Encoding>
{
private:
    // disable copying of functional class
    two_point_crossover_impl(const two_point_crossover_impl& that);
    two_point_crossover_impl& operator=(const two_point_crossover_impl& that);

public:
    two_point_crossover_impl();
    virtual ~two_point_crossover_impl();

    virtual void crossover(const Chromosome<Encoding>& p1,
                           const Chromosome<Encoding>& p2,
                           Chromosome<Encoding>& c1,
                           Chromosome<Encoding>& c2) const;
};

/*!
 * \class two_point_crossover
 *
 * interface to two_point_crossover_impl
 */
template <template <typename> class Chromosome, typename Encoding> class two_point_crossover;
template <template <typename> class Chromosome>
class two_point_crossover<Chromosome,boolean_encoding> :
    public two_point_crossover_impl<Chromosome,boolean_encoding>
{
};

/*!
 * \class two_point_crossover
 *
 * interface to two_point_crossover_impl
 */
template <template <typename> class Chromosome>
class two_point_crossover<Chromosome,binary_encoding> :
    public two_point_crossover_impl<Chromosome,binary_encoding>
{
};

/*!
 * \class two_point_crossover
 *
 * interface to two_point_crossover_impl
 */
template <template <typename> class Chromosome>
class two_point_crossover<Chromosome,real_encoding> :
    public two_point_crossover_impl<Chromosome,real_encoding>
{
};

/*!
 * \class two_point_crossover
 *
 * interface to two_point_crossover_impl
 */
template <template <typename> class Chromosome>
class two_point_crossover<Chromosome,integer_encoding> :
    public two_point_crossover_impl<Chromosome,integer_encoding>
{
};

/*!
 * \class two_point_crossover
 *
 * interface to two_point_crossover_impl
 */
template <template <typename> class Chromosome>
class two_point_crossover<Chromosome,gap_encoding> :
    public two_point_crossover_impl<Chromosome,gap_encoding>
{
};

/*!
 * \class two_point_crossover
 *
 * interface to two_point_crossover_impl
 */
template <template <typename> class Chromosome>
class two_point_crossover<Chromosome,gsap_encoding> :
    public two_point_crossover_impl<Chromosome,gsap_encoding>
{
};

/*!
 * \class hux_crossover_impl
 *
 * implementation class for hux crossover (half-uniform crossover)
 */
template <template <typename> class Chromosome, typename Encoding>
class hux_crossover_impl : public crossover_operator<Chromosome,Encoding>
{
private:
    // disable copying of functional class
    hux_crossover_impl(const hux_crossover_impl& that);
    hux_crossover_impl& operator=(const hux_crossover_impl& that);

public:
    hux_crossover_impl();
    virtual ~hux_crossover_impl();

    virtual void crossover(const Chromosome<Encoding>& p1,
                           const Chromosome<Encoding>& p2,
                           Chromosome<Encoding>& c1,
                           Chromosome<Encoding>& c2) const;
};

/*!
 * \class hux_crossover
 *
 * interface to hux_crossover_impl
 */
template <template <typename> class Chromosome, typename Encoding> class hux_crossover;
template <template <typename> class Chromosome>
class hux_crossover<Chromosome,boolean_encoding> :
    public hux_crossover_impl<Chromosome,boolean_encoding>
{
};

/*!
 * \class hux_crossover
 *
 * interface to hux_crossover_impl
 */
template <template <typename> class Chromosome>
class hux_crossover<Chromosome,binary_encoding> :
    public hux_crossover_impl<Chromosome,binary_encoding>
{
};

/*!
 * \class sbx_crossover_impl
 *
 * implementation class for simulated binary crossover
 */
template <template <typename> class Chromosome, typename Encoding>
class sbx_crossover_impl : public crossover_operator<Chromosome,Encoding>
{
private:
    double m_eta;
    
public:
    sbx_crossover_impl();
    virtual ~sbx_crossover_impl();

    virtual void initialize();
    virtual void crossover(const Chromosome<Encoding>& p1,
                           const Chromosome<Encoding>& p2,
                           Chromosome<Encoding>& c1,
                           Chromosome<Encoding>& c2) const;
};

/*!
 * \class sbx_crossover
 *
 * interface to sbx_crossover_impl
 */
template <template <typename> class Chromosome, typename Encoding> class sbx_crossover;
template <template <typename> class Chromosome>
class sbx_crossover<Chromosome,real_encoding> :
    public sbx_crossover_impl<Chromosome,real_encoding>
{
};

/*!
 * \class cycle_crossover_impl
 *
 * implementation class for cycle crossover (CX).  CX works by finding
 * subsets of the indices such that the allele values for the two parents
 * form a complete "cycle" amongst the indices.  it then swaps the alleles
 * between the two parents.
 *
 * CX is useful for permutation problems in which absolute position within
 * the genome is the most important factor in determining fitness.
 */
template <template <typename> class Chromosome, typename Encoding>
class cycle_crossover_impl : public crossover_operator<Chromosome,Encoding>
{
private:
    // disable copying of functional class
    cycle_crossover_impl(const cycle_crossover_impl& that);
    cycle_crossover_impl& operator=(const cycle_crossover_impl& that);

public:
    cycle_crossover_impl();
    virtual ~cycle_crossover_impl();

    virtual void crossover(const Chromosome<Encoding>& p1,
                           const Chromosome<Encoding>& p2,
                           Chromosome<Encoding>& c1,
                           Chromosome<Encoding>& c2) const;
};

/*!
 * \class cycle_crossover
 *
 * interface to cycle_crossover_impl
 */
template <template <typename> class Chromosome, typename Encoding> class cycle_crossover;
template <template <typename> class Chromosome>
class cycle_crossover<Chromosome,permutation_encoding> :
    public cycle_crossover_impl<Chromosome,permutation_encoding>
{
};

/*!
 * \class order_crossover_impl
 *
 * implementation class for order crossover (OX).  OX is useful for
 * permutation problems in which relative ordering within the genome
 * is more important than absolute position.
 */
template <template <typename> class Chromosome, typename Encoding>
class order_crossover_impl : public crossover_operator<Chromosome,Encoding>
{
private:
    // disable copying of functional class
    order_crossover_impl(const order_crossover_impl& that);
    order_crossover_impl& operator=(const order_crossover_impl& that);

public:
    order_crossover_impl();
    virtual ~order_crossover_impl();

    virtual void crossover(const Chromosome<Encoding>& p1,
                           const Chromosome<Encoding>& p2,
                           Chromosome<Encoding>& c1,
                           Chromosome<Encoding>& c2) const;
};

/*!
 * \class order_crossover
 *
 * interface to order_crossover_impl
 */
template <template <typename> class Chromosome, typename Encoding> class order_crossover;
template <template <typename> class Chromosome>
class order_crossover<Chromosome,permutation_encoding> :
    public order_crossover_impl<Chromosome,permutation_encoding>
{
};

/*!
 * \class bit_vector_crossover_operator_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class bit_vector_crossover_operator_factory
{
public:
    static crossover_operator<Chromosome,Encoding>* construct();
};

/*!
 * \class real_crossover_operator_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class real_crossover_operator_factory
{
public:
    static crossover_operator<Chromosome,Encoding>* construct();
};

/*!
 * \class permutation_crossover_operator_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class permutation_crossover_operator_factory
{
public:
    static crossover_operator<Chromosome,Encoding>* construct();
};

/*!
 * \class integer_crossover_operator_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class integer_crossover_operator_factory
{
public:
    static crossover_operator<Chromosome,Encoding>* construct();
};

/*!
 * \class gap_crossover_operator_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class gap_crossover_operator_factory
{
public:
    static crossover_operator<Chromosome,Encoding>* construct();
};

/*!
 * \class gap_crossover_operator_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class gsap_crossover_operator_factory
{
public:
    static crossover_operator<Chromosome,Encoding>* construct();
};

/*!
 * \class crossover_operator_factory
 *
 * empty base class -- specializations handle actual logic
 */
template <template <typename> class Chromosome, typename Encoding> class crossover_operator_factory;

/*!
 * \class crossover_operator_factory
 *
 * creates crossover operators for binary encodings
 */
template <template <typename> class Chromosome>
class crossover_operator_factory<Chromosome,binary_encoding> :
    public bit_vector_crossover_operator_factory<Chromosome,binary_encoding>
{
};

/*!
 * \class crossover_operator_factory
 *
 * creates crossover operators for boolean encodings
 */
template <template <typename> class Chromosome>
class crossover_operator_factory<Chromosome,boolean_encoding> :
    public bit_vector_crossover_operator_factory<Chromosome,boolean_encoding>
{
};

/*!
 * \class crossover_operator_factory
 *
 * creates crossover operators for real encodings
 */
template <template <typename> class Chromosome>
class crossover_operator_factory<Chromosome,real_encoding> :
    public real_crossover_operator_factory<Chromosome,real_encoding>
{
};

/*!
 * \class crossover_operator_factory
 *
 * creates crossover operators for permutation encodings
 */
template <template <typename> class Chromosome>
class crossover_operator_factory<Chromosome,permutation_encoding> :
    public permutation_crossover_operator_factory<Chromosome,permutation_encoding>
{
};

/*!
 * \class crossover_operator_factory
 *
 * creates crossover operators for integer encodings
 */
template <template <typename> class Chromosome>
class crossover_operator_factory<Chromosome,integer_encoding> :
    public integer_crossover_operator_factory<Chromosome,integer_encoding>
{
};

/*!
 * \class crossover_operator_factory
 *
 * creates crossover operators for gap encodings
 */
template <template <typename> class Chromosome>
class crossover_operator_factory<Chromosome,gap_encoding> :
    public gap_crossover_operator_factory<Chromosome,gap_encoding>
{
};

/*!
 * \class crossover_operator_factory
 *
 * creates crossover operators for gap encodings
 */
template <template <typename> class Chromosome>
class crossover_operator_factory<Chromosome,gsap_encoding> :
    public gsap_crossover_operator_factory<Chromosome,gsap_encoding>
{
};

#include "crossover.cpp"

#endif
    
