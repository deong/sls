/*!
 * \file mutation.h
 *
 * evolutionary algorithm mutation operators
 *
 * Deon Garrett
 * deong@acm.org
 */

#ifndef _MUTATION_H_
#define _MUTATION_H_

#include "chromosome.h"
#include "encoding.h"

/*!
 * \class mutation_operator
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class mutation_operator
{
private:
    // disable copying of functional classes
    mutation_operator(const mutation_operator& that);
    mutation_operator& operator=(const mutation_operator& that);

public:
    mutation_operator();
    virtual ~mutation_operator();

    virtual void initialize();
    virtual void mutate(Chromosome<Encoding>& sol) const = 0;
};

/*!
 * \class bitwise_mutation_impl
 *
 * implementation for bitflip mutation operator
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class bitwise_mutation_impl : public mutation_operator<Chromosome,Encoding>
{
private:
    // disable copying for functional class
    bitwise_mutation_impl(const bitwise_mutation_impl& that);
    bitwise_mutation_impl& operator=(const bitwise_mutation_impl& that);

protected:
    double m_rate;
    
public:
    bitwise_mutation_impl();
    virtual ~bitwise_mutation_impl();

    virtual void initialize();
    virtual void mutate(Chromosome<Encoding>& sol) const;
};

/*!
 * \class bitwise_mutation
 *
 * interface to bitwise_mutation_impl
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding> class bitwise_mutation;
template <template <typename> class Chromosome>
class bitwise_mutation<Chromosome,boolean_encoding> :
    public bitwise_mutation_impl<Chromosome,boolean_encoding>
{
};

/*!
 * \class bitwise_mutation
 *
 * interface to bitwise_mutation_impl
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome>
class bitwise_mutation<Chromosome,binary_encoding> :
    public bitwise_mutation_impl<Chromosome,binary_encoding>
{
};

/*!
 * \class swap_mutation_impl
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class swap_mutation_impl : public mutation_operator<Chromosome,Encoding>
{
private:
    // disable copying of functional classes
    swap_mutation_impl(const swap_mutation_impl& that);
    swap_mutation_impl& operator=(const swap_mutation_impl& that);

protected:
    double m_rate;
    
public:
    swap_mutation_impl();
    virtual ~swap_mutation_impl();

    virtual void initialize();
    virtual void mutate(Chromosome<Encoding>& sol) const;
};

/*!
 * \class swap_mutation
 *
 * interface to swap_mutation_impl
 * 
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding> class swap_mutation;
template <template <typename> class Chromosome>
class swap_mutation<Chromosome,permutation_encoding> :
    public swap_mutation_impl<Chromosome,permutation_encoding>
{
};

/*!
 * \class swap_mutation
 *
 * interface to swap_mutation_impl
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome>
class swap_mutation<Chromosome,gap_encoding> :
    public swap_mutation_impl<Chromosome,gap_encoding>
{
};

/*!
 * \class swap_mutation
 *
 * interface to swap_mutation_impl
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome>
class swap_mutation<Chromosome,gsap_encoding> :
    public swap_mutation_impl<Chromosome,gsap_encoding>
{
};

/*!
 * \class gaussian_mutation_impl
 *
 * mutate real valued chromosome via gaussian randomness
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class gaussian_mutation_impl : public mutation_operator<Chromosome,Encoding>
{
private:
    // disable copying of functional class
    gaussian_mutation_impl(const gaussian_mutation_impl& that);
    gaussian_mutation_impl& operator=(const gaussian_mutation_impl& that);

protected:
    double m_rate;
    double m_mu;
    double m_sigma;

public:
    gaussian_mutation_impl();
    virtual ~gaussian_mutation_impl();

    virtual void initialize();
    virtual void mutate(Chromosome<Encoding>& sol) const;
};

/*!
 * \class gaussian_mutation
 *
 * interface to gaussian_mutation_impl
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding> class gaussian_mutation;
template <template <typename> class Chromosome>
class gaussian_mutation<Chromosome,real_encoding> :
    public gaussian_mutation_impl<Chromosome,real_encoding>
{
};

/*!
 * \class polynomial_mutation_impl
 *
 * polynomial mutation for real valued chromosomes
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class polynomial_mutation_impl : public mutation_operator<Chromosome,Encoding>
{
private:
    // disable copying of functional class
    polynomial_mutation_impl(const polynomial_mutation_impl& that);
    polynomial_mutation_impl& operator=(const polynomial_mutation_impl& that);

protected:
    double m_rate;
    double m_eta;

public:
    polynomial_mutation_impl();
    virtual ~polynomial_mutation_impl();

    virtual void initialize();
    virtual void mutate(Chromosome<Encoding>& sol) const;
};

/*!
 * \class polynomial_mutation
 *
 * interface to polynomial_mutation_impl
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding> class polynomial_mutation;
template <template <typename> class Chromosome>
class polynomial_mutation<Chromosome,real_encoding> :
    public polynomial_mutation_impl<Chromosome,real_encoding>
{
};

/*!
 * \class shift_mutation_impl
 *
 * for gap chromosomes, shift a job to another agent
 * 
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class shift_mutation_impl : public mutation_operator<Chromosome,Encoding>
{
private:
    // disable copying
    shift_mutation_impl(const shift_mutation_impl& that);
    shift_mutation_impl& operator=(const shift_mutation_impl& that);

protected:
    double m_rate;

public:
    shift_mutation_impl();
    virtual ~shift_mutation_impl();

    virtual void initialize();
    virtual void mutate(Chromosome<Encoding>& sol) const;
};

/*!
 * \class shift_mutation
 *
 * interface to shift_mutation_impl
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding> class shift_mutation;
template <template <typename> class Chromosome>
class shift_mutation<Chromosome,gap_encoding> :
    public shift_mutation_impl<Chromosome,gap_encoding>
{
};

/*!
 * \class shift_mutation
 *
 * interface to shift_mutation_impl
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome>
class shift_mutation<Chromosome,gsap_encoding> :
    public shift_mutation_impl<Chromosome,gsap_encoding>
{
};

/*!
 * \class sss_mutation_impl
 *
 * for gap chromosomes, randomly choose between shift and swap mutations
 *
 * \author deong
 * \date 06/11/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class sss_mutation_impl : public mutation_operator<Chromosome,Encoding>
{
private:
    // disable copying
    sss_mutation_impl(const sss_mutation_impl& that);
    sss_mutation_impl& operator=(const sss_mutation_impl& that);

protected:
    double m_rate;

public:
    sss_mutation_impl();
    virtual ~sss_mutation_impl();
    virtual void initialize();
    virtual void mutate(Chromosome<Encoding>& sol) const;
};

/*!
 * \class sss_mutation
 *
 * interface to sss_mutation_impl
 *
 * \author deong
 * \date 06/11/2007
 */
template <template <typename> class Chromosome, typename Encoding> class sss_mutation;
template <template <typename> class Chromosome>
class sss_mutation<Chromosome,gap_encoding> :
    public sss_mutation_impl<Chromosome,gap_encoding>
{
};

/*!
 * \class sss_mutation
 *
 * interface to sss_mutation_impl
 *
 * \author deong
 * \date 06/11/2007
 */
template <template <typename> class Chromosome>
class sss_mutation<Chromosome,gsap_encoding> :
    public sss_mutation_impl<Chromosome,gsap_encoding>
{
};

/*!
 * \class bit_mutation_operator_factory
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class bit_mutation_operator_factory
{
public:
    static mutation_operator<Chromosome,Encoding>* construct();
};

/*!
 * \class permutation_mutation_operator_factory
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class permutation_mutation_operator_factory
{
public:
    static mutation_operator<Chromosome,Encoding>* construct();
};

/*!
 * \class real_mutation_operator_factory
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class real_mutation_operator_factory
{
public:
    static mutation_operator<Chromosome,Encoding>* construct();
};

/*!
 * \class integer_mutation_operator_factory
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class integer_mutation_operator_factory
{
public:
    static mutation_operator<Chromosome,Encoding>* construct();
};

/*!
 * \class gap_mutation_operator_factory
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class gap_mutation_operator_factory
{
public:
    static mutation_operator<Chromosome,Encoding>* construct();
};

/*!
 * \class gap_mutation_operator_factory
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding>
class gsap_mutation_operator_factory
{
public:
    static mutation_operator<Chromosome,Encoding>* construct();
};

/*!
 * \class mutation_operator_factory
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome, typename Encoding> class mutation_operator_factory;

/*!
 * \class mutation_operator_factory
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome>
class mutation_operator_factory<Chromosome,boolean_encoding> :
    public bit_mutation_operator_factory<Chromosome,boolean_encoding>
{
};

/*!
 * \class mutation_operator_factory
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome>
class mutation_operator_factory<Chromosome,binary_encoding> :
    public bit_mutation_operator_factory<Chromosome,binary_encoding>
{
};

/*!
 * \class mutation_operator_factory
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome>
class mutation_operator_factory<Chromosome,real_encoding> :
    public real_mutation_operator_factory<Chromosome,real_encoding>
{
};

/*!
 * \class mutation_operator_factory
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome>
class mutation_operator_factory<Chromosome,permutation_encoding> :
    public permutation_mutation_operator_factory<Chromosome,permutation_encoding>
{
};

/*!
 * \class mutation_operator_factory
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome>
class mutation_operator_factory<Chromosome,integer_encoding> :
    public integer_mutation_operator_factory<Chromosome,integer_encoding>
{
};

/*!
 * \class mutation_operator_factory
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome>
class mutation_operator_factory<Chromosome,gap_encoding> :
    public gap_mutation_operator_factory<Chromosome,gap_encoding>
{
};

/*!
 * \class mutation_operator_factory
 *
 * \author deong
 * \date 05/10/2007
 */
template <template <typename> class Chromosome>
class mutation_operator_factory<Chromosome,gsap_encoding> :
    public gsap_mutation_operator_factory<Chromosome,gsap_encoding>
{
};

#include "mutation.cpp"

#endif
