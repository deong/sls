/**
 * \file repair.h
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _REPAIR_H_
#define _REPAIR_H_

#include "chromosome.h"
#include "neighborhood.h"
#include "encoding.h"

/**
 * \class repair_operator
 * \brief abstract base class for all repair operators 
 */
template <template <typename> class Chromosome, typename Encoding>
class repair_operator
{
private:
    repair_operator(const repair_operator& that);
    repair_operator& operator=(const repair_operator& that);

public:
    repair_operator();
    virtual ~repair_operator();
    virtual void repair(Chromosome<Encoding>& chr,
                        const typename Encoding::ProblemType* prob) const=0;
};

/**
 * \class gap_repair_impl
 * \brief implementation class for gap_repair_nd operator
 */
template <template <typename> class Chromosome, typename Encoding>
class gap_repair_nd_impl : public repair_operator<Chromosome,Encoding>
{
private:
    gap_repair_nd_impl(const gap_repair_nd_impl& that);
    gap_repair_nd_impl& operator=(const gap_repair_nd_impl& that);

public:
    gap_repair_nd_impl();
    virtual ~gap_repair_nd_impl();
    virtual void repair(Chromosome<Encoding>& chr,
                        const typename Encoding::ProblemType* prob) const;

protected:
    neighborhood<Chromosome,Encoding>* m_nf;
};

/**
 * \class gap_repair_nd
 * \brief public interface to gap_repair_nd operator
 */
template <template <typename> class Chromosome, typename Encoding> class gap_repair_nd;
template <template <typename> class Chromosome>
class gap_repair_nd<Chromosome,gap_encoding> :
    public gap_repair_nd_impl<Chromosome,gap_encoding>
{
};

/**
 * \class gsap_repair_impl
 * \brief implementation class for gsap_repair operator
 */
template <template <typename> class Chromosome, typename Encoding>
class gsap_repair_impl : public repair_operator<Chromosome,Encoding>
{
private:
    gsap_repair_impl(const gsap_repair_impl& that);
    gsap_repair_impl& operator=(const gsap_repair_impl& that);

public:
    gsap_repair_impl();
    virtual ~gsap_repair_impl();
    virtual void repair(Chromosome<Encoding>& chr,
                        const typename Encoding::ProblemType* prob) const;

protected:
    neighborhood<Chromosome,Encoding>* m_nf;
};

/**
 * \class gap_repair_nd
 * \brief public interface to gap_repair_nd operator
 */
template <template <typename> class Chromosome, typename Encoding> class gsap_repair;
template <template <typename> class Chromosome>
class gsap_repair<Chromosome,gsap_encoding> :
    public gsap_repair_impl<Chromosome,gsap_encoding>
{
};

/**
 * \class gap_repair_sd_impl
 * \brief implementation class for gap_repair_sd operator
 */
template <template <typename> class Chromosome, typename Encoding>
class gap_repair_sd_impl : public repair_operator<Chromosome,Encoding>
{
private:
    gap_repair_sd_impl(const gap_repair_sd_impl& that);
    gap_repair_sd_impl& operator=(const gap_repair_sd_impl& that);

public:
    gap_repair_sd_impl();
    virtual ~gap_repair_sd_impl();
    virtual void repair(Chromosome<Encoding>& chr,
                        const typename Encoding::ProblemType* prob) const;

protected:
    neighborhood<Chromosome,Encoding>* m_nf;
};

/**
 * \class gap_repair_sd
 * \brief public interface to gap_repair_sd operator
 */
template <template <typename> class Chromosome, typename Encoding> class gap_repair_sd;
template <template <typename> class Chromosome>
class gap_repair_sd<Chromosome,gap_encoding> :
    public gap_repair_sd_impl<Chromosome,gap_encoding>
{
};

/*!
 * \class bit_vector_repair_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class bit_vector_repair_factory : public factory
{
public:
    repair_operator<Chromosome,Encoding>* construct();
};

/*!
 * \class real_repair_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class real_repair_factory : public factory
{
public:
    repair_operator<Chromosome,Encoding>* construct();
};

/*!
 * \class permutation_repair_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class permutation_repair_factory : public factory
{
public:
    repair_operator<Chromosome,Encoding>* construct();
};

/*!
 * \class integer_repair_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class integer_repair_factory : public factory
{
public:
    repair_operator<Chromosome,Encoding>* construct();
};

/*!
 * \class gap_repair_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class gap_repair_factory : public factory
{
public:
    repair_operator<Chromosome,Encoding>* construct();
};

/*!
 * \class gsap_repair_factory
 */
template <template <typename> class Chromosome, typename Encoding>
class gsap_repair_factory : public factory
{
public:
    repair_operator<Chromosome,Encoding>* construct();
};

/*!
 * \class repair_factory
 *
 * empty base class -- specializations handle actual logic
 */
template <template <typename> class Chromosome, typename Encoding> class repair_factory;

/*!
 * \class repair_factory
 *
 * creates repair operators for binary encodings
 */
template <template <typename> class Chromosome>
class repair_factory<Chromosome,binary_encoding> :
    public bit_vector_repair_factory<Chromosome,binary_encoding>
{
};

/*!
 * \class repair_factory
 *
 * creates repair operators for boolean encodings
 */
template <template <typename> class Chromosome>
class repair_factory<Chromosome,boolean_encoding> :
    public bit_vector_repair_factory<Chromosome,boolean_encoding>
{
};

/*!
 * \class repair_factory
 *
 * creates repair operators for real encodings
 */
template <template <typename> class Chromosome>
class repair_factory<Chromosome,real_encoding> :
    public real_repair_factory<Chromosome,real_encoding>
{
};

/*!
 * \class repair_factory
 *
 * creates repair operators for permutation encodings
 */
template <template <typename> class Chromosome>
class repair_factory<Chromosome,permutation_encoding> :
    public permutation_repair_factory<Chromosome,permutation_encoding>
{
};

/*!
 * \class repair_factory
 *
 * creates repair operators for integer encodings
 */
template <template <typename> class Chromosome>
class repair_factory<Chromosome,integer_encoding> :
    public integer_repair_factory<Chromosome,integer_encoding>
{
};

/*!
 * \class repair_factory
 *
 * creates repair operators for gap encodings
 */
template <template <typename> class Chromosome>
class repair_factory<Chromosome,gap_encoding> :
    public gap_repair_factory<Chromosome,gap_encoding>
{
};

/*!
 * \class repair_factory
 *
 * creates repair operators for gap encodings
 */
template <template <typename> class Chromosome>
class repair_factory<Chromosome,gsap_encoding> :
    public gsap_repair_factory<Chromosome,gsap_encoding>
{
};

#include "repair.cpp"

#endif
