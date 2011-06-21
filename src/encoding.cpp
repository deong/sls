/*!
 * \file encoding.cpp
 *
 * for stochastic local search, we need the notion of a candidate
 * solution to the problem.  the solution may be encoded in some form
 * (such as bit strings in genetic algorithms).  these classes provide
 * the representations for several different types of encodings
 *
 * Deon Garrett
 * deong@acm.org
 */

#include <iostream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <cassert>
#include "encoding.h"
#include "mtrandom.h"
#include "configuration.h"
#include "keywords.h"
#include "utilities.h"

using namespace std;

// initialize static instances
vector<int> binary_encoding::m_bpp;
bool binary_encoding::m_gray;
unsigned int binary_encoding::m_len;
vector<pair<double,double> > numeric_encoding::m_range;
vector<vector<int> > integer_encoding::m_legal_values;
unsigned int gap_encoding::m_agents;
unsigned int gap_encoding::m_tasks;
vector<double> gap_encoding::m_alpha;
vector<gsap_problem::task_element> gsap_encoding::elements;
vector<vector<int> > gsap_encoding::agents_for_task;
int gsap_encoding::m_unass_pen;
int gsap_encoding::m_cap_pen;

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename G, typename P>
encoding<G,P>::encoding()
{
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename G, typename P>
encoding<G,P>::encoding(const problem* p) :
    m_phenotype(p->dimensions(), 0)
{
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename G, typename P>
encoding<G,P>::~encoding()
{
}

/*!
 * \brief assignment operator
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename G, typename P>
encoding<G,P>& encoding<G,P>::operator=(const encoding<G,P>& that)
{
    m_genotype = that.m_genotype;
    m_phenotype = that.m_phenotype;
    return *this;
}

/*!
 * \brief equality operator
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename G, typename P>
bool encoding<G,P>::operator==(const encoding<G,P>& that)
{
    return this->m_genotype == that.m_genotype;
}

/*!
 * \brief return a given element of the encoded genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename G, typename P>
inline G& encoding<G,P>::operator[](unsigned int index)
{
    return this->m_genotype[index];
}

/*!
 * \brief return a given element of the encoded genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename G, typename P>
inline const G& encoding<G,P>::operator[](unsigned int index) const
{
    return this->m_genotype[index];
}

/*!
 * \brief return the genotype vector
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename G, typename P>
vector<G>& encoding<G,P>::genotype()
{
    return m_genotype;
}

/*!
 * \brief return the genotype vector
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename G, typename P>
const vector<G>& encoding<G,P>::genotype() const
{
    return this->m_genotype;
}

/*!
 * \brief return the phenotype vector
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename G, typename P>
vector<P>& encoding<G,P>::phenotype()
{
    return this->m_phenotype;
}

/*!
 * \brief return the phenotype vector
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename G, typename P>
const vector<P>& encoding<G,P>::phenotype() const
{
    return this->m_phenotype;
}

/*!
 * \brief return the length of the genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename G, typename P>
inline unsigned int encoding<G,P>::length() const
{
    return this->m_genotype.size();
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename P>
bit_vector_encoding<P>::bit_vector_encoding()
{
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename P>
bit_vector_encoding<P>::bit_vector_encoding(const problem* p) :
    encoding<int,P>::encoding(p)
{
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename P>
bit_vector_encoding<P>::~bit_vector_encoding()
{
}

/*!
 * \brief return an iterator to the front of the genotype vector
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename P>
typename bit_vector_encoding<P>::iterator bit_vector_encoding<P>::begin()
{
    return this->m_genotype.begin();
}

/*!
 * \brief return an iterator to the front of the genotype vector
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename P>
typename bit_vector_encoding<P>::const_iterator bit_vector_encoding<P>::begin() const
{
    return this->m_genotype.begin();
}

/*!
 * \brief return an iterator to the end of the genotype vector
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename P>
typename bit_vector_encoding<P>::iterator bit_vector_encoding<P>::end()
{
    return this->m_genotype.end();
}

/*!
 * \brief return an iterator to the end of the genotype vector
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename P>
typename bit_vector_encoding<P>::const_iterator bit_vector_encoding<P>::end() const
{
    return this->m_genotype.end();
}

/*!
 * \brief randomize the encoded solution
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename P>
void bit_vector_encoding<P>::randomize()
{
    mtrandom mt;
    for(unsigned int i=0; i<this->m_genotype.size(); i++)
    {
        this->m_genotype[i] = mt.random(2);
    }
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
boolean_encoding::boolean_encoding()
{
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
boolean_encoding::boolean_encoding(const bit_string_problem* p) :
    bit_vector_encoding<int>::bit_vector_encoding(p)
{
    m_genotype.resize(p->dimensions());
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
boolean_encoding::~boolean_encoding()
{
}

/*!
 * \brief initialize encoding parameters
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void boolean_encoding::initialize_parameters(const bit_string_problem* p)
{
}

/*!
 * \brief clear encoding parameters
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void boolean_encoding::clear_parameters()
{
}

/*!
 * \brief encode an integer into a boolean genotype vector
 *
 * \author deong@acm.org
 * \date 10/10/2008
 */
inline void boolean_encoding::encode(unsigned int p)
{
    unsigned int num_bits=this->length();
    
    // convert the integer to binary
    vector<int> result(num_bits);
    for(unsigned int b=0; b<num_bits; ++b)
    {
        result[b]=p%2;
        p/=2;
    }
    
    // insert the bits into the correct spot in the encoding
    copy(result.begin(),result.end(),m_genotype.begin());
}

/*!
 * \brief decode the genotype into a phenotype vector
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
inline void boolean_encoding::decode()
{
    for(unsigned int i=0; i<m_genotype.size(); i++)
    {
        m_phenotype[i] = m_genotype[i];
    }
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
numeric_encoding::numeric_encoding()
{
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
numeric_encoding::~numeric_encoding()
{
}

/*!
 * \brief initialize the encoding parameters
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void numeric_parameters::initialize_parameters(const numeric_problem* p)
{
    unsigned int d = p->dimensions();
    for(unsigned int i=0; i<d; i++)
    {
        pair<double,double> range = p->parameter_range(i);
        numeric_encoding::m_range.push_back(range);
    }
}

/*!
 * \brief erase the encoding parameters
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void numeric_parameters::cleanup()
{
    numeric_encoding::m_range.clear();
}

/*!
 * \brief return the range of valid values for an allele
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
const pair<double,double>& numeric_encoding::parameter_range(int pnum)
{
    return m_range[pnum];
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
binary_encoding::binary_encoding()
{
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
binary_encoding::binary_encoding(const numeric_problem* p) :
    bit_vector_encoding<double>::bit_vector_encoding(p)
{
    m_genotype.resize(m_len);
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
binary_encoding::~binary_encoding()
{
}

/*!
 * \brief initialize the encoding parameters
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void binary_encoding::initialize_parameters(const numeric_problem* p)
{
    numeric_parameters::initialize_parameters(p);
    binary_parameters::initialize_parameters(p);
}

/*!
 * \brief clear the encoding parameters
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void binary_encoding::clear_parameters()
{
    numeric_parameters::cleanup();
    binary_parameters::cleanup();
}

inline void binary_encoding::encode(const vector<double>& params)
{
    unsigned int start=0;

    for(unsigned int param=0; param<params.size(); ++param)
    {
        unsigned int num_bits = m_bpp[param];
        
        // map the double onto the appropriate integer range
        pair<double,double> prange=m_range[param];
        double range=prange.second-prange.first;
        double max_bit_val=pow(2.0,static_cast<double>(num_bits))-1;
        int int_val=static_cast<int>((params[param]-prange.first)*max_bit_val/range+0.5);
        
        // convert the integer to binary
        vector<int> result(m_bpp[param]);
        for(unsigned int b=0; b<num_bits; ++b)
        {
            result[b]=int_val%2;
            int_val/=2;
        }

        if(m_gray)
        {
            for(unsigned int b=0; b<num_bits-1; ++b)
            {
                result[b]=!(result[b]==result[b+1]);
            }
        }
        
        // insert the bits into the correct spot in the encoding
        copy(result.begin(),result.end(),m_genotype.begin()+start);
        start+=num_bits;
    }
}

/*!
 * \brief decode the genotype into a phenotype vector
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
inline void binary_encoding::decode()
{
    unsigned int start = 0;
    
    // for each parameter
    for(unsigned int param=0; param<m_bpp.size(); param++)
    {
        unsigned int num_bits = m_bpp[param];
        unsigned int intval = 0;
        if(m_gray)
        {
            // convert from gray to binary
            vector<int> binary(num_bits);
            binary[num_bits-1] = m_genotype[start+num_bits-1];
            intval = binary[num_bits-1];
            for(int i=num_bits-2; i>=0; i--)
            {
                binary[i] = !(binary[i+1] == m_genotype[start+i]);
                intval += intval + binary[i];
            }
        }
        else
        {
            // convert from binary encoding to integer
            for(int i=num_bits-1; i>=0; i--)
            {
                intval += intval + m_genotype[start+i];
            }
        }
        
        // convert from integer to double in the appropriate range
        pair<double,double> prange = m_range[param];
        double range = prange.second - prange.first;
        double m = range / (pow(2.0,double(num_bits)) - 1.0);
        m_phenotype[param] = m * double(intval) + prange.first;
        
        start += num_bits;
    }
}

/*!
 * \brief initialize the encoding parameters
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void binary_parameters::initialize_parameters(const numeric_problem* p)
{
    unsigned int dim = p->dimensions();
    unsigned int totallen = 0;
    if(configuration::keyword_exists("parameter_length"))
    {
        unsigned int plen;
        configuration::unsigned_integer_parameter("parameter_length", plen, true);
        for(unsigned int i=0; i<dim; i++)
        {
            binary_encoding::m_bpp.push_back(plen);
        }
        totallen = plen * dim;
    }
    else
    {
        ostringstream s;
        for(unsigned int i=0; i<dim; i++)
        {
            s << "parameter_" << i << "_length";
            string pname = s.str();
            unsigned int plen = 0;
            configuration::unsigned_integer_parameter(pname, plen, true);
            totallen += plen;
            binary_encoding::m_bpp.push_back(plen);
        }
    }

    binary_encoding::m_gray = false;
    string enc;
    configuration::string_parameter("parameter_encoding", enc, false);
    if(enc == "gray")
    {
        binary_encoding::m_gray = true;
    }
    else if(enc == "binary")
    {
        binary_encoding::m_gray = false;
    }
    else
    {
        cerr << "illegal value for parameter_encoding: " << enc << endl;
        cerr << "must be binary or gray" << endl;
        exit(1);
    }

    binary_encoding::m_len = totallen;
}

/*!
 * \brief erase the encoding parameters
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void binary_parameters::cleanup()
{
    binary_encoding::m_bpp.clear();
    binary_encoding::m_len = 0;
    binary_encoding::m_gray = false;
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
real_encoding::real_encoding()
{
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
real_encoding::real_encoding(const numeric_problem* p) :
    encoding<double,double>::encoding(p)
{
    m_genotype.resize(p->dimensions());
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
real_encoding::~real_encoding()
{
}

/*!
 * \brief return an iterator to the beginning of the genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
real_encoding::iterator real_encoding::begin()
{
    return this->m_genotype.begin();
}

/*!
 * \brief return an iterator to the beginning of the genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
real_encoding::const_iterator real_encoding::begin() const
{
    return this->m_genotype.begin();
}

/*!
 * \brief return an iterator to the end of the genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
real_encoding::iterator real_encoding::end()
{
    return this->m_genotype.end();
}

/*!
 * \brief return an iterator to the end of the genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
real_encoding::const_iterator real_encoding::end() const
{
    return this->m_genotype.end();
}

/*!
 * \brief initialize the encoding parameters
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void real_encoding::initialize_parameters(const numeric_problem* p)
{
    numeric_parameters::initialize_parameters(p);
}

/*!
 * \brief clear the encoding parameters
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void real_encoding::clear_parameters()
{
    numeric_parameters::cleanup();
}

/*!
 * \brief randomize the encoded solution
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void real_encoding::randomize()
{
    mtrandom mt;
    for(unsigned int i=0; i<m_genotype.size(); i++)
    {
        pair<double,double> rng = m_range[i];
        m_genotype[i] = mt.random(rng.first, rng.second);
    }
}

/*!
 * \brief decode the genotype into a phenotype vector
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
inline void real_encoding::decode()
{
    for(unsigned int i=0; i<m_genotype.size(); i++)
    {
        m_phenotype[i] = m_genotype[i];
    }
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
permutation_encoding::permutation_encoding()
{
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
permutation_encoding::permutation_encoding(const permutation_problem* p) :
    encoding<int,int>::encoding(p)
{
    m_genotype.resize(p->dimensions());
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
permutation_encoding::~permutation_encoding()
{
}

/*!
 * \brief return an iterator to the beginning of the genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
permutation_encoding::iterator permutation_encoding::begin()
{
    return this->m_genotype.begin();
}

/*!
 * \brief return an iterator to the beginning of the genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
permutation_encoding::const_iterator permutation_encoding::begin() const
{
    return this->m_genotype.begin();
}

/*!
 * \brief return an iterator to the end of the genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
permutation_encoding::iterator permutation_encoding::end()
{
    return this->m_genotype.end();
}

/*!
 * \brief return an iterator to the end of the genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
permutation_encoding::const_iterator permutation_encoding::end() const
{
    return this->m_genotype.end();
}

/*!
 * \brief initialize the encoding parameters
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void permutation_encoding::initialize_parameters(const permutation_problem* p)
{
}

/*!
 * \brief clear the encoding parameters
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void permutation_encoding::clear_parameters()
{
}

/*!
 * \brief randomize the encoded solution
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void permutation_encoding::randomize()
{
    mtrandom mt;
    m_genotype = mt.permutation(m_genotype.size());
}

/*!
 * \brief decode the genotype into a phenotype vector
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
inline void permutation_encoding::decode()
{
    for(unsigned int i=0; i<m_genotype.size(); i++)
    {
        m_phenotype[i] = m_genotype[i];
    }
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
integer_encoding::integer_encoding()
{
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
integer_encoding::integer_encoding(const integer_problem* p) :
    encoding<int,int>::encoding(p)
{
    m_genotype.resize(p->dimensions());
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
integer_encoding::~integer_encoding()
{
}

/*!
 * \brief return an iterator to the beginning of the genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
integer_encoding::iterator integer_encoding::begin()
{
    return this->m_genotype.begin();
}

/*!
 * \brief return an iterator to the beginning of the genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
integer_encoding::const_iterator integer_encoding::begin() const
{
    return this->m_genotype.begin();
}

/*!
 * \brief return an iterator to the end of the genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
integer_encoding::iterator integer_encoding::end()
{
    return this->m_genotype.end();
}

/*!
 * \brief return an iterator to the end of the genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
integer_encoding::const_iterator integer_encoding::end() const
{
    return this->m_genotype.end();
}

/*!
 * \brief initialize the encoding parameters
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void integer_encoding::initialize_parameters(const integer_problem* p)
{
    for(unsigned int i=0; i<p->dimensions(); i++)
    {
        vector<int> vals;
        p->legal_values(i, vals);
        m_legal_values.push_back(vals);
    }
}

/*!
 * \brief clear the encoding parameters
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void integer_encoding::clear_parameters()
{
    for(unsigned int i=0; i<m_legal_values.size(); i++)
    {
        m_legal_values[i].clear();
    }
    m_legal_values.clear();
}

/*!
 * \brief return a vector of valid allele values
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
const vector<int>& integer_encoding::legal_values(unsigned int index) const
{
    return m_legal_values[index];
}

/*!
 * \brief randomize the encoded solution
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void integer_encoding::randomize()
{
    mtrandom mt;
    for(unsigned int i=0; i<m_genotype.size(); i++)
    {
        m_genotype[i] = m_legal_values[i][mt.random(m_legal_values[i].size())];
    }
}

/*!
 * \brief decode the genotype into a phenotype vector
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void integer_encoding::decode()
{
    for(unsigned int i=0; i<m_genotype.size(); i++)
    {
        m_phenotype[i] = m_genotype[i];
    }
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
gap_encoding::gap_encoding()
{
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
gap_encoding::gap_encoding(const gap_problem* p) :
    encoding<int,int>::encoding(p)
{
    m_genotype.resize(p->dimensions());
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
gap_encoding::~gap_encoding()
{
}

/*!
 * \brief return the number of agents in the problem
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
unsigned int gap_encoding::agents() const
{
    return m_agents;
}

/*!
 * \brief return the number of tasks in the problem
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
unsigned int gap_encoding::tasks() const
{
    return m_tasks;
}

/*!
 * \brief return an iterator to the beginning of the genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
gap_encoding::iterator gap_encoding::begin()
{
    return this->m_genotype.begin();
}

/*!
 * \brief return an iterator to the beginning of the genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
gap_encoding::const_iterator gap_encoding::begin() const
{
    return this->m_genotype.begin();
}

/*!
 * \brief return an iterator to the end of the genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
gap_encoding::iterator gap_encoding::end()
{
    return this->m_genotype.end();
}

/*!
 * \brief return an iterator to the end of the genotype
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
gap_encoding::const_iterator gap_encoding::end() const
{
    return this->m_genotype.end();
}

/*!
 * \brief initialize the encoding parameters
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void gap_encoding::initialize_parameters(const gap_problem* p)
{
    m_agents = p->agents;
    m_tasks = p->tasks;

    double pen;
    configuration::double_parameter(keywords::PENALIZATION_FACTOR, pen, true);
    m_alpha.assign(p->dimensions(), pen);
}

/*!
 * \brief clear the encoding parameters
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void gap_encoding::clear_parameters()
{
    m_agents = m_tasks = 0;
}

/*!
 * \brief randomize the encoded solution
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void gap_encoding::randomize()
{
    mtrandom mt;
    for(unsigned int i=0; i<m_genotype.size(); i++)
    {
        m_genotype[i] = mt.random(m_agents);
    }
}

/*!
 * \brief decode the genotype into a phenotype vector
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
void gap_encoding::decode()
{
    for(unsigned int i=0; i<m_genotype.size(); i++)
    {
        m_phenotype[i] = m_genotype[i];
    }
}

/**
 * \brief constructor
 *
 * \author deong
 * \date 01/22/2009
 */
gsap_encoding::gsap_encoding()
{
}

/**
 * \brief constructor
 *
 * \author deong
 * \date 01/22/2009
 */
gsap_encoding::gsap_encoding(const ProblemType* p) :
    encoding<int,int>::encoding(p)
{
    m_genotype.resize(p->dimensions());
}

/**
 * \brief destructor
 *
 * \author deong
 * \date 01/22/2009
 */
gsap_encoding::~gsap_encoding()
{
}

/**
 * \brief return iterator to beginning of encoding
 *
 * \author deong
 * \date 01/22/2009
 */
gsap_encoding::iterator gsap_encoding::begin()
{
    return m_genotype.begin();
}

/**
 * \brief return iterator to end of encoding
 *
 * \author deong
 * \date 01/22/2009
 */
gsap_encoding::iterator gsap_encoding::end()
{
    return m_genotype.end();
}

/**
 * \brief return iterator to beginning of encoding
 *
 * \author deong
 * \date 01/22/2009
 */
gsap_encoding::const_iterator gsap_encoding::begin() const
{
    return m_genotype.begin();
}

/**
 * \brief return iterator to end of encoding
 *
 * \author deong
 * \date 01/22/2009
 */
gsap_encoding::const_iterator gsap_encoding::end() const
{
    return m_genotype.end();
}

/**
 * \brief randomize solution
 *
 * \todo make randomization attempt to honor constraints
 * 
 * \author deong
 * \date 01/22/2009
 */
void gsap_encoding::randomize()
{
    mtrandom mt;
    for(unsigned int i=0; i<m_genotype.size(); ++i)
    {
        int taskid = elements[i].task;
        m_genotype[i] = agents_for_task[taskid][mt.random(mt.random(agents_for_task[taskid].size()))];
    }
}

/**
 * \brief trivially decode the solution
 *
 * \author deong
 * \date 01/22/2009
 */
void gsap_encoding::decode()
{
    copy(m_genotype.begin(),m_genotype.end(),m_phenotype.begin());
}

/**
 * \brief initialize encoding parameters
 *
 * \author deong
 * \date 01/22/2009
 */
void gsap_encoding::initialize_parameters(const ProblemType* p)
{
    configuration::integer_parameter(keywords::UNASSIGNED_PENALTY,m_unass_pen,true);
    configuration::integer_parameter(keywords::CAPACITY_PENALTY,m_cap_pen,true);
    elements=p->get_elements();
    agents_for_task = p->get_agent_task_map();
}

/**
 * \brief clear encoding parameters
 *
 * \author deong
 * \date 01/22/2009
 */
void gsap_encoding::clear_parameters()
{
    elements.clear();
    agents_for_task.clear();
}

/*!
 * \brief print the encoded data
 *
 * \author deong
 * \date 05/09/2007
 */
ostream& operator<<(ostream& s, const boolean_encoding& e)
{
    for(unsigned int i=0; i<e.m_genotype.size(); i++)
    {
        s << e.m_genotype[i];
    }
    
    return s;
}

/*!
 * \brief read in an encoded chromosome
 *
 * \author deong@acm.org
 * \date 10/14/2008
 */
istream& operator>>(istream& s, boolean_encoding& e)
{
    for(unsigned int i=0; i<e.m_genotype.size(); i++)
    {
        char c;
        s >> c;
        e.m_genotype[i]=(c=='0'?0:1);
    }
    
    return s;
}

/*!
 * \brief print the encoded data
 *
 * \author deong
 * \date 05/09/2007
 */
ostream& operator<<(ostream& s, const binary_encoding& e)
{
    for(unsigned int i=0; i<e.m_genotype.size(); i++)
    {
        s << e.m_genotype[i];
    }
    
    return s;
}

/*!
 * \brief read in an encoded value
 *
 * \author deong@acm.org
 * \date 10/14/2008
 */
istream& operator>>(istream& s, binary_encoding& e)
{
    for(unsigned int i=0; i<e.m_genotype.size(); i++)
    {
        char c;
        s >> c;
        e.m_genotype[i]=(c=='0'?0:1);
    }
    
    return s;
}


/*!
 * \brief print the encoded data
 *
 * \author deong
 * \date 05/09/2007
 */
ostream& operator<<(ostream& s, const real_encoding& e)
{
    for(unsigned int i=0; i<e.m_genotype.size()-1; i++)
    {
        s << e.m_genotype[i] << " ";
    }
    s << e.m_genotype[e.m_genotype.size()-1];
    
    return s;
}

/*!
 * \brief read in an encoded value
 *
 * \author deong@acm.org
 * \date 10/14/2008
 */
istream& operator>>(istream& s, real_encoding& e)
{
    for(unsigned int i=0; i<e.m_genotype.size(); i++)
    {
        s >> e.m_genotype[i];
    }

    return s;
}

/*!
 * \brief print the encoded data
 *
 * \author deong
 * \date 05/09/2007
 */
ostream& operator<<(ostream& s, const permutation_encoding& e)
{
    for(unsigned int i=0; i<e.m_genotype.size()-1; i++)
    {
        s << e.m_genotype[i] << " ";
    }
    s << e.m_genotype[e.m_genotype.size()-1];
    
    return s;
}

/*!
 * \brief read in an encoded value
 *
 * \author deong@acm.org
 * \date 10/14/2008
 */
istream& operator>>(istream& s, permutation_encoding& e)
{
    for(unsigned int i=0; i<e.m_genotype.size(); i++)
    {
        s >> e.m_genotype[i];
    }

    return s;
}

/*!
 * \brief print the encoded data
 *
 * \author deong
 * \date 05/09/2007
 */
ostream& operator<<(ostream& s, const integer_encoding& e)
{
    for(unsigned int i=0; i<e.m_genotype.size()-1; i++)
    {
        s << e.m_genotype[i] << " ";
    }
    s << e.m_genotype[e.m_genotype.size()-1];

    return s;
}

/*!
 * \brief read in an encoded value
 *
 * \author deong@acm.org
 * \date 10/14/2008
 */
istream& operator>>(istream& s, integer_encoding& e)
{
    for(unsigned int i=0; i<e.m_genotype.size(); ++i)
    {
        s >> e.m_genotype[i];
    }
    return s;
}

/*!
 * \brief print the encoded data
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
ostream& operator<<(ostream& s, const gap_encoding& e)
{
    for(unsigned int i=0; i<e.m_genotype.size()-1; i++)
    {
        s << e.m_genotype[i] << " ";
    }
    s << e.m_genotype[e.m_genotype.size()-1];

    return s;
}

/*!
 * \brief read in an encoded value
 *
 * \author deong@acm.org
 * \date 10/14/2008
 */
istream& operator>>(istream& s, gap_encoding& e)
{
    for(unsigned int i=0; i<e.m_genotype.size(); ++i)
    {
        s >> e.m_genotype[i];
    }
    return s;
}

/*!
 * \brief print the encoded data
 *
 * \author deong
 * \date 01/22/2009
 */
ostream& operator<<(ostream& s, const gsap_encoding& e)
{
    for(unsigned int i=0; i<e.m_genotype.size()-1; i++)
    {
        s << e.m_genotype[i] << " ";
    }
    s << e.m_genotype[e.m_genotype.size()-1];

    return s;
}

/*!
 * \brief read in an encoded value
 *
 * \author deong@acm.org
 * \date 01/22/2009
 */
istream& operator>>(istream& s, gsap_encoding& e)
{
    for(unsigned int i=0; i<e.m_genotype.size(); ++i)
    {
        s >> e.m_genotype[i];
    }
    return s;
}
