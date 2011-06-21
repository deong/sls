/*!
 * \file chromosome.cpp
 *
 * Deon Garrett
 * deong@acm.org
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <iterator>
#include "chromosome.h"
#include "comparator.h"
#include "encoding.h"
#include "problems.h"
#include "configuration.h"
#include "keywords.h"

using namespace std;

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
basic_chromosome<Encoding>::basic_chromosome()
{
    m_parameters = new Encoding;
    m_feasible = true;
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
basic_chromosome<Encoding>::basic_chromosome(const typename Encoding::ProblemType* p)
{
    m_parameters = new Encoding(p);
    fitness.resize(p->objectives());
    m_feasible = true;
}

/*!
 * \brief copy constructor
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
basic_chromosome<Encoding>::basic_chromosome(const basic_chromosome<Encoding>& that)
{
    this->m_parameters = new Encoding;
    *(this->m_parameters) = *(that.m_parameters);
    this->fitness = that.fitness;
    this->m_feasible = that.m_feasible;
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
basic_chromosome<Encoding>::~basic_chromosome()
{
    delete m_parameters;
}

/*!
 * \brief assignment operator
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
basic_chromosome<Encoding>& basic_chromosome<Encoding>::operator=(const basic_chromosome& that)
{
    *(this->m_parameters) = *(that.m_parameters);
    this->fitness = that.fitness;
    this->m_feasible = that.m_feasible;
    return *this;
}

/*!
 * \brief inequality operator
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
bool basic_chromosome<Encoding>::operator!=(const basic_chromosome& that) const
{
    for(unsigned int i=0; i<this->m_parameters->length(); i++)
    {
        if((*m_parameters)[i] != (*(that.m_parameters))[i])
        {
            return true;
        }
    }
    return false;
}

/*!
 * \brief inequality operator specialized for real number encodings
 *
 * treats chromosomes as equal if no alleles differ by more than 1/10^8
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <>
bool basic_chromosome<real_encoding>::operator!=(const basic_chromosome& that) const
{
    for(unsigned int i=0; i<this->m_parameters->length(); i++)
    {
        if(fabs((*m_parameters)[i] - (*(that.m_parameters))[i]) > 1e-8)
        {
            return true;
        }
    }
    return false;
}

/*!
 * \brief equality comparison
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
bool basic_chromosome<Encoding>::operator==(const basic_chromosome& that) const
{
    for(unsigned int i=0; i<this->m_parameters->length(); i++)
    {
        if((*m_parameters)[i] != (*(that.m_parameters))[i])
        {
            return false;
        }
    }
    return true;
}

/*!
 * \brief less than comparison (needed to allow chromosomes to function as map keys)
 *
 * NOTE: don't use this -- it implements a mostly meaningless comparison
 * 
 * \author deong@acm.org
 * \date 10/15/2008
 */
template <typename Encoding>
bool basic_chromosome<Encoding>::operator<(const basic_chromosome& that) const
{
    for(unsigned int i=0; i<this->m_parameters->length(); ++i)
    {
        if((*m_parameters)[i] < (*(that.m_parameters))[i])
        {
            return true;
        }
        else if((*m_parameters)[i] > (*(that.m_parameters))[i])
        {
            return false;
        }
    }
    return false;
}

/*!
 * \brief return a given element from the encoded string
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
typename Encoding::Genotype& basic_chromosome<Encoding>::operator[](int index)
{
    return (*m_parameters)[index];
}

/*!
 * \brief return a const element from the encoded string
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
const typename Encoding::Genotype& basic_chromosome<Encoding>::operator[](int index) const
{
    return (*m_parameters)[index];
}

/*!
 * \brief iterator functionality
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
typename basic_chromosome<Encoding>::iterator basic_chromosome<Encoding>::begin()
{
    return this->m_parameters->begin();
}

/*!
 * \brief const_iterator functionality
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
typename basic_chromosome<Encoding>::const_iterator basic_chromosome<Encoding>::begin() const
{
    return this->m_parameters->begin();
}

/*!
 * \brief iterator functionality
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
typename basic_chromosome<Encoding>::iterator basic_chromosome<Encoding>::end()
{
    return this->m_parameters->end();
}

/*!
 * \brief const_iterator functionality
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
typename basic_chromosome<Encoding>::const_iterator basic_chromosome<Encoding>::end() const
{
    return this->m_parameters->end();
}

/*!
 * \brief compute the distance between the fitness vectors of two chromosomes
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
double basic_chromosome<Encoding>::fitness_distance(const basic_chromosome& that) const
{
    double d = 0.0;
    for(unsigned int i=0; i<fitness.size(); i++)
    {
        d += (static_cast<double>(fitness[i]) - static_cast<double>(that.fitness[i])) *
            (static_cast<double>(fitness[i]) - static_cast<double>(that.fitness[i]));
    }
    return sqrt(d);
}

/*!
 * \brief return the length of the encoded string
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
unsigned int basic_chromosome<Encoding>::length() const
{
    return m_parameters->length();
}

/*!
 * \brief return the encoding data for the chromosome
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
vector<typename Encoding::Genotype>& basic_chromosome<Encoding>::genotype()
{
    return m_parameters->genotype();
}

/*!
 * \brief return a const version of the chromosome's encoded data
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
const vector<typename Encoding::Genotype>& basic_chromosome<Encoding>::genotype() const
{
    return m_parameters->genotype();
}

/*!
 * \brief check chromosome for feasibility
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
bool basic_chromosome<Encoding>::feasible() const
{
    return m_feasible;
}

/*!
 * \brief randomize the chromosome
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
void basic_chromosome<Encoding>::randomize()
{
    m_parameters->randomize();
}

/*!
 * \brief evaluate the chromosome using the given fitness function
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
void basic_chromosome<Encoding>::evaluate(const typename Encoding::ProblemType* prob)
{
    this->m_parameters->decode();
    m_feasible = prob->evaluate(this->m_parameters->phenotype(), this->fitness);
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
chromosome<Encoding>::chromosome() :
    basic_chromosome<Encoding>::basic_chromosome()
{
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
chromosome<Encoding>::chromosome(const typename Encoding::ProblemType* p) :
    basic_chromosome<Encoding>::basic_chromosome(p)
{
}

/*!
 * \brief copy constructor
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
chromosome<Encoding>::chromosome(const chromosome& that) :
    basic_chromosome<Encoding>::basic_chromosome(that)
{
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
chromosome<Encoding>::~chromosome()
{
}

/*!
 * \brief assignment operator
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
chromosome<Encoding>& chromosome<Encoding>::operator=(const chromosome& that)
{
    basic_chromosome<Encoding>::operator=(that);
    return *this;
}

/*!
 * \brief equality comparison
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
bool chromosome<Encoding>::operator==(const chromosome& that) const
{
    return basic_chromosome<Encoding>::operator==(that);
}

/*!
 * \brief inequality comparison
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <typename Encoding>
bool chromosome<Encoding>::operator!=(const chromosome& that) const
{
    return basic_chromosome<Encoding>::operator!=(that);
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
chromosome<integer_encoding>::chromosome() :
    basic_chromosome<integer_encoding>::basic_chromosome()
{
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
chromosome<integer_encoding>::chromosome(const integer_encoding::ProblemType* p) :
    basic_chromosome<integer_encoding>::basic_chromosome(p)
{
}

/*!
 * \brief copy constructor
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
chromosome<integer_encoding>::chromosome(const chromosome& that) :
    basic_chromosome<integer_encoding>::basic_chromosome(that)
{
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
chromosome<integer_encoding>::~chromosome()
{
}

/*!
 * \brief assignment operator
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
chromosome<integer_encoding>& chromosome<integer_encoding>::operator=(const chromosome& that)
{
    basic_chromosome<integer_encoding>::operator=(that);
    return *this;
}

/*!
 * \brief equality comparison
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
bool chromosome<integer_encoding>::operator==(const chromosome& that) const
{
    return basic_chromosome<integer_encoding>::operator==(that);
}

/*!
 * \brief inequality comparison
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
bool chromosome<integer_encoding>::operator!=(const chromosome& that) const
{
    return basic_chromosome<integer_encoding>::operator!=(that);
}

/*!
 * \brief return a vector of allowed values for a given allele
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
const vector<int>& chromosome<integer_encoding>::legal_values(unsigned int index) const
{
    return this->m_parameters->legal_values(index);
}

/*!
 * \brief evaluate the fitness of the chromosome
 *
 * \author deong
 * \date 05/08/2007
 */
void chromosome<integer_encoding>::evaluate(const integer_encoding::ProblemType* p)
{
    this->m_parameters->decode();
    this->m_feasible = p->evaluate(this->m_parameters->phenotype(), this->fitness);
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/08/2007
 */
chromosome<gap_encoding>::chromosome()
{
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/08/2007
 */
chromosome<gap_encoding>::chromosome(const gap_encoding::ProblemType* p) :
    basic_chromosome<gap_encoding>::basic_chromosome(p)
{
    m_cap_used.assign(p->agents, 0);
}

/*!
 * \brief copy constructor
 *
 * \author deong
 * \date 05/08/2007
 */
chromosome<gap_encoding>::chromosome(const chromosome& that) :
    basic_chromosome<gap_encoding>::basic_chromosome(that)
{
    m_cap_used = that.m_cap_used;
    m_agt_map = that.m_agt_map;
    m_feasible = that.m_feasible;
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/08/2007
 */
chromosome<gap_encoding>::~chromosome()
{
}

/*!
 * \brief assignment operator
 *
 * \author deong
 * \date 05/08/2007
 */
chromosome<gap_encoding>& chromosome<gap_encoding>::operator=(const chromosome& that)
{
    basic_chromosome<gap_encoding>::operator=(that);
    m_cap_used = that.m_cap_used;
    m_agt_map = that.m_agt_map;
    m_feasible = that.m_feasible;
    return *this;
}

/*!
 * \brief equality comparator
 *
 * \author deong
 * \date 05/08/2007
 */
bool chromosome<gap_encoding>::operator==(const chromosome& that) const
{
    return basic_chromosome<gap_encoding>::operator==(that);
}

/*!
 * \brief inequality comparator
 *
 * \author deong
 * \date 05/08/2007
 */
bool chromosome<gap_encoding>::operator!=(const chromosome& that) const
{
    return basic_chromosome<gap_encoding>::operator!=(that);
}

/*!
 * \brief return the number of agents for the problem
 *
 * \author deong
 * \date 05/08/2007
 */
unsigned int chromosome<gap_encoding>::agents() const
{
    return this->m_parameters->agents();
}

/*!
 * \brief return the number of tasks for the problem
 *
 * \author deong
 * \date 05/08/2007
 */
unsigned int chromosome<gap_encoding>::tasks() const
{
    return this->m_parameters->tasks();
}

/*!
 * \brief compute the feasibility of the chromosome
 *
 * \author deong
 * \date 05/08/2007
 */
void chromosome<gap_encoding>::compute_infeasibility(const gap_encoding::ProblemType* p)
{
    // set up the agent mapping and set the capacity used to zero
    m_agt_map.clear();
    for(unsigned int i=0; i<p->agents; i++)
    {
        vector<int> v;
        m_agt_map.push_back(v);
        m_cap_used[i] = 0;
    }

    // compute the resources used by each agent
    this->m_feasible=true;
    for(unsigned int i=0; i<p->tasks; i++)
    {
        m_cap_used[(*(this->m_parameters))[i]] += p->resources[(*(this->m_parameters))[i]][i];
        if(m_cap_used[(*(this->m_parameters))[i]] > p->capacity[(*(this->m_parameters))[i]])
            this->m_feasible=false;
        m_agt_map[(*(this->m_parameters))[i]].push_back(i);
    }
}


/*!
 * \brief penalize an infeasible solution
 *
 * \author deong
 * \date 05/08/2007
 */
void chromosome<gap_encoding>::penalize(const gap_encoding::ProblemType* p)
{
    for(unsigned int k=0; k<this->fitness.size(); k++)
    {
        gap_encoding::FitnessType pi = 0;
        for(unsigned int i=0; i<this->m_parameters->length(); i++)
        {
            pi += static_cast<gap_encoding::FitnessType>
                (this->m_parameters->m_alpha[i] * max(0, m_cap_used[(*(this->m_parameters))[i]] -
                                                      p->capacity[(*(this->m_parameters))[i]]));
        }
        this->fitness[k] += pi;
    }
}

/*!
 * \brief attempt to repair an infeasible solution
 *
 * \author deong
 * \date 05/08/2007
 */
void chromosome<gap_encoding>::repair(const gap_encoding::ProblemType* p)
{
    // try to reach feasibility by shifting tasks to other sailors
    mtrandom mt;
    bool done=false;
    while(!done)
    {
        done=true;
        
        vector<int> order = mt.permutation(p->agents);
        for(unsigned int i=0; i<p->agents; i++)
        {
            int agt = order[i];
            if(m_cap_used[agt] > p->capacity[agt])
            {
                // too many resources used, shift job to another sailor
                vector<int> torder = mt.permutation(m_agt_map[agt].size());
                for(unsigned int j=0; j<m_agt_map[agt].size(); j++)
                {
                    int tsk = m_agt_map[agt][torder[j]];

                    // check if there is an agent that can take on the current task
                    vector<int> destorder = mt.permutation(p->agents);
                    bool change_made=false;
                    for(unsigned int k=0; k<p->agents; k++)
                    {
                        int dest = destorder[k];
                        if(dest==agt)
                        {
                            continue;
                        }
                                        
                        if(m_cap_used[dest] + p->resources[dest][tsk] <= p->capacity[dest])
                        {
                            (*(this->m_parameters))[tsk] = dest;
                            m_cap_used[agt] -= p->resources[agt][tsk];
                            m_cap_used[dest] += p->resources[dest][tsk];
                            m_agt_map[agt].erase(find(m_agt_map[agt].begin(), m_agt_map[agt].end(), tsk));
                            m_agt_map[dest].push_back(tsk);
                            change_made = true;
                            done=false;
                            break;
                        }
                        else
                        {
                            // try to swap with any task that is currently assigned to
                            // the destination agent
                            vector<int> morder = mt.permutation(m_agt_map[dest].size());
                            for(unsigned int m=0; m<m_agt_map[dest].size(); m++)
                            {
                                // OK...time for a quick glossary refresher
                                // agt = the original agent that was over capacity
                                // dest = the agent we're trying to swap with
                                // tsk = the task currently assigned to agt
                                // dsttsk = m_agt_map[dest][m] = the task currently assigned to dest
                                int dsttsk = m_agt_map[dest][morder[m]];
                                if((m_cap_used[agt] +
                                    p->resources[agt][dsttsk] -
                                    p->resources[agt][tsk] <= p->capacity[agt]) &&
                                   (m_cap_used[dest] +
                                    p->resources[dest][tsk] -
                                    p->resources[dest][dsttsk] <= p->capacity[dest]))
                                {
                                    swap((*(this->m_parameters))[tsk],(*(this->m_parameters))[dsttsk]);
                                    m_cap_used[agt] += (p->resources[agt][dsttsk] - p->resources[agt][tsk]);
                                    m_cap_used[dest] += (p->resources[dest][tsk] - p->resources[dest][dsttsk]);
                                    m_agt_map[agt].erase(find(m_agt_map[agt].begin(),m_agt_map[agt].end(),tsk));
                                    m_agt_map[agt].push_back(dsttsk);
                                    m_agt_map[dest].erase(find(m_agt_map[dest].begin(),m_agt_map[dest].end(),dsttsk));
                                    m_agt_map[dest].push_back(tsk);
                                    change_made = true;
                                    done=false;
                                    break;
                                }
                            }
                            // if swap made, do another break to reset dest agent search
                            if(change_made)
                            {
                                break;
                            }
                        }
                    }

                    // change_made breaks put us here; in the middle of searching for
                    // a task to shuffle about.  If we've already found a change to
                    // make, we can break out of this loop as well
                    if(change_made)
                    {
                        cout << "change made...\n";
                        cout << "genotype: " << *(this->m_parameters) << endl;
                        cout << "max_capacities: ";
                        for(unsigned int i=0; i<p->agents; i++)
                        {
                            cout << setw(4) << p->capacity[i];
                        }
                        cout << "\nresources_used: ";
                        for(unsigned int i=0; i<p->agents; i++)
                        {
                            cout << setw(4) << m_cap_used[i];
                        }
                        cout << endl;
                        break;
                    }
                }
                // if this agent is within capacity constraints now, no need to keep
                // mucking about with him
                if(m_cap_used[agt] <= p->capacity[agt])
                {
                    break;
                }
            }
        }
    }// while(!done)
    
    this->m_feasible = true;
    for(unsigned int i=0; i<p->agents; i++)
    {
        if(m_cap_used[i] > p->capacity[i])
        {
            this->m_feasible = false;
            break;
        }
    }
    exit(0);
}

/*!
 * \brief evaluate the fitness of the chromosome
 *
 * \author deong
 * \date 05/08/2007
 */
void chromosome<gap_encoding>::evaluate(const gap_encoding::ProblemType* p)
{
    this->m_parameters->decode();
    this->compute_infeasibility(p);
    p->evaluate(this->m_parameters->phenotype(), this->fitness);
    if(!this->m_feasible)
    {
        penalize(p);
    }
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 01/22/2009
 */
chromosome<gsap_encoding>::chromosome()
{
}

/*!
 * \brief constructor
 *
 * \author deong
 * \date 01/22/2009
 */
chromosome<gsap_encoding>::chromosome(const gsap_encoding::ProblemType* p) :
    basic_chromosome<gsap_encoding>::basic_chromosome(p)
{
    m_cap_used.assign(p->agents(), 0);
}

/*!
 * \brief copy constructor
 *
 * \author deong
 * \date 01/22/2009
 */
chromosome<gsap_encoding>::chromosome(const chromosome& that) :
    basic_chromosome<gsap_encoding>::basic_chromosome(that)
{
    m_cap_used = that.m_cap_used;
    m_feasible = that.m_feasible;
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 01/22/2009
 */
chromosome<gsap_encoding>::~chromosome()
{
}

/*!
 * \brief assignment operator
 *
 * \author deong
 * \date 01/22/2009
 */
chromosome<gsap_encoding>& chromosome<gsap_encoding>::operator=(const chromosome& that)
{
    basic_chromosome<gsap_encoding>::operator=(that);
    m_cap_used = that.m_cap_used;
    m_feasible = that.m_feasible;
    return *this;
}

/*!
 * \brief equality comparator
 *
 * \author deong
 * \date 01/22/2009
 */
bool chromosome<gsap_encoding>::operator==(const chromosome& that) const
{
    return basic_chromosome<gsap_encoding>::operator==(that);
}

/*!
 * \brief inequality comparator
 *
 * \author deong
 * \date 01/22/2009
 */
bool chromosome<gsap_encoding>::operator!=(const chromosome& that) const
{
    return basic_chromosome<gsap_encoding>::operator!=(that);
}

/*!
 * \brief return the number of agents for the problem
 *
 * \author deong
 * \date 01/22/2009
 */
unsigned int chromosome<gsap_encoding>::agents() const
{
    return 0;
}

/*!
 * \brief return the number of tasks for the problem
 *
 * \author deong
 * \date 01/22/2009
 */
unsigned int chromosome<gsap_encoding>::tasks() const
{
        return 0;
}

/*!
 * \brief compute the feasibility of the chromosome
 *
 * \author deong
 * \date 01/22/2009
 */
void chromosome<gsap_encoding>::compute_infeasibility(const gsap_encoding::ProblemType* p)
{
        this->m_feasible=true;
        m_excess=0;
        n_unassigned=0;
        m_cap_used.assign(p->agents(),0);
        for(unsigned int i=0; i<m_parameters->length(); ++i)
        {
                if((*(this->m_parameters))[i]==-1)
                {
                        ++n_unassigned;
                        continue;
                }
                
        m_cap_used[(*(this->m_parameters))[i]] += p->resources()[(*(this->m_parameters))[i]]
                        [p->get_elements()[i].task];
        if(m_cap_used[(*(this->m_parameters))[i]] > p->capacities()[(*(this->m_parameters))[i]])
                {
                        m_excess+=m_cap_used[(*(this->m_parameters))[i]] > p->capacities()[(*(this->m_parameters))[i]];
            this->m_feasible=false;
                }
        }
        if(n_unassigned>0)
                this->m_feasible=false;
//      cout << "n_unassigned = " << n_unassigned << endl;
//      cout << "m_excess = " << m_excess << endl;
//      cout << "m_cap_used = ";
//      copy(m_cap_used.begin(),m_cap_used.end(),ostream_iterator<unsigned int>(cout," "));
//      cout << endl;
}

/*!
 * \brief penalize an infeasible solution
 *
 * \author deong
 * \date 01/22/2009
 */
void chromosome<gsap_encoding>::penalize(const gsap_encoding::ProblemType* p)
{
    for(unsigned int i=0; i<this->fitness.size(); ++i)
    {
        this->fitness[i] += (n_unassigned*this->m_parameters->m_unass_pen);
        this->fitness[i] += (m_excess*this->m_parameters->m_cap_pen);
    }
}

/*!
 * \brief attempt to repair an infeasible solution
 *
 * \author deong
 * \date 01/22/2009
 */
void chromosome<gsap_encoding>::repair(const gsap_encoding::ProblemType* p)
{
}

/*!
 * \brief evaluate the fitness of the chromosome
 *
 * \author deong
 * \date 01/22/2009
 */
void chromosome<gsap_encoding>::evaluate(const gsap_encoding::ProblemType* p)
{
    this->m_parameters->decode();
    this->compute_infeasibility(p);
        p->evaluate(this->m_parameters->phenotype(), this->fitness);
    if(!this->m_feasible)
    {
        penalize(p);
    }
}

/*!
 * \brief print a human readable representation of the chromosome
 *
 * \author deong
 * \date 05/08/2007
 */
template <typename Encoding>
ostream& operator<<(ostream& ostr, const chromosome<Encoding>& sol)
{
    ostr << "genotype: " << *(sol.m_parameters) << endl;
    ostr << "phenotype:";
    for(unsigned int i=0; i<sol.m_parameters->phenotype().size(); i++)
    {
        ostr << " " << sol.m_parameters->phenotype()[i];
    }
    ostr << endl;
    ostr << "fitness:";
    for(unsigned int i=0; i<sol.fitness.size(); i++)
    {
        ostr << " " << setprecision(11) << sol.fitness[i];
    }
    return ostr;
}

/*!
 * \brief read in a chromosome from a textual representation
 *
 * \author deong@acm.org
 * \date 10/14/2008
 */
template <typename Encoding>
istream& operator>>(istream& istr, chromosome<Encoding>& sol)
{
    string junk;
    istr >> junk;
    if(junk!="genotype:")
    {
        cout << "error reading chromosome: missing genotype token" << endl;
        return istr;
    }
    istr >> *(sol.m_parameters);
    istr >> junk;
    if(junk!="phenotype:")
    {
        cout << "error reading chromosome: missing phenotype token" << endl;
        return istr;
    }
    for(unsigned int i=0; i<sol.m_parameters->phenotype().size(); i++)
    {
        istr >> sol.m_parameters->phenotype()[i];
    }
    istr >> junk;
    if(junk!="fitness:")
    {
        cout << "error reading chromosome: missing fitness token" << endl;
        return istr;
    }
    for(unsigned int i=0; i<sol.fitness.size(); i++)
    {
        istr >> sol.fitness[i];
    }
    return istr;
}
    
/*!
 * \brief print a human readable representation of the chromosome
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <>
ostream& operator<<(ostream& ostr, const chromosome<boolean_encoding>& sol)
{
    ostr << "genotype: " << *(sol.m_parameters) << endl;
    ostr << "fitness:";
    for(unsigned int i=0; i<sol.fitness.size(); i++)
    {
        ostr << " " << setprecision(11) << sol.fitness[i];
    }
    return ostr;
}

/*!
 * \brief read from a textual representation
 *
 * \author deong@acm.org
 * \date 10/14/2008
 */
template <>
istream& operator>>(istream& istr, chromosome<boolean_encoding>& sol)
{
    string junk;
    istr >> junk;
    if(junk!="genotype:")
    {
        cout << "error reading chromosome: missing genotype token (" << junk << ")" << endl;
        return istr;
    }
    istr >> *(sol.m_parameters);
    istr >> junk;
    if(junk!="fitness:")
    {
        cout << "error reading chromosome: missing fitness token (" << junk << ")" << endl;
        return istr;
    }
    for(unsigned int i=0; i<sol.fitness.size(); i++)
    {
        istr >> sol.fitness[i];
    }
    return istr;
}

/*!
 * \brief print a human readable representation of the chromosome
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <>
ostream& operator<<(ostream& ostr, const chromosome<permutation_encoding>& sol)
{
    ostr << "genotype: " << *(sol.m_parameters) << endl;
    ostr << "fitness:";
    for(unsigned int i=0; i<sol.fitness.size(); i++)
    {
        ostr << " " << setprecision(11) << sol.fitness[i];
    }
    return ostr;
}

/*!
 * \brief read from a textual representation
 *
 * \author deong@acm.org
 * \date 10/14/2008
 */
template <>
istream& operator>>(istream& istr, chromosome<permutation_encoding>& sol)
{
    string junk;
    istr >> junk;
    if(junk!="genotype:")
    {
        cout << "error reading chromosome: missing genotype token" << endl;
        return istr;
    }
    istr >> *(sol.m_parameters);
    istr >> junk;
    if(junk!="fitness:")
    {
        cout << "error reading chromosome: missing fitness token" << endl;
        return istr;
    }
    for(unsigned int i=0; i<sol.fitness.size(); i++)
    {
        istr >> sol.fitness[i];
    }
    return istr;
}

/*!
 * \brief print a human readable representation of the chromosome
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <>
ostream& operator<<(ostream& ostr, const chromosome<real_encoding>& sol)
{
    ostr << "genotype: " << *(sol.m_parameters) << endl;
    ostr << "fitness:";
    for(unsigned int i=0; i<sol.fitness.size(); i++)
    {
        ostr << " " << setprecision(11) << sol.fitness[i];
    }
    return ostr;
}

/*!
 * \brief read from a textual representation
 *
 * \author deong@acm.org
 * \date 10/14/2008
 */
template <>
istream& operator>>(istream& istr, chromosome<real_encoding>& sol)
{
    string junk;
    istr >> junk;
    if(junk!="genotype:")
    {
        cout << "error reading chromosome: missing genotype token" << endl;
        return istr;
    }
    istr >> *(sol.m_parameters);
    istr >> junk;
    if(junk!="fitness:")
    {
        cout << "error reading chromosome: missing fitness token" << endl;
        return istr;
    }
    for(unsigned int i=0; i<sol.fitness.size(); i++)
    {
        istr >> sol.fitness[i];
    }
    return istr;
}
 
/*!
 * \brief print a human readable representation of the chromosome
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <>
ostream& operator<<(ostream& ostr, const chromosome<integer_encoding>& sol)
{
    ostr << "genotype: " << *(sol.m_parameters) << endl;
    ostr << "fitness:";
    for(unsigned int i=0; i<sol.fitness.size(); i++)
    {
        ostr << " " << setprecision(11) << sol.fitness[i];
    }
    return ostr;
}

/*!
 * \brief read from a textual representation
 *
 * \author deong@acm.org
 * \date 10/14/2008
 */
template <>
istream& operator>>(istream& istr, chromosome<integer_encoding>& sol)
{
    string junk;
    istr >> junk;
    if(junk!="genotype:")
    {
        cout << "error reading chromosome: missing genotype token" << endl;
        return istr;
    }
    istr >> *(sol.m_parameters);
    istr >> junk;
    if(junk!="fitness:")
    {
        cout << "error reading chromosome: missing fitness token" << endl;
        return istr;
    }
    for(unsigned int i=0; i<sol.fitness.size(); i++)
    {
        istr >> sol.fitness[i];
    }
    return istr;
}

/*!
 * \brief print a human readable representation of the chromosome
 *
 * \author deong
 * \date 05/08/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <>
ostream& operator<<(ostream& ostr, const chromosome<gap_encoding>& sol)
{
    ostr << "genotype: " << *(sol.m_parameters) << endl;
    ostr << "fitness:";
    for(unsigned int i=0; i<sol.fitness.size(); i++)
    {
        ostr << " " << setprecision(11) << sol.fitness[i];
    }
    ostr << endl;
    ostr << "feasible: ";
    if(sol.m_feasible)
    {
        ostr << "true";
    }
    else
    {
        ostr << "false" << endl;
        ostr << "resource_usage: " << sol.m_cap_used;
    }
    ostr << endl;
    return ostr;
}

/*!
 * \brief read from a textual representation
 *p
 * \author deong@acm.org
 * \date 10/14/2008
 */
template <>
istream& operator>>(istream& istr, chromosome<gap_encoding>& sol)
{
    string junk;
    istr >> junk;
    if(junk!="genotype:")
    {
        cout << "error reading chromosome: missing genotype token" << endl;
        return istr;
    }
    istr >> *(sol.m_parameters);
    istr >> junk;
    if(junk!="fitness:")
    {
        cout << "error reading chromosome: missing fitness token" << endl;
        return istr;
    }
    for(unsigned int i=0; i<sol.fitness.size(); i++)
    {
        istr >> sol.fitness[i];
    }
    istr >> junk;
    if(junk!="feasible:")
    {
        cout << "error reading chromosome: missing feasible token" << endl;
        return istr;
    }
    istr >> junk;
    if(junk=="true")
        sol.m_feasible=true;
    else
    {
        sol.m_feasible=false;
        istr >> junk;
        if(junk!="resource_usage:")
        {
            cout << "error reading chromosome: missing resource_usage token" << endl;
            return istr;
        }
        for(unsigned int i=0; i<sol.m_cap_used.size(); i++)
        {
            istr >> sol.m_cap_used[i];
        }
    }
    return istr;
}

/*!
 * \brief print a human readable representation of the chromosome
 *
 * \author deong
 * \date 01/22/2009
 */
template <>
ostream& operator<<(ostream& ostr, const chromosome<gsap_encoding>& sol)
{
    ostr << "genotype: " << *(sol.m_parameters) << endl;
    ostr << "fitness:";
    for(unsigned int i=0; i<sol.fitness.size(); i++)
    {
        ostr << " " << setprecision(11) << sol.fitness[i];
    }
    ostr << endl;
    ostr << "feasible: ";
    if(sol.m_feasible)
    {
        ostr << "true";
    }
    else
    {
        ostr << "false" << endl;
        ostr << "resource_usage: " << sol.m_cap_used;
    }
    ostr << endl;
    return ostr;
}

/*!
 * \brief read from a textual representation
 *
 * \author deong@acm.org
 * \date 01/22/2009
 */
template <>
istream& operator>>(istream& istr, chromosome<gsap_encoding>& sol)
{
    string junk;
    istr >> junk;
    if(junk!="genotype:")
    {
        cout << "error reading chromosome: missing genotype token" << endl;
        return istr;
    }
    istr >> *(sol.m_parameters);
    istr >> junk;
    if(junk!="fitness:")
    {
        cout << "error reading chromosome: missing fitness token" << endl;
        return istr;
    }
    for(unsigned int i=0; i<sol.fitness.size(); i++)
    {
        istr >> sol.fitness[i];
    }
    istr >> junk;
    if(junk!="feasible:")
    {
        cout << "error reading chromosome: missing feasible token" << endl;
        return istr;
    }
    istr >> junk;
    if(junk=="true")
        sol.m_feasible=true;
    else
    {
        sol.m_feasible=false;
        istr >> junk;
        if(junk!="resource_usage:")
        {
                        cout << "error reading chromosome: missing resource_usage token" << endl;
            return istr;
        }
        for(unsigned int i=0; i<sol.m_cap_used.size(); i++)
        {
            istr >> sol.m_cap_used[i];
        }
    }
    return istr;
}

