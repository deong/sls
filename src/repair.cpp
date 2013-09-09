/**
 * \file repair.cpp
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <iostream>
#include <iterator>
#include "repair.h"
#include "chromosome.h"
#include "neighborhood.h"
#include "move.h"

using namespace std;

/**
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
repair_operator<Chromosome,Encoding>::repair_operator()
{
}

/**
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
repair_operator<Chromosome,Encoding>::~repair_operator()
{
}

/**
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
gap_repair_nd_impl<Chromosome,Encoding>::gap_repair_nd_impl()
{
    m_nf=new sss_neighborhood<Chromosome,Encoding>;
}

/**
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
gap_repair_nd_impl<Chromosome,Encoding>::~gap_repair_nd_impl()
{
    delete m_nf;
}

/**
 * \brief repair the individiual
 *
 * This uses a next-descent local search to find individuals that violate
 * capacity contraints by less and less until it either reaches a local
 * optima or finds an individual that does not violate constraints
 */
template <template <typename> class Chromosome, typename Encoding>
void gap_repair_nd_impl<Chromosome,Encoding>::repair(Chromosome<Encoding>& chr,
                                                     const typename Encoding::ProblemType* prob) const
{
    // first compute feasibility
    vector<int> used;
    used.assign(prob->agents,0);
    for(unsigned int i=0; i<chr.length(); i++)
    {
        used[chr[i]]+=prob->resources[chr[i]][i];
    }
    int excess=0;
    for(unsigned int i=0; i<used.size(); i++)
    {
        excess+=max(0,used[i]-prob->capacity[i]);
    }

    // already feasible
    if(excess==0)
        return;
    
    bool at_local_opt = false;
    while(!at_local_opt && excess>0)
    {
        at_local_opt = true;
        this->m_nf->initialize(chr);
        while(this->m_nf->has_more_neighbors())
        {
            move<Chromosome,Encoding> m=this->m_nf->next_neighbor();

            // analyze move to see if it would decrease excess resource usage
            vector<int> new_used=used;
            for(unsigned int i=0; i<m.length(); i++)
            {
                unsigned int task=m[i].first;
                int src_agt=chr[task];
                int dst_agt=m[i].second;
                new_used[src_agt]-=prob->resources[src_agt][task];
                new_used[dst_agt]+=prob->resources[dst_agt][task];
            }

            int new_excess=0;
            for(unsigned int i=0; i<new_used.size(); i++)
            {
                new_excess+=max(0,new_used[i]-prob->capacity[i]);
            }
            
            // check to see if the overall move decreased the excess and apply
            // any moves that move us closer to feasibility
            if(new_excess<excess)
            {
                for(unsigned int i=0; i<m.length(); i++)
                {
                    int task=m[i].first;
                    int src_agt=chr[task];
                    int dst_agt=m[i].second;
                    used[src_agt]-=prob->resources[src_agt][task];
                    used[dst_agt]+=prob->resources[dst_agt][task];
                    chr[task]=dst_agt;
                }
                excess=new_excess;
                at_local_opt=false;
                break;
            }
        }
    }
    // responsibility of client to evaluate new chromosome
    // chr.evaluate(prob);
}

/**
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
gap_repair_sd_impl<Chromosome,Encoding>::gap_repair_sd_impl()
{
    m_nf=new sss_neighborhood<Chromosome,Encoding>;
}

/**
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
gap_repair_sd_impl<Chromosome,Encoding>::~gap_repair_sd_impl()
{
    delete m_nf;
}

/**
 * \brief repair the individiual
 *
 * This uses a steepest-descent local search to find individuals that violate
 * capacity contraints by less and less until it either reaches a local
 * optima or finds an individual that does not violate constraints
 */
template <template <typename> class Chromosome, typename Encoding>
void gap_repair_sd_impl<Chromosome,Encoding>::repair(Chromosome<Encoding>& chr,
                                                     const typename Encoding::ProblemType* prob) const
{
    // first compute feasibility
    vector<int> used;
    used.assign(prob->agents,0);
    for(unsigned int i=0; i<chr.length(); i++)
    {
        used[chr[i]]+=prob->resources[chr[i]][i];
    }
    int excess=0;
    for(unsigned int i=0; i<used.size(); i++)
    {
        excess+=max(0,used[i]-prob->capacity[i]);
    }

    // already feasible
    if(excess==0)
        return;
    
    bool at_local_opt = false;
    while(!at_local_opt && excess>0)
    {
        move<Chromosome,Encoding> best_move;
        int best_move_excess=excess;

        at_local_opt = true;
        this->m_nf->initialize(chr);
        while(this->m_nf->has_more_neighbors())
        {
            move<Chromosome,Encoding> m=this->m_nf->next_neighbor();

            // analyze move to see if it would decrease excess resource usage
            vector<int> new_used=used;
            for(unsigned int i=0; i<m.length(); i++)
            {
                unsigned int task=m[i].first;
                int src_agt=chr[task];
                int dst_agt=m[i].second;
                new_used[src_agt]-=prob->resources[src_agt][task];
                new_used[dst_agt]+=prob->resources[dst_agt][task];
            }

            int new_excess=0;
            for(unsigned int i=0; i<new_used.size(); i++)
            {
                new_excess+=max(0,new_used[i]-prob->capacity[i]);
            }

            if(new_excess<best_move_excess)
            {
                best_move_excess=new_excess;
                best_move=m;
            }
        }
    
        // check to see if the overall move decreased the excess and apply
        // any moves that move us closer to feasibility
        if(best_move_excess<excess)
        {
            for(unsigned int i=0; i<best_move.length(); i++)
            {
                int task=best_move[i].first;
                int src_agt=chr[task];
                int dst_agt=best_move[i].second;
                used[src_agt]-=prob->resources[src_agt][task];
                used[dst_agt]+=prob->resources[dst_agt][task];
                chr[task]=dst_agt;
            }
            excess=best_move_excess;
            at_local_opt=false;
        }
    }
}

/**
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
gsap_repair_impl<Chromosome,Encoding>::gsap_repair_impl()
{
    m_nf=new sss_neighborhood<Chromosome,Encoding>;
}

/**
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
gsap_repair_impl<Chromosome,Encoding>::~gsap_repair_impl()
{
    delete m_nf;
}

/**
 * \brief repair the individiual
 */
template <template <typename> class Chromosome, typename Encoding>
void gsap_repair_impl<Chromosome,Encoding>::repair(Chromosome<Encoding>& chr,
                                                   const typename Encoding::ProblemType* prob) const
{
    mtrandom mt;

    const vector<gsap_problem::task_element>& elements=prob->get_elements();
    const vector<unsigned int>& tsb=prob->time_slot_boundaries();

    // set up the time slot boundaries
    for(unsigned int i=0; i<tsb.size()-1; ++i)
    {
        vector<int> indices(tsb[i+1]-tsb[i]);
        for(unsigned int j=tsb[i]; j<tsb[i+1]; ++j)
        {
            indices[j-tsb[i]]=j;
        }
        mt.shuffle(indices);

        // first, try to find agents for any unassigned tasks
        for(unsigned int j=0; j<indices.size(); ++j)
        {
            if(chr[indices[j]]==-1)
            {
                vector<int> others;
                prob->legal_values(indices[j],others);
                mt.shuffle(others);
                for(unsigned int k=0; k<others.size(); ++k)
                {
                    if(find(chr.begin()+tsb[i],chr.begin()+tsb[i+1],others[k])==chr.begin()+tsb[i+1])
                    {
                        chr[indices[j]]=others[k];
                        break;
                    }
                }
            }
        }
                
        // at this point, indices is a randomly shuffled vector of the positions
        // in the chromosome corresponding to time slot i

        // repair any violations caused by agents assigned to multiple tasks per
        // time slot
        vector<int> used;
        for(unsigned int j=0; j<indices.size(); ++j)
        {
            if(chr[indices[j]]==-1)
                continue;
                        
            if(find(used.begin(),used.end(),chr[indices[j]])!=used.end())
            {
                // chr[indices[i]] is repeated within the time slot, try to
                // repair the chromosome by shifting this task to another
                // agent
                vector<int> others;
                prob->legal_values(indices[j],others);
                mt.shuffle(others);

                bool fixed=false;
                for(unsigned int k=0; k<others.size(); ++k)
                {
                    // don't try and assign same sailor we already have
                    if(others[k]==chr[indices[j]])
                        continue;
                                        
                    // check if others[k] is already assigned in this time slot
                    // if not, make the switch
                    if(find(used.begin(),used.end(),others[k])==used.end())
                    {
                        chr[indices[j]]=others[k];
                        fixed=true;
                        break;
                    }
                }
                // if no change made, set chr[indices[j]]=-1
                if(!fixed)
                    chr[indices[j]]=-1;
            }
            if(chr[indices[j]]!=-1)
                used.push_back(chr[indices[j]]);
        }
    }
        
        
    // detect any capacity constraint violations
    vector<unsigned int> used;
    used.assign(prob->agents(),0);
    for(unsigned int i=0; i<chr.length(); i++)
    {
        if(chr[i]==-1)
            continue;
        used[chr[i]]+=prob->resources()[chr[i]][elements[i].task];
    }
    vector<int> excess;
    excess.assign(chr.length(),0);
    for(unsigned int i=0; i<used.size(); i++)
    {
        if(used[i]>prob->capacities()[i])
        {
            excess[i]+=used[i]-prob->capacities()[i];
        }
    }

    // now try to repair any capacity constraint violations
    vector<int> indices(chr.length());
    for(unsigned int i=0; i<chr.length(); ++i)
        indices[i]=i;
    mt.shuffle(indices);

    // go through every position in random order. if the agent at each
    // position is over capacity, try to reassign that task to an agent
    // that is under capacity
    for(unsigned int i=0; i<indices.size(); ++i)
    {
        if(chr[indices[i]]==-1 || excess[chr[indices[i]]]<=0)
            continue;

        // get a list of the other agents who can perform this task
        vector<int> others;
        prob->legal_values(indices[i],others);
        mt.shuffle(others);

        for(unsigned int j=0; j<others.size(); ++j)
        {
            // skip possible null assignment of current agent 
            if(others[j]==chr[indices[i]])
                continue;

            // what would happen to agent1 and agent2's resource usage if
            // I took this task from agent1 and gave it to agent2?

            // this condition makes the switch if it puts both agents under their limits
            if((used[others[j]]+prob->resources()[chr[indices[i]]][elements[indices[i]].task]<
                prob->capacities()[others[j]]) &&
               (used[chr[indices[i]]]-prob->resources()[chr[indices[i]]][elements[indices[i]].task]<
                prob->capacities()[chr[indices[i]]]))
            {
                chr[indices[i]]=others[j];
                break;
            }

            // this condition makes the switch if it decreases the excess of agent1
            // without making agent2 go over the limit, even if agent1 is still infeasible
            // if(used[others[j]]+prob->resources()[others[j]][elements[indices[i]].task]<
            //    prob->capacities()[others[j]])
            // {
            //     excess[chr[indices[i]]]-=prob->resources()[chr[indices[i]]][elements[indices[i]].task];
            //     excess[others[j]]+=prob->resources()[others[j]][elements[indices[i]].task];
            //     chr[indices[i]]=others[j];
            //     break;
            // }
        }
    }

    for(unsigned int i=0; i<excess.size(); ++i)
    {
        if(excess[i]>0)
        {
            break;
        }
    }
}

                        
/*!
 * \brief create an initialized repair operator for bit-vector encodings
 */
template <template <typename> class Chromosome, typename Encoding>
repair_operator<Chromosome,Encoding>* bit_vector_repair_factory<Chromosome,Encoding>::construct()
{
    return 0;
}

/*!
 * \brief create an initialized repair operator for real encodings
 */
template <template <typename> class Chromosome, typename Encoding>
repair_operator<Chromosome,Encoding>* real_repair_factory<Chromosome,Encoding>::construct()
{
    return 0;
}

/*!
 * \brief create an initialized repair operator for permutation encodings
 */
template <template <typename> class Chromosome, typename Encoding>
repair_operator<Chromosome,Encoding>* permutation_repair_factory<Chromosome,Encoding>::construct() 
{
    return 0;
}

/*!
 * \brief create an initialized repair operator for integer encodings
 */
template <template <typename> class Chromosome, typename Encoding>
repair_operator<Chromosome,Encoding>* integer_repair_factory<Chromosome,Encoding>::construct()
{
    return 0;
}

/*!
 * \brief create an initialized repair operator for gap encodings
 */
template <template <typename> class Chromosome, typename Encoding>
repair_operator<Chromosome,Encoding>* gap_repair_factory<Chromosome,Encoding>::construct()
{
    string repname;
    configuration::string_parameter(this->m_prefix+keywords::REPAIR_OPERATOR, repname, true);
    if(repname == keywords::GAP_REPAIR_ND)
    {
        repair_operator<Chromosome,Encoding>* rep=new gap_repair_nd<Chromosome,Encoding>;
        return rep;
    }
    else if(repname == keywords::GAP_REPAIR_SD)
    {
        repair_operator<Chromosome,Encoding>* rep=new gap_repair_sd<Chromosome,Encoding>;
        return rep;
    }
    else
    {
        cerr << "illegal value for parameter " << this->m_prefix+keywords::REPAIR_OPERATOR << ": "
             << repname << " specified" << endl;
        exit(0);
        return 0;
    }
}

/*!
 * \brief create an initialized repair operator for gsap encodings
 */
template <template <typename> class Chromosome, typename Encoding>
repair_operator<Chromosome,Encoding>* gsap_repair_factory<Chromosome,Encoding>::construct()
{
    string repname;
    configuration::string_parameter(this->m_prefix+keywords::REPAIR_OPERATOR, repname, true);
    if(repname == keywords::GSAP_REPAIR)
    {
        repair_operator<Chromosome,Encoding>* rep=new gsap_repair<Chromosome,Encoding>;
        return rep;
    }
    else
    {
        cerr << "illegal value for parameter " << this->m_prefix+keywords::REPAIR_OPERATOR << ": "
             << repname << " specified" << endl;
        exit(0);
        return 0;
    }
}

