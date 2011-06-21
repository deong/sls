/*
 * enumerate.cc
 *
 * enumerate the space finding all pareto optimal solutions
 *
 * Deon Garrett
 * deong@acm.org
 */

#include <iostream>
#include <algorithm>
#include "chromosome.h"
#include "encoding.h"
#include "problems.h"
#include "pfront.h"

using namespace std;

void enumerate_bit_strings(const string&);
void enumerate_qap(string);
void enumerate_gap(string);

int main(int argc, char** argv)
{
    if(argc != 3)
    {
        cerr << "usage: enumerate <problem> <config>" << endl;
        exit(EXIT_FAILURE);
    }

    string probtype = argv[1];
    string conffile = argv[2];
    
    configuration::read_configuration_file(conffile);

    if(probtype == "qap")
    {
        enumerate_qap(argv[2]);
    }
    else if(probtype == "gap")
    {
        enumerate_gap(argv[2]);
    }
    else if(probtype=="knapsack")
    {
        enumerate_bit_strings(argv[2]);
    }
    else
    {
        cerr << "illegal problem type" << endl;
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}

void enumerate_bit_strings(const string& conffile)
{
    boolean_encoding::ProblemType* prob=boolean_encoding::ProblemFactoryType::construct();
    chromosome<boolean_encoding> p(prob);
    unsigned int maxval=static_cast<unsigned int>(pow(2.0,p.length())-1);
    for(unsigned int i=0; i<maxval; ++i)
    {
        unsigned int int_val=i;
        for(unsigned int b=0; b<p.length(); ++b)
        {
            p[b]=int_val%2;
            int_val/=2;
        }

        p.evaluate(prob);
        if(p.feasible())
            cout << p << endl;
    }
}

void enumerate_qap(string conffile)
{
    permutation_encoding::ProblemType* prob = permutation_encoding::ProblemFactoryType::construct();

    chromosome<permutation_encoding> p(prob);
    for(unsigned int i=0; i<p.length(); i++)
    {
        p[i] = i;
    }
    p.evaluate(prob);
    
//     pareto_front<chromosome,permutation_encoding> front;
//     front.add(p);
    
//     while(next_permutation(p.genotype().begin(), p.genotype().end()))
//     {
//         p.evaluate(prob);
//         front.add(p);
//     }

//     cout << front << endl;

    cout << p << endl;
    while(next_permutation(p.genotype().begin(),p.genotype().end()))
    {
        p.evaluate(prob);
        cout << p << endl;
    }
}

void enumerate_gap(string conffile)
{
    gap_encoding::ProblemType* prob = gap_encoding::ProblemFactoryType::construct();

    gap_encoding::initialize_parameters(prob);
    
    chromosome<gap_encoding> p(prob);
    pareto_front<chromosome,gap_encoding> front;

    unsigned int num_solutions = static_cast<unsigned int>(pow(static_cast<double>(p.agents()),
                                                               static_cast<double>(p.tasks())));

    for(unsigned int i=0; i<num_solutions; i++)
    {
        unsigned int val = i;
        for(unsigned int j=0; j<p.length(); j++)
        {
            p[j] = val % p.agents();
            val = val / p.agents();
        }

        p.evaluate(prob);
        if(p.feasible())
        {
            front.add(p);
        }
    }
    
    cout << front << endl;
}

    
