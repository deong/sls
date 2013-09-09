/*!
 * \file pfront.cpp
 *
 * class representing a pareto front
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include "pfront.h"
#include "chromosome.h"
#include "population.h"
#include "comparator.h"
#include "problems.h"

using namespace std;

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
pareto_front<Chromosome,Encoding>::pareto_front() :
	population<Chromosome,Encoding>::population()
{
}

/*!
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
pareto_front<Chromosome,Encoding>::pareto_front(unsigned int sz) :
	population<Chromosome,Encoding>::population(sz)
{
}

/*!
 * \brief copy constructor
 */
template <template <typename> class Chromosome, typename Encoding>
pareto_front<Chromosome,Encoding>::pareto_front(const pareto_front& pf)
{
	this->m_individuals = pf.m_individuals;
	this->m_count = pf.m_count;
}

/*!
 * \brief population constructor
 *
 * fills front with all nondominated individuals from the population
 */
template <template <typename> class Chromosome, typename Encoding>
pareto_front<Chromosome,Encoding>::pareto_front(const population<Chromosome,Encoding>& pop)
{
	for(unsigned int i=0; i<pop.size(); i++) {
		add(pop[i]);
	}
}

/*!
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
pareto_front<Chromosome,Encoding>::~pareto_front()
{
}

/*!
 * \brief add a new chromosome to the front
 *
 * verifies that the chromosome is nondominated and removes any members
 * from the front that become dominated as a result of the new chromosome
 */
template <template <typename> class Chromosome, typename Encoding>
bool pareto_front<Chromosome,Encoding>::add(const Chromosome<Encoding>& chr)
{
	weak_dominance_comparator<Chromosome,Encoding> pdc;

	bool addi = true;
	for(int j=int(this->m_count-1); j>=0; j--) {
		int compres = pdc.compare(chr, this->m_individuals[j]);
		if(compres == 1) {
			addi = false;
			break;
		} else if(compres == -1) {
			this->remove_at(j);
		}
	}
	if(addi) {
		population<Chromosome,Encoding>::add(chr);
	}

	return addi;
}

/*!
 * \brief construct a pareto front from a given population
 */
template <template <typename> class Chromosome, typename Encoding>
void pareto_front<Chromosome,Encoding>::construct_front(const population<Chromosome,Encoding>& pop)
{
	this->clear();
	for(unsigned int i=0; i<pop.size(); i++) {
		this->add(pop[i]);
	}
}

/*!
 * \brief read in a pareto front from a given file
 */
template <template <typename> class Chromosome, typename Encoding>
void pareto_front<Chromosome,Encoding>::construct_front(string filename, const typename Encoding::ProblemType* prob)
{
	ifstream in(filename.c_str());
	if(!in) {
		cerr << "error opening file: " << filename << endl;
		exit(1);
	}

	while(!in.eof()) {
		char line[1024];
		in.getline(line,1024);
		if(strncmp(line,"genotype",8)==0) {
			// read the genotype
			vector<typename Encoding::Genotype> gene;
			string strline(line);
			istringstream istr(strline);
			string tok;
			istr >> tok;        // eat the initial token (genotype:)
			while(istr >> tok) {
				// read as a float and cast to more specific type
				typename Encoding::Genotype itok = static_cast<typename Encoding::Genotype>(atof(tok.c_str()));
				gene.push_back(itok);
			}

			// read the fitness
			vector<typename Encoding::FitnessType> fit;
			in.getline(line,1024);
			if(strncmp(line,"fitness",7)!=0) {
				cerr << "error in data file: genotype not followed by fitness" << endl;
				exit(1);
			}
			istringstream istr2(line);
			istr2 >> tok;       // eat the initial token (fitness:)
			while(istr2 >> tok) {
				int itok = atoi(tok.c_str());
				fit.push_back(itok);
			}

			// create the chromosome and insert it into the front
			Chromosome<Encoding> chr(prob);
			chr.fitness = fit;
			for(unsigned int i=0; i<chr.length(); i++) {
				chr[i] = gene[i];
			}
			this->add(chr);
		}
	}
	in.close();
}

/*!
 * \brief display the solutions in the pareto front
 */
template <template <typename> class Chromosome, typename Encoding>
ostream& operator<<(ostream& s, const pareto_front<Chromosome,Encoding>& pf)
{
	s << "Pareto Front" << endl;
	for(unsigned int i=0; i<pf.size(); i++) {
		//s << setw(11) << pf[i].fitness << endl;
		s << pf[i] << endl;
	}
	return s;
}

