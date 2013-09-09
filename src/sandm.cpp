/**
 * \file sandm.cpp
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include "sandm.h"
#include "chromosome.h"
#include "mutation.h"
#include "pfront.h"
#include "configuration.h"
#include "keywords.h"
#include "mtrandom.h"

using namespace std;

/**
 * \brief constructor
 */
template <template <typename> class Chromosome, typename Encoding>
sandm<Chromosome,Encoding>::sandm() :
	sls<Chromosome,Encoding>(),
	m_sample_size(0),
	m_mut_op(0)
{
}

/**
 * \brief destructor
 */
template <template <typename> class Chromosome, typename Encoding>
sandm<Chromosome,Encoding>::~sandm()
{
	if(m_mut_op) {
		delete m_mut_op;
	}
}

/**
 * \brief initialize the algorithm components
 */
template <template <typename> class Chromosome, typename Encoding>
void sandm<Chromosome,Encoding>::initialize()
{
	// initialize the base sls algorithm
	sls<Chromosome,Encoding>::initialize();

	// initialize the mutation operator
	mutation_operator_factory<Chromosome,Encoding> mf;
	m_mut_op=mf.construct();
}

/**
 * \brief run the algorithm
 */
template <template <typename> class Chromosome, typename Encoding>
void sandm<Chromosome,Encoding>::run()
{
	//! generate a big initial sample
	for(unsigned int i=0; i<m_sample_size; i++) {
		Chromosome<Encoding> chr(this->m_fitfunc);
		chr.randomize();
		chr.evaluate(this->m_fitfunc);
		this->chromosome_evaluated(chr);
		cout << chr << endl;
		m_front.add(chr);
	}
	cout << "initial front: \n" << m_front << endl << endl;

	//! now take some fraction of the points on the front and mutate
	//! them, keeping anything that is nondominated
	mtrandom mt;
	while(!this->terminate()) {
		Chromosome<Encoding> child=m_front[mt.random(0,m_front.size())];
		m_mut_op->mutate(child);
		child.evaluate(this->m_fitfunc);
		this->chromosome_evaluated(child);
		m_front.add(child);
	}

	cout << "Final Pareto Front: \n" << m_front << endl << endl;
}
