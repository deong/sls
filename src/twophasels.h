/*!
 * \file twophasels.h
 *
 * two-phase local search
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _TWOPHASELS_H_
#define _TWOPHASELS_H_

#include "sls.h"
#include "localsearch.h"
#include "pfront.h"
#include "comparator.h"

using namespace std;

/*!
 * \class twophasels
 *
 * multiobjective random restarts local search
 */
template <template <typename> class Chromosome, typename Encoding>
class twophasels : virtual public sls<Chromosome,Encoding>
{
protected:
	unsigned int m_phase_two_runs;
	local_search<Chromosome,Encoding>* m_hc;
	local_search<Chromosome,Encoding>* m_hc2;
	scalarizing_comparator<Chromosome,Encoding> m_comp;
	pareto_front<Chromosome,Encoding> m_front;

private:
	// disable copying
	twophasels(const twophasels& that);
	twophasels& operator=(const twophasels& that);

public:
	twophasels();
	virtual ~twophasels();

	virtual void generation_completed();
	virtual void initialize();
	virtual void run();
};

#include "twophasels.cpp"

#endif
