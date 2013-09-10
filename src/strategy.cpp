/*!
 * \file strategy.cpp
 *
 * how to apply a local search operator
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <cstdlib>
#include "strategy.h"
#include "kvparse/kvparse.h"
#include "keywords.h"

using namespace std;

/*!
 * \brief return a local search strategy
 */
strategy strategy_factory::construct()
{
	string strname;
	kvparse::parameter_value(keywords::HC_STRATEGY, strname, true);
	if(strname == keywords::STRATEGY_ALL) {
		return STRATEGY_ALL;
	} else if(strname == keywords::STRATEGY_RANDOM) {
		return STRATEGY_RANDOM;
	} else if(strname == keywords::STRATEGY_BEST) {
		return STRATEGY_BEST;
	} else if(strname == keywords::STRATEGY_WORST) {
		return STRATEGY_WORST;
	} else if(strname == keywords::STRATEGY_NONE) {
		return STRATEGY_NONE;
	} else {
		cerr << "invalid strategy specified: " << strname << endl;
		exit(1);
	}
}


