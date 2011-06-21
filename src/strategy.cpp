/*!
 * \file strategy.cpp
 *
 * how to apply a local search operator
 *
 * Deon Garrett
 * deong@acm.org
 */

#include <cstdlib>
#include "strategy.h"
#include "configuration.h"
#include "keywords.h"

using namespace std;

/*!
 * \brief return a local search strategy
 *
 * \author deong
 * \date 05/12/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
strategy strategy_factory::construct()
{
    string strname;
    configuration::string_parameter(keywords::HC_STRATEGY, strname, true);
    if(strname == keywords::STRATEGY_ALL)
    {
        return STRATEGY_ALL;
    }
    else if(strname == keywords::STRATEGY_RANDOM)
    {
        return STRATEGY_RANDOM;
    }
    else if(strname == keywords::STRATEGY_BEST)
    {
        return STRATEGY_BEST;
    }
    else if(strname == keywords::STRATEGY_WORST)
    {
        return STRATEGY_WORST;
    }
    else if(strname == keywords::STRATEGY_NONE)
    {
        return STRATEGY_NONE;
    }
    else
    {
        cerr << "invalid strategy specified: " << strname << endl;
        exit(1);
    }
}

    
