/*!
 * \file strategy.h
 *
 * how to apply a local search operator
 *
 * Deon Garrett
 * deong@acm.org
 */

#ifndef _STRATEGY_H_
#define _STRATEGY_H_

typedef enum
{
    STRATEGY_ALL,
    STRATEGY_RANDOM,
    STRATEGY_BEST,
    STRATEGY_WORST,
    STRATEGY_NONE
} strategy;

/*!
 * \class strategy_factory
 *
 * \author deong
 * \date 05/12/2007
 */
class strategy_factory
{
public:
    static strategy construct();
};

#endif


