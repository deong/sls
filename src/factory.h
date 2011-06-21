/*!
 * \file factory.h
 *
 * Deon Garrett
 * University of Memphis
 * deong@acm.org
 */

#ifndef _FACTORY_H_
#define _FACTORY_H_

#include <string>

using std::string;

/*!
 * \class factory
 *
 * basic class for all factory classes; provides ability
 * to customize configuration variables by specifying a
 * prefix
 *
 * \author deong
 * \date 06/27/2007
 */
class factory
{
public:
    factory();
    virtual ~factory();
    void set_prefix(const string& prefix);

protected:
    string m_prefix;
};

#endif
