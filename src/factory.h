/*!
 * \file factory.h
 *
 * Deon Garrett
 * jdgarrett@gmail.com
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
