/*!
 * \file factory.cpp
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include "factory.h"

/*!
 * \brief constructor
 */
factory::factory()
{
	m_prefix="";
}

/*!
 * \brief destructor
 */
factory::~factory()
{
}

/*!
 * \brief set the prefix string
 */
void factory::set_prefix(const string& prefix)
{
	m_prefix=prefix;
}
