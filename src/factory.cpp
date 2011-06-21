/*!
 * \file factory.cpp
 *
 * Deon Garrett
 * University of Memphis
 * deong@acm.org
 */

#include "factory.h"

/*!
 * \brief constructor
 *
 * \author deong
 * \date 06/27/2007
 */
factory::factory()
{
    m_prefix="";
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 06/27/2007
 */
factory::~factory()
{
}

/*!
 * \brief set the prefix string
 *
 * \author deong
 * \date 06/27/2007
 */
void factory::set_prefix(const string& prefix)
{
    m_prefix=prefix;
}
