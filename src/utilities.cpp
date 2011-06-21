/*!
 * \file utilities.cpp
 *
 * Defines some useful utility routines for performing tasks
 * such as printing warnings and errors.
 *
 * Deon Garrett
 * deong@acm.org
 */

#define _DEBUG_

#include <iostream>
#include <cstdlib>
#include <string>
#include <utility>
#include <cmath>
#include "utilities.h"

using namespace std;

/*!
 * \brief print a warning message
 *
 * \author deong
 * \date 05/12/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
void warning(const string& msg)
{
    cerr << msg << endl;
}

/*!
 * \brief print an error message and exit
 *
 * \author deong
 * \date 05/12/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
void error(const string& msg)
{
    cerr << msg << endl;
    exit(1);
}

/*!
 * \brief print a debug message if debugging is enabled
 *
 * \author deong
 * \date 05/12/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
void debug(const string& msg)
{
#ifdef _DEBUG_
    cerr << msg << endl;
#endif
}

/*!
 * \brief compare items for sorting
 *
 * \author deong
 * \date 05/12/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
int sort_compare(const double* item1, const double* item2)
{
    if(*item1 < *item2)
    {
        return -1;
    }
    else if(*item1 > *item2)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

/*!
 * \brief compare two vectors for equality
 *
 * \author deong
 * \date 05/12/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
bool operator==(const vector<int>& v1, const vector<int>& v2)
{
    if(v1.size() != v2.size())
    {
        return false;
    }
    
    for(unsigned int i=0; i<v1.size(); i++)
    {
        if(v1[i] != v2[i])
        {
            return false;
        }
    }
    return true;
}

/*!
 * \brief compute correlation coefficient
 *
 * \author deong
 * \date 05/12/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY	DESCRIPTION
 * \endcode
 */
double compute_pearson_correlation(const vector<pair<double,double> >& points)
{
    double meanx = 0;
    double meany = 0;
    int num_points = static_cast<int>(points.size());

    for(int i=0; i<num_points; i++)
    {
        meanx += points[i].first;
        meany += points[i].second;
    }
    meanx /= num_points;
    meany /= num_points;

    double numerator = 0;
    double denom1 = 0;
    double denom2 = 0;

    for(int i=0; i<num_points; i++)
    {
        double xdev = points[i].first - meanx;
        double ydev = points[i].second - meany;
        numerator += xdev * ydev;
        denom1 += xdev * xdev;
        denom2 += ydev * ydev;
    }

    return numerator / (sqrt(denom1) * sqrt(denom2));
}
