/*!
 * \file utilities.h
 *
 * Defines some useful utility routines for performing tasks
 * such as printing warnings and errors.
 *
 * Deon Garrett
 * University of Memphis
 * deong@acm.org
 */

#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <iostream>
#include <string>
#include <vector>
#include <utility>

using namespace std;

void error(const string& msg);
void warning(const string& msg);
void debug(const string& msg);
int sort_compare(const double* item1, const double* item2);

double compute_pearson_correlation(const vector<pair<double,double> >& points);

template <typename T>
ostream& operator<<(ostream& ostr, const vector<T>& v)
{
    for(int i=0; i<int(v.size()); i++)
    {
        ostr << v[i];
        if(i!=int(v.size())-1)
        {
            ostr << " ";
        }
    }
    return ostr;
}

bool operator==(const vector<int>& v1, const vector<int>& v2);

#endif
