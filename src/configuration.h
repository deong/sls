/*!
 * \file configuration.h
 *
 * Reads a configuration file into a parameter database.  To use it,
 * call the globalInstance method returning a pointer to the single
 * object in existence.
 *
 * Deon Garrett
 * deong@acm.org
 */

#ifndef _CONFIGURATION_H_
#define _CONFIGURATION_H_

#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include "utilities.h"

using namespace std;

/*!
 * \class configuration
 *
 * This class implements the singleton design
 * pattern.  It simply loads, stores, and provides access to
 * keyword-value pairs, which represent configuration parameters.
 *
 * \author deong
 * \date 05/08/2007
 */
class configuration
{
private:
    // my current collection of configuration parameters,
    // represented as keyword,value pairs.
    static map<string,list<string> > m_db;

    // disable copying for singleton classes
    configuration();
    configuration(const configuration&);
    configuration &operator=(const configuration&);

public:
    virtual ~configuration();

    static void clear();
    static bool read_configuration_file(const string &fileName);
    static int add_value(const string &keyword,const string &value);
    static int remove_value(const string &keyword,const string &value);
    static bool keyword_exists(const string &keyword);
    static bool has_unique_value(const string &keyword);
    static list<string> values(const string &keyword);
    static string value(const string &keyword);

    static bool string_parameter(const string &keyword, string& value, bool required = false);
    static bool integer_parameter(const string &keyword, int& value, bool required = false);
    static bool unsigned_integer_parameter(const string &keyword, unsigned int& value, bool required = false);
    static bool double_parameter(const string &keyword, double& value, bool required = false);
    static bool boolean_parameter(const string& keyword, bool& value, bool required = false);
    static bool list_parameter(const string& keyword, list<string>& value, bool required = false);

    /*!
     * \brief retrieve parameter values as a vector of a specified type
     *
     * \author deong
     * \date 05/08/2007
     *
     * \code
     * Modification History
     * MM/DD/YYYY	DESCRIPTION
     * \endcode
     */
    template <class T>
    static bool vector_parameter(const string& keyword, vector<T>& v, bool required = false)
    {
        if(!keyword_exists(keyword))
        {
            if(required)
            {
                error("required keyword: " + keyword + " not specified.");
            }
            else
            {
                return false;
            }
        }
        else
        {
            v.clear();
            string vecvals = value(keyword);
            istringstream istr(vecvals);
            while(!istr.eof())
            {
                T x;
                istr >> x;
                v.push_back(x);
            }
        }
        return true;
    }
    
    static void dump_contents(ostream &ostr);
};

#endif
