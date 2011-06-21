/*!
 * \file gapgen.cpp
 *
 * problem generator for the multiobjective generalized assignment problem
 *
 * Deon Garrett
 * deong@acm.org
 *
 *
 * parameters
 * ----------
 * --agents (-n): number of agents
 * --tasks  (-m): number of tasks
 * --obj    (-k): number of objectives
 * --type   (-t): type of instance {A, C, D, or E}
 * --frac   (-f): fraction of correlated costs
 * --coef   (-r): correlation coefficient for costs
 * --slack  (-s): slack in resource constraints
 *
 * output format
 * -------------
 * agents tasks objectives
 *
 * x_1 x_2 x_3 ... x_n  (capacity constraints for each agent)
 *
 * r_11 r_12 r_13 ... r_1m  (resources required for agent 1 to do each task)
 * r_21 r_22 r_23 ... r_2m  (resources required for agent 2 to do each task)
 * ...
 * r_n1 r_n2 r_n3 ... r_nm  (resources required for agent n to do each task)
 *
 * c_111 c_112 c_113 ... c_11m (cost 1 for agent 1 to do each task)
 * c_121 c_122 c_123 ... c_12m (cost 1 for agent 2 to do each task)
 * ...
 * c_1n1 c_1n2 c_1n3 ... c_1nm (cost 1 for agent n to do each task)
 *
 * c_211 c_212 c_213 ... c_21m (cost 2 for agent 1 to do each task)
 * c_221 c_222 c_223 ... c_22m (cost 2 for agent 2 to do each task)
 * ...
 * c_2n1 c_2n2 c_2n3 ... c_2nm (cost 2 for agent n to do each task)
 *
 * ...
 * 
 * c_k11 c_k12 c_k13 ... c_k1m (cost k for agent 1 to do each task)
 * c_k21 c_k22 c_k23 ... c_k2m (cost k for agent 2 to do each task)
 * ...
 * c_kn1 c_kn2 c_kn3 ... c_knm (cost k for agent n to do each task)
 *
 */

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <ctime>
#include <unistd.h>
#include "mtrandom.h"

using namespace std;

void print_usage();
void generate_type_a(unsigned int n, unsigned int m, unsigned int k, double slack, double r, double f);
void generate_type_c(unsigned int n, unsigned int m, unsigned int k, double slack, double r, double f);
void generate_type_d(unsigned int n, unsigned int m, unsigned int k, double slack, double r, double f);
void generate_type_e(unsigned int n, unsigned int m, unsigned int k, double slack, double r, double f);
double correl_val(double x, double r);

int main(int argc, char** argv)
{
    if(argc == 1)
    {
        print_usage();
        exit(1);
    }

    unsigned int n = 0;
    unsigned int m = 0;
    unsigned int k = 0;
    char type = 'C';
    bool n_spec = false;
    bool m_spec = false;
    bool k_spec = false;
    bool type_spec = false;
    double       f = 0.0;
    double       r = 0.0;
    double   slack = 1.0;
    
    //mtrandom::initialize();
    mtrandom rng;
    
    char curropt;
    while((curropt = getopt(argc, argv, "hn:m:k:t:s:f::r::")) != -1)
    {
        switch(curropt)
        {
        case 'n':
            n = static_cast<unsigned int>(atoi(optarg));
            n_spec = true;
            break;
        case 'm':
            m = static_cast<unsigned int>(atoi(optarg));
            m_spec = true;
            break;
        case 'k':
            k = static_cast<unsigned int>(atoi(optarg));
            k_spec = true;
            break;
        case 't':
            type = optarg[0];
            type_spec = true;
            break;
	case 's':
	    slack=atof(optarg);
	    break;
	case 'f':
            f = atof(optarg);
            break;
        case 'r':
            r = atof(optarg);
            break;
        case 'h':
            print_usage();
            break;
        default:
            cout << "invalid option specified" << endl;
            exit(1);
        }
    }

    if(!n_spec || !m_spec || !k_spec || !type_spec)
    {
        print_usage();
        exit(1);
    }
    
    cout << ";; type: " << type << "\t";
    cout << "agents: " << n << "\t";
    cout << "tasks: " << m << "\t";
    cout << "objectives: " << k << "\t";
    cout << "slack: " << slack << "\t";
    cout << "corrfrac: " << f << "\t";
    cout << "corrcoef: " << r << endl << endl;

    cout << n << " " << m << " " << k << endl << endl;
    
    switch(type)
    {
    case 'A':
	generate_type_a(n,m,k,slack,r,f);
	break;
    case 'C':
        generate_type_c(n,m,k,slack,r,f);
        break;
    case 'D':
        generate_type_d(n,m,k,slack,r,f);
        break;
    case 'E':
        generate_type_e(n,m,k,slack,r,f);
        break;
    default:
        cout << "invalid instance type specified." << endl;
        exit(1);
    }

    return 0;
}

void print_usage()
{
    cout << "gapgen: generates instances of multiobjective generalized assignment problem" << endl << endl;
    cout << "parameters: [default values, * if required]" << endl;
    cout << "\t-h: print this help screen" << endl;
    cout << "\t-n <agents> [*]" << endl;
    cout << "\t-m <tasks> [*]" << endl;
    cout << "\t-k <objectives> [*]" << endl;
    cout << "\t-t <type (C, D, or E)> [*]" << endl;
    cout << "\t-s <slack>" << endl;
    cout << "\t-f <frac> fraction of costs correlated [0.0]" << endl;
    cout << "\t-r <coef> correlation coefficient [0.0]" << endl << endl;
    exit(1);
}

void generate_type_c(unsigned int n, unsigned int m, unsigned int k, double slack, double r, double f)
{
    mtrandom rng;
    
    // generate the resource allocations
    unsigned int** res = new unsigned int*[n];
    for(unsigned int i=0; i<n; i++)
    {
        res[i] = new unsigned int[m];
        for(unsigned int j=0; j<m; j++)
        {
            res[i][j] = rng.random(5, 25);
        }
    }

    // generate the cost matrices
    unsigned int*** costs = new unsigned int**[k];
    for(unsigned int obj=0; obj<k; obj++)
    {
        costs[obj] = new unsigned int*[n];
        for(unsigned int i=0; i<n; i++)
        {
            costs[obj][i] = new unsigned int[m];
        }
    }
    
    for(unsigned int i=0; i<n; i++)
    {
        for(unsigned int j=0; j<m; j++)
        {
            for(unsigned int obj=0; obj<k; obj++)
            {
                if(obj==0)
                {
                    costs[obj][i][j] = rng.random(10, 50);
                }
                else if(rng.random() < f) 
                {
                    double x = static_cast<double>(costs[0][i][j] - 10) / 40.0;
                    double y = (r >= 0) ? correl_val(x, r) : 1.0 - correl_val(x, -r);
                    costs[obj][i][j] = static_cast<unsigned int>(40 * y + 10);
                }
                else
                {
                    costs[obj][i][j] = rng.random(10, 50);
                }
            }
        }
    }

    // generate the capacity constraints
    unsigned int* cap = new unsigned int[n];
    for(unsigned int i=0; i<n; i++)
    {
        unsigned int totalres = 0;
        for(unsigned int j=0; j<m; j++)
        {
            totalres += res[i][j];
        }
        cap[i] = static_cast<unsigned int>(slack * totalres / n);
    }

    // write out the capacity constraints
    for(unsigned int i=0; i<n; i++)
    {
        cout << cap[i] << " ";
    }
    cout << endl << endl;

    // write out the resource allocations
    for(unsigned int i=0; i<n; i++)
    {
        for(unsigned int j=0; j<m; j++)
        {
            cout << setw(4) << res[i][j];
        }
        cout << endl;
    }
    cout << endl;

    // write out the cost matrices
    for(unsigned int obj=0; obj<k; obj++)
    {
        for(unsigned int i=0; i<n; i++)
        {
            for(unsigned int j=0; j<m; j++)
            {
                cout << setw(4) << costs[obj][i][j];
            }
            cout << endl;
        }
        cout << endl;
    }

    // free the memory from the data
    delete[] cap;
    for(unsigned int i=0; i<n; i++)
    {
        delete[] res[i];
    }
    delete[] res;
    for(unsigned int i=0; i<k; i++)
    {
        for(unsigned int j=0; j<n; j++)
        {
            delete[] costs[i][j];
        }
        delete[] costs[i];
    }
    delete[] costs;
}

void generate_type_d(unsigned int n, unsigned int m, unsigned int k, double slack, double r, double f)
{
    mtrandom rng;
    
    // generate the resource allocations
    unsigned int** res = new unsigned int*[n];
    for(unsigned int i=0; i<n; i++)
    {
        res[i] = new unsigned int[m];
        for(unsigned int j=0; j<m; j++)
        {
            res[i][j] = rng.random(1, 100);
        }
    }

    // generate the cost matrices
    unsigned int*** costs = new unsigned int**[k];
    for(unsigned int obj=0; obj<k; obj++)
    {
        costs[obj] = new unsigned int*[n];
        for(unsigned int i=0; i<n; i++)
        {
            costs[obj][i] = new unsigned int[m];
        }
    }
    
    for(unsigned int i=0; i<n; i++)
    {
        for(unsigned int j=0; j<m; j++)
        {
            for(unsigned int obj=0; obj<k; obj++)
            {
                if(obj==0)
                {
                    costs[obj][i][j] = 111 - res[i][j] + rng.random(-10, 10);
                }
                else if(rng.random() < f) 
                {
                    double x = static_cast<double>(costs[0][i][j] - 1) / 119.0;
                    double y = (r >= 0) ? correl_val(x, r) : 1.0 - correl_val(x, -r);
                    costs[obj][i][j] = static_cast<unsigned int>(119 * y + 1);
                }
                else
                {
                    costs[obj][i][j] = 111 - res[i][j] + rng.random(-10, 10);
                }
            }
        }
    }

    // generate the capacity constraints
    unsigned int* cap = new unsigned int[n];
    for(unsigned int i=0; i<n; i++)
    {
        unsigned int totalres = 0;
        for(unsigned int j=0; j<m; j++)
        {
            totalres += res[i][j];
        }
        cap[i] = static_cast<unsigned int>(slack * totalres / n);
    }

    // write out the capacity constraints
    for(unsigned int i=0; i<n; i++)
    {
        cout << cap[i] << " ";
    }
    cout << endl << endl;

    // write out the resource allocations
    for(unsigned int i=0; i<n; i++)
    {
        for(unsigned int j=0; j<m; j++)
        {
            cout << setw(4) << res[i][j];
        }
        cout << endl;
    }
    cout << endl;

    // write out the cost matrices
    for(unsigned int obj=0; obj<k; obj++)
    {
        for(unsigned int i=0; i<n; i++)
        {
            for(unsigned int j=0; j<m; j++)
            {
                cout << setw(4) << costs[obj][i][j];
            }
            cout << endl;
        }
        cout << endl;
    }

    // free the memory from the data
    delete[] cap;
    for(unsigned int i=0; i<n; i++)
    {
        delete[] res[i];
    }
    delete[] res;
    for(unsigned int i=0; i<k; i++)
    {
        for(unsigned int j=0; j<n; j++)
        {
            delete[] costs[i][j];
        }
        delete[] costs[i];
    }
    delete[] costs;
}

void generate_type_e(unsigned int n, unsigned int m, unsigned int k, double slack, double r, double f)
{
    mtrandom rng;
    
    // generate the resource allocations
    unsigned int** res = new unsigned int*[n];
    for(unsigned int i=0; i<n; i++)
    {
        res[i] = new unsigned int[m];
        for(unsigned int j=0; j<m; j++)
        {
            res[i][j] = static_cast<unsigned int>(1 - 10.0*log(rng.random()));
        }
    }

    // generate the cost matrices
    unsigned int*** costs = new unsigned int**[k];
    for(unsigned int obj=0; obj<k; obj++)
    {
        costs[obj] = new unsigned int*[n];
        for(unsigned int i=0; i<n; i++)
        {
            costs[obj][i] = new unsigned int[m];
        }
    }
    
    for(unsigned int i=0; i<n; i++)
    {
        for(unsigned int j=0; j<m; j++)
        {
            for(unsigned int obj=0; obj<k; obj++)
            {
                if(obj==0)
                {
                    costs[obj][i][j] = static_cast<unsigned int>(1000.0 / res[i][j] - 10.0 * rng.random());
                }
                else if(rng.random() < f) 
                {
                    double x = static_cast<double>(costs[0][i][j] - 10) / 990.0;
                    double y = (r >= 0) ? correl_val(x, r) : 1.0 - correl_val(x, -r);
                    costs[obj][i][j] = static_cast<unsigned int>(990 * y + 10);
                }
                else
                {
                    costs[obj][i][j] = static_cast<unsigned int>(1000.0 / res[i][j] - 10.0 * rng.random());
                }
            }
        }
    }

    // generate the capacity constraints
    unsigned int* cap = new unsigned int[n];
    for(unsigned int i=0; i<n; i++)
    {
        unsigned int totalres = 0;
        for(unsigned int j=0; j<m; j++)
        {
            totalres += res[i][j];
        }
        cap[i] = static_cast<unsigned int>(slack * totalres / n);
    }

    // write out the capacity constraints
    for(unsigned int i=0; i<n; i++)
    {
        cout << cap[i] << " ";
    }
    cout << endl << endl;

    // write out the resource allocations
    for(unsigned int i=0; i<n; i++)
    {
        for(unsigned int j=0; j<m; j++)
        {
            cout << setw(4) << res[i][j];
        }
        cout << endl;
    }
    cout << endl;

    // write out the cost matrices
    for(unsigned int obj=0; obj<k; obj++)
    {
        for(unsigned int i=0; i<n; i++)
        {
            for(unsigned int j=0; j<m; j++)
            {
                cout << setw(4) << costs[obj][i][j];
            }
            cout << endl;
        }
        cout << endl;
    }

    // free the memory from the data
    delete[] cap;
    for(unsigned int i=0; i<n; i++)
    {
        delete[] res[i];
    }
    delete[] res;
    for(unsigned int i=0; i<k; i++)
    {
        for(unsigned int j=0; j<n; j++)
        {
            delete[] costs[i][j];
        }
        delete[] costs[i];
    }
    delete[] costs;
}

double correl_val(double x, double r)
{
    double q;
    double w;
    double diff;
    double a;
    mtrandom rng;
    
    if(r==1)
    {
        return x;
    }
    if(r==0)
    {
        return rng.random();
    }
    
    do
    {
        q= rng.random();
        diff = q - x;
        if(diff<0)
        {
            diff = -diff;
        }
        w = exp(-(diff*diff) / (2.0 * (1.0-pow(r,0.5)) * (1.0-pow(r,0.5)))) / (1.0-pow(r,0.5)) * 2.506;
        
        a = rng.random() * 1.0 / (1.0-pow(r,0.5)) * 2.506;
    } while(w<a);
    
    return(q);
}

void generate_type_a(unsigned int n, unsigned int m, unsigned int k,
		     double slack, double r, double f)
{
    mtrandom rng;
    
    // generate the resource allocations
    unsigned int** res = new unsigned int*[n];
    for(unsigned int i=0; i<n; i++)
    {
        res[i] = new unsigned int[m];
        for(unsigned int j=0; j<m; j++)
        {
            res[i][j] = rng.random(5, 25);
        }
    }

    // generate the cost matrices
    unsigned int*** costs = new unsigned int**[k];
    for(unsigned int obj=0; obj<k; obj++)
    {
        costs[obj] = new unsigned int*[n];
        for(unsigned int i=0; i<n; i++)
        {
            costs[obj][i] = new unsigned int[m];
        }
    }
    
    for(unsigned int i=0; i<n; i++)
    {
        for(unsigned int j=0; j<m; j++)
        {
            for(unsigned int obj=0; obj<k; obj++)
            {
                if(obj==0)
                {
                    costs[obj][i][j] = rng.random(10, 100);
                }
                else if(rng.random() < f) 
                {
                    double x = static_cast<double>(costs[0][i][j] - 10) / 90.0;
                    double y = (r >= 0) ? correl_val(x, r) : 1.0 - correl_val(x, -r);
                    costs[obj][i][j] = static_cast<unsigned int>(90 * y + 10);
                }
                else
                {
                    costs[obj][i][j] = rng.random(10, 100);
                }
            }
        }
    }

    // generate the capacity constraints
    unsigned int* cap = new unsigned int[n];
    for(unsigned int i=0; i<n; i++)
    {
        unsigned int totalres = 0;
        for(unsigned int j=0; j<m; j++)
        {
            totalres += res[i][j];
        }
        cap[i] = static_cast<unsigned int>(slack * totalres / n);
    }

    // write out the capacity constraints
    for(unsigned int i=0; i<n; i++)
    {
        cout << cap[i] << " ";
    }
    cout << endl << endl;

    // write out the resource allocations
    for(unsigned int i=0; i<n; i++)
    {
        for(unsigned int j=0; j<m; j++)
        {
            cout << setw(4) << res[i][j];
        }
        cout << endl;
    }
    cout << endl;

    // write out the cost matrices
    for(unsigned int obj=0; obj<k; obj++)
    {
        for(unsigned int i=0; i<n; i++)
        {
            for(unsigned int j=0; j<m; j++)
            {
                cout << setw(4) << costs[obj][i][j];
            }
            cout << endl;
        }
        cout << endl;
    }

    // free the memory from the data
    delete[] cap;
    for(unsigned int i=0; i<n; i++)
    {
        delete[] res[i];
    }
    delete[] res;
    for(unsigned int i=0; i<k; i++)
    {
        for(unsigned int j=0; j<n; j++)
        {
            delete[] costs[i][j];
        }
        delete[] costs[i];
    }
    delete[] costs;
}

