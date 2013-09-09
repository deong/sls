/*!
 * \file ksgen.cpp
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

static const int MAX_WEIGHT = 50;
static const int MIN_WEIGHT = 1;
static const int MAX_VALUE = 100;
static const int MIN_VALUE = 1;

void print_usage();
void generate_instance(unsigned int n, unsigned int k);

int main(int argc, char** argv)
{
	if(argc == 1) {
		print_usage();
		exit(1);
	}

	unsigned int n = 0;
	unsigned int k = 0;
	unsigned int slack = 0;
	bool n_spec = false;
	bool k_spec = false;
	bool slack_spec = false;

	mtrandom::initialize();

	char curropt;
	while((curropt = getopt(argc, argv, "hn:k:s:")) != -1) {
		switch(curropt) {
		case 'n':
			n = static_cast<unsigned int>(atoi(optarg));
			n_spec = true;
			break;
		case 'k':
			k = static_cast<unsigned int>(atoi(optarg));
			k_spec = true;
			break;
		case 's':
			slack=atoi(optarg);
			slack_spec=true;
			break;
		case 'h':
			print_usage();
			break;
		default:
			cout << "invalid option specified" << endl;
			exit(1);
		}
	}

	if(!n_spec || !k_spec || !slack_spec) {
		print_usage();
		exit(1);
	}

	cout << ";; items: " << n << "\t";
	cout << "objectives: " << k << "\t";
	cout << "slack: " << slack << "\t";
	cout << endl;

	// knapsack capacity equal to average weight * slack factor
	unsigned int cap = static_cast<unsigned int>((double)(MAX_WEIGHT+MIN_WEIGHT)/2)*slack;
	cout << n << " " << k << " " << cap << endl;

	generate_instance(n,k);
	return 0;
}

void generate_instance(unsigned int n, unsigned int k)
{
	mtrandom rng;
	for(unsigned int i=0; i<n; ++i) {
		// generate random weight
		cout << rng.random(MIN_WEIGHT,MAX_WEIGHT) << " ";

		// generate random value for each objective
		for(unsigned int j=0; j<k; ++j) {
			cout << rng.random(MIN_VALUE,MAX_VALUE) << " ";
		}
		cout << endl;
	}
}

void print_usage()
{
	cout << "ksgen: generates instances of multiobjective 0-1 knapsack problem" << endl << endl;
	cout << "parameters: [default values, * if required]" << endl;
	cout << "\t-h: print this help screen" << endl;
	cout << "\t-n <items> [*]" << endl;
	cout << "\t-k <objectives> [*]" << endl;
	cout << "\t-s <slack> [*]" << endl;
	exit(1);
}

