/*!
 * \file problems.cpp
 *
 * defines the objective functions to be optimized by the search
 * algorithms
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <cstdlib>
#include <utility>
#include <cmath>
#include <iomanip>
#include <iterator>
#include "problems.h"
#include "mtrandom.h"
#include "kvparse/kvparse.h"
#include "keywords.h"
#include "utilities.h"

using namespace std;

/*!
 * \brief constructor
 */
problem::problem()
{
}

/*!
 * \brief destructor
 */
problem::~problem()
{
}

/*!
 * \brief override to perform problem specific initialization
 */
void problem::initialize()
{
}

/*!
 * \brief constructor
 */
numeric_problem::numeric_problem()
{
}

/*!
 * \brief destructor
 */
numeric_problem::~numeric_problem()
{
}

/*!
 * \brief constructor
 */
bit_string_problem::bit_string_problem()
{
}

/*!
 * \brief destructor
 */
bit_string_problem::~bit_string_problem()
{
}

/*!
 * \brief constructor
 */
permutation_problem::permutation_problem()
{
}

/*!
 * \brief destructor
 */
permutation_problem::~permutation_problem()
{
}

/*!
 * \brief constructor
 */
integer_problem::integer_problem()
{
}

/*!
 * \brief destructor
 */
integer_problem::~integer_problem()
{
}

/*!
 * \brief constructor
 */
qap_problem::qap_problem() :
	mqap_flow(0),
	mqap_dist(0),
	mqap_n(0),
	mqap_obj(0),
	mqap_delta(0)
{
}

/*!
 * \brief destructor
 */
qap_problem::~qap_problem()
{
	for(int i=0; i<mqap_n; i++) {
		delete[] mqap_dist[i];
	}
	delete[] mqap_dist;

	for(int k=0; k<mqap_obj; k++) {
		for(int i=0; i<mqap_n; i++) {
			delete[] mqap_flow[k][i];
			delete[] mqap_delta[k][i];
		}
		delete[] mqap_flow[k];
		delete[] mqap_delta[k];
	}
	delete[] mqap_flow;
	delete[] mqap_delta;
}

/*!
 * \brief read qap problem instance from data file
 */
void qap_problem::initialize()
{
	string filename;
	kvparse::parameter_value(keywords::PROBLEM_DATA, filename, true);

	ifstream in(filename.c_str());
	if(!in) {
		cerr << "could not open QAP data file: " << filename << endl;
		exit(1);
	}

	// check to see if first line is a comment
	char firstchar = static_cast<char>(in.peek());
	if(firstchar == ';' || firstchar == '#') {
		char comment[1024];
		in.getline(comment, 1024);
	}

	// get the problem size and the number of flow matrices
	string firstline;
	getline(in, firstline);
	istringstream iss(firstline);

	// read the tokens "facilities = n"
	string junk1, junk2;
	iss >> junk1;
	iss >> junk2;
	iss >> mqap_n;

	// read the tokens "objectives = k"
	iss >> junk1;
	iss >> junk2;
	iss >> mqap_obj;

	// create the distance and flow matrices
	mqap_dist = new int*[mqap_n];
	for(int i=0; i<mqap_n; i++) {
		mqap_dist[i] = new int[mqap_n];
		for(int j=0; j<mqap_n; j++) {
			in >> mqap_dist[i][j];
		}
	}

	// read the flow matrices
	mqap_flow = new int**[mqap_obj];
	for(int k=0; k<mqap_obj; k++) {
		mqap_flow[k] = new int*[mqap_n];
		for(int i=0; i<mqap_n; i++) {
			mqap_flow[k][i] = new int[mqap_n];
			for(int j=0; j<mqap_n; j++) {
				in >> mqap_flow[k][i][j];
			}
		}
	}

	// construct the delta matrix
	mqap_delta = new permutation_problem::FitnessType**[mqap_obj];
	for(int i=0; i<mqap_obj; i++) {
		mqap_delta[i] = new permutation_problem::FitnessType*[mqap_n];
		for(int j=0; j<mqap_n; j++) {
			mqap_delta[i][j] = new permutation_problem::FitnessType[mqap_n];
		}
	}
}

/*!
 * \brief return the size of the problem instance
 */
unsigned int qap_problem::dimensions() const
{
	return mqap_n;
}

/*!
 * \brief return the number of objectives
 */
unsigned int qap_problem::objectives() const
{
	return mqap_obj;
}

/*!
 * \brief compute delta_f (change in fitness) given a potential swap
 */
void qap_problem::compute_delta(vector<int>& p, int i, int j)
{
	long d;
	int k;
	for(int x=0; x<mqap_obj; x++) {
		d = (mqap_dist[i][i]-mqap_dist[j][j]) *
		    (mqap_flow[x][p[j]][p[j]]-mqap_flow[x][p[i]][p[i]]) +
		    (mqap_dist[i][j]-mqap_dist[j][i]) *
		    (mqap_flow[x][p[j]][p[i]]-mqap_flow[x][p[i]][p[j]]);
		for(k=0; k<mqap_n; k++) {
			if(k!=i && k!=j) {
				d = d + (mqap_dist[k][i]-mqap_dist[k][j]) *
				    (mqap_flow[x][p[k]][p[j]]-mqap_flow[x][p[k]][p[i]]) +
				    (mqap_dist[i][k]-mqap_dist[j][k]) *
				    (mqap_flow[x][p[j]][p[k]]-mqap_flow[x][p[i]][p[k]]);
			}
		}
		mqap_delta[x][i][j] = d;
	}
}

/*!
 * \brief compute fast delta
 */
void qap_problem::compute_delta_part(vector<int>& p, int i, int j, int r, int s)
{
	for(int x=0; x<mqap_obj; x++) {
		mqap_delta[x][i][j] = mqap_delta[x][i][j] +
		                      (mqap_dist[r][i]- mqap_dist[r][j] + mqap_dist[s][j] - mqap_dist[s][i]) *
		                      (mqap_flow[x][p[s]][p[i]] - mqap_flow[x][p[s]][p[j]] +
		                       mqap_flow[x][p[r]][p[j]] - mqap_flow[x][p[r]][p[i]]) +
		                      (mqap_dist[i][r] - mqap_dist[j][r] + mqap_dist[j][s] - mqap_dist[i][s]) *
		                      (mqap_flow[x][p[i]][p[s]] - mqap_flow[x][p[j]][p[s]] +
		                       mqap_flow[x][p[j]][p[r]] - mqap_flow[x][p[i]][p[r]]);
	}
}

/*!
 * \brief compute the fitness of the given permutation
 */
bool qap_problem::evaluate(const vector<int>& p, vector<permutation_problem::FitnessType>& fit) const
{
	for(int i=0; i<mqap_obj; i++) {
		fit[i] = 0;
		for(int j=0; j<mqap_n; j++) {
			for(int k=0; k<mqap_n; k++) {
				fit[i] += mqap_dist[j][k] * mqap_flow[i][p[j]][p[k]];
			}
		}
	}
	return true;
}

/*!
 * \brief update the delta matrix after accepting a move
 */
void qap_problem::update(vector<int>& p, int pos1, int pos2)
{
	for(int i=0; i<mqap_n-1; i++) {
		for(int j=i+1; j<mqap_n; j++) {
			if(i!=pos1 && i!=pos2 && j!=pos1 && j!=pos2) {
				compute_delta_part(p,i,j,pos1,pos2);
			} else {
				compute_delta(p,i,j);
			}
		}
	}
}

/*!
 * \brief return the delta_f given a potential move
 */
void qap_problem::get_delta(int i, int j, vector<permutation_problem::FitnessType>& del)
{
	for(int k=0; k<mqap_obj; k++) {
		del[k] = mqap_delta[k][i][j];
	}
}

/*!
 * \brief constructor
 */
gap_problem::gap_problem() :
	agents(0),
	tasks(0),
	mgap_obj(0),
	capacity(0),
	resources(0),
	cost(0)
{
}

/*!
 * \brief destructor
 */
gap_problem::~gap_problem()
{
	for(unsigned int k=0; k<mgap_obj; k++) {
		for(unsigned int i=0; i<agents; i++) {
			delete[] cost[k][i];
		}
		delete[] cost[k];
	}
	delete[] cost;

	for(unsigned int i=0; i<agents; i++) {
		delete[] resources[i];
	}
	delete[] resources;
	delete[] capacity;
}

/*!
 * \brief initialize the gap instance from a data file
 */
void gap_problem::initialize()
{
	string filename;
	kvparse::parameter_value(keywords::PROBLEM_DATA, filename, true);

	ifstream in(filename.c_str());
	if(!in) {
		cerr << "could not open GAP data file: " << filename << endl;
		exit(1);
	}

	// check to see if first line is a comment
	char firstchar = static_cast<char>(in.peek());
	if(firstchar == ';' || firstchar == '#') {
		string comment;
		getline(in, comment);
	}

	// get the numbers of agents, tasks, and objectives
	in >> agents >> tasks >> mgap_obj;

	// read in the resource capacity of each agent
	capacity = new int[agents];
	for(unsigned int i=0; i<agents; i++) {
		in >> capacity[i];
	}

	// get the resource requirements for each agent/task pair
	resources = new int*[agents];
	for(unsigned int i=0; i<agents; i++) {
		resources[i] = new int[tasks];
		for(unsigned int j=0; j<tasks; j++) {
			in >> resources[i][j];
		}
	}

	// get the cost matrix for each objective
	cost = new int**[mgap_obj];
	for(unsigned int k=0; k<mgap_obj; k++) {
		cost[k] = new int*[agents];
		for(unsigned int i=0; i<agents; i++) {
			cost[k][i] = new int[tasks];
			for(unsigned int j=0; j<tasks; j++) {
				in >> cost[k][i][j];
			}
		}
	}
}

/*!
 * \brief return the number of tasks in the instance
 */
unsigned int gap_problem::dimensions() const
{
	return tasks;
}

/*!
 * \brief return the number of objectives
 */
unsigned int gap_problem::objectives() const
{
	return mgap_obj;
}

/*!
 * \brief get list of acceptable values for a given allele
 */
void gap_problem::legal_values(unsigned int index, vector<int>& val) const
{
	val.clear();
	for(int i=0; i<int(agents); i++) {
		val.push_back(i);
	}
}

/*!
 * \brief evaluate the fitness of a candidate solution
 */
bool gap_problem::evaluate(const vector<int>& p, vector<int>& fit) const
{
	bool feas = true;

	// set up the agent mapping and set the capacity used to zero
	vector<int> capused;
	capused.assign(agents, 0);

	// compute the resources used by each agent
	for(unsigned int i=0; i<tasks; i++) {
		capused[p[i]] += resources[p[i]][i];

		if(capused[p[i]] > capacity[p[i]]) {
			feas = false;
			break;
		}
	}

	for(unsigned int k=0; k<mgap_obj; k++) {
		fit[k] = 0;
		for(unsigned int i=0; i<tasks; i++) {
			fit[k] += cost[k][p[i]][i];
		}
	}

	return feas;
}

/*!
 * \brief constructor
 */
gsap_problem::gsap_problem() :
	gsap_obj(0),
	gsap_dim(0),
	n_agents(0),
	n_days(0),
	n_time_slots(0)
{
}

/*!
 * \brief destructor
 */
gsap_problem::~gsap_problem()
{
	m_capacities.clear();

	for(unsigned int i=0; i<n_days; ++i) {
		for(unsigned int j=0; j<n_time_slots; ++j) {
			m_tasks[i][j].clear();
		}
		m_tasks[i].clear();
	}
	m_tasks.clear();

	for(unsigned int i=0; i<n_agents; ++i) {
		m_resources[i].clear();
	}
	m_resources.clear();

	for(unsigned int k=0; k<gsap_obj; ++k) {
		for(unsigned int i=0; i<n_agents; ++i) {
			m_costs[k][i].clear();
		}
		m_costs[k].clear();
	}
	m_costs.clear();
}

/**
 * \brief return number of task elements
 */
unsigned int gsap_problem::dimensions() const
{
	return gsap_dim;
}

/**
 * \brief return the number of objectives
 */
unsigned int gsap_problem::objectives() const
{
	return gsap_obj;
}

/**
 * \brief return the number of agents in the problem
 */
unsigned int gsap_problem::agents() const
{
	return n_agents;
}

/**
 * \brief return the number of tasks in the problem
 */
unsigned int gsap_problem::tasks() const
{
	return n_tasks;
}

/*!
 * \brief initialize the problem from a given data file
 */
void gsap_problem::initialize()
{
	string filename;
	kvparse::parameter_value(keywords::PROBLEM_DATA, filename, true);

	ifstream in(filename.c_str());
	if(!in) {
		cerr << "could not open GAP data file: " << filename << endl;
		exit(1);
	}

	// check to see if first line is a comment
	char firstchar = static_cast<char>(in.peek());
	if(firstchar == ';' || firstchar == '#') {
		string comment;
		getline(in, comment);
	}

	// get basic problem size information
	in >> n_agents >> n_tasks >> gsap_obj >> n_days >> n_time_slots;

	//! read in capacity of each agent
	m_capacities.resize(n_agents);
	for(unsigned int i=0; i<n_agents; ++i) {
		in >> m_capacities[i];
	}

	//! read in the task requirements for each time slot
	gsap_dim = 0;
	m_elements.clear();
	m_tasks.resize(n_days);
	for(unsigned int i=0; i<n_days; ++i) {
		m_tasks[i].resize(n_time_slots);
		for(unsigned int j=0; j<n_time_slots; ++j) {
			string line;
			do {
				getline(in, line);
			} while(line.find_first_of("0123456789")==string::npos);

			tsb.push_back(gsap_dim);
			vector<unsigned int> tmp;
			istringstream iss(line);
			while(true) {
				int v;
				iss >> v;
				if(iss.eof()) {
					break;
				}
				m_tasks[i][j].push_back(v);
				m_elements.push_back(task_element(i,j,v));
				++gsap_dim;
			}
		}
	}
	tsb.push_back(gsap_dim);

	//! read in the resource matrix
	m_agents_for_task.resize(n_tasks);
	m_resources.resize(n_agents);
	for(unsigned int i=0; i<n_agents; ++i) {
		m_resources[i].resize(n_tasks);
		for(unsigned int j=0; j<n_tasks; ++j) {
			in >> m_resources[i][j];
			if(m_resources[i][j]>0) {
				m_agents_for_task[j].push_back(i);
			}
		}
	}

	//! read in each of the $k$ cost matrices
	m_costs.resize(gsap_obj);
	for(unsigned int k=0; k<gsap_obj; ++k) {
		m_costs[k].resize(n_agents);
		for(unsigned int i=0; i<n_agents; ++i) {
			m_costs[k][i].resize(n_tasks);
			for(unsigned int j=0; j<n_tasks; ++j) {
				in >> m_costs[k][i][j];
			}
		}
	}
}

/**
 * \brief get the list of task elements
 */
const vector<gsap_problem::task_element>& gsap_problem::get_elements() const
{
	return m_elements;
}

/**
 * \brief return the lists of agents for each task
 */
const vector<vector<int> >& gsap_problem::get_agent_task_map() const
{
	return m_agents_for_task;
}

/**
 * \brief return a vector of time slot boundaries
 */
const vector<unsigned int>& gsap_problem::time_slot_boundaries() const
{
	return tsb;
}

/*!
 * \brief return the capacity vector
 */
const vector<unsigned int>& gsap_problem::capacities() const
{
	return m_capacities;
}

/*!
 * \brief return the resource matrix
 */
const vector<vector<unsigned int> >& gsap_problem::resources() const
{
	return m_resources;
}

/**
 * \brief return list of agents who can perform a given task
 */
void gsap_problem::legal_values(unsigned int index, vector<int>& agts) const
{
	copy(m_agents_for_task[m_elements[index].task].begin(),
	     m_agents_for_task[m_elements[index].task].end(),
	     back_inserter(agts));
}

/**
 * \brief evaluate a candidate solution to the gsap problem
 *
 * \todo handle constraints
 */
bool gsap_problem::evaluate(const vector<int>& p, vector<int>& fit) const
{
	vector<unsigned int> used;
	used.assign(n_agents,0);
	int unmatched=0;
	for(unsigned int i=0; i<p.size(); ++i) {
		if(p[i]==-1) {
			unmatched++;
			continue;
		}
		used[p[i]]+=m_resources[p[i]][m_elements[i].task];
	}

	for(unsigned int k=0; k<fit.size(); ++k) {
		fit[k] = 0;
		for(unsigned int i=0; i<p.size(); ++i) {
			if(p[i]==-1) {
				continue;
			}
			fit[k] += m_costs[k][p[i]][m_elements[i].task];
		}
	}

	for(unsigned int i=0; i<used.size(); ++i) {
		if(used[i]>m_capacities[i]) {
			return false;
		}
	}
	return true;
}

/*!
 * \brief constructor
 */
knapsack_problem::knapsack_problem() :
	ks_obj(0),
	ks_dim(0),
	ks_capacity(0)
{
}

/*!
 * \brief destructor
 */
knapsack_problem::~knapsack_problem()
{
}

/*!
 * \brief initialize the knapsack problem instance
 */
void knapsack_problem::initialize()
{
	string filename;
	kvparse::parameter_value(keywords::PROBLEM_DATA, filename, true);

	ifstream in(filename.c_str());
	if(!in) {
		cerr << "could not open knapsack data file: " << filename << endl;
		exit(1);
	}

	// check to see if first line is a comment
	char firstchar = static_cast<char>(in.peek());
	if(firstchar == ';' || firstchar == '#') {
		string comment;
		getline(in, comment);
	}

	in >> ks_dim >> ks_obj >> ks_capacity;
	ks_weights.resize(ks_dim);
	ks_values.resize(ks_dim);
	for(unsigned int i=0; i<ks_dim; ++i) {
		in >> ks_weights[i];
		ks_values[i].resize(ks_obj);
		for(unsigned int j=0; j<ks_obj; ++j) {
			in >> ks_values[i][j];
		}
	}
}

/*!
 * \brief return the number of objectives
 */
unsigned int knapsack_problem::objectives() const
{
	return ks_obj;
}

/*!
 * \brief return the number of possible items
 */
unsigned int knapsack_problem::dimensions() const
{
	return ks_dim;
}

/*!
 * \brief evaluate a candidate solution to the knapsack problem
 *
 * Note that all objective function values are negated so that minimization
 * is appropriate.
 */
bool knapsack_problem::evaluate(const vector<int>& p, vector<int>& fit) const
{
	int total_weight=0;
	for(unsigned int i=0; i<ks_obj; ++i) {
		fit[i]=0;
		for(unsigned int j=0; j<ks_dim; ++j) {
			fit[i]+=p[j]*ks_values[j][i];
			total_weight+=p[j]*ks_weights[j];
		}
		// negate the fitness value for minimization
		fit[i] = -fit[i];
	}
	return (static_cast<unsigned int>(total_weight)<=ks_capacity);
}

/*!
 * \brief constructor
 */
onemax_problem::onemax_problem() :
	onemax_dim(0)
{
}

/*!
 * \brief destructor
 */
onemax_problem::~onemax_problem()
{
}

/*!
 * \brief initialize the number of bits in the problem
 */
void onemax_problem::initialize()
{
	kvparse::parameter_value(keywords::DIMENSIONS, onemax_dim, true);
}

/*!
 * \brief return the problem size
 */
unsigned int onemax_problem::dimensions() const
{
	return onemax_dim;
}

/*!
 * \brief onemax is a single objective problem
 */
unsigned int onemax_problem::objectives() const
{
	return 1;
}

/*!
 * \brief evaluate the fitness of a candidate solution
 */
bool onemax_problem::evaluate(const vector<int>& p, vector<int>& fit) const
{
	int ones = 0;
	for(unsigned int i=0; i<p.size(); i++) {
		ones += p[i];
	}
	fit[0] = ones;
	return true;
}

/*!
 * \brief constructor
 */
lotz_problem::lotz_problem() :
	lotz_dim(0)
{
}

/*!
 * \brief destructor
 */
lotz_problem::~lotz_problem()
{
}

/*!
 * \brief initialize the problem size
 */
void lotz_problem::initialize()
{
	kvparse::parameter_value(keywords::DIMENSIONS, lotz_dim, true);
}

/*!
 * \brief return the number of bits in the problem
 */
unsigned int lotz_problem::dimensions() const
{
	return lotz_dim;
}

/*!
 * \brief lotz is a bi-objective problem
 */
unsigned int lotz_problem::objectives() const
{
	return 2;
}

/*!
 * \brief evaluate the fitness of a candidate solution
 */
bool lotz_problem::evaluate(const vector<int>& p, vector<int>& fit) const
{
	fit[0] = 0;
	for(unsigned int i=0; i<p.size(); i++) {
		if(p[i] == 1) {
			fit[0]++;
		} else {
			break;
		}
	}
	fit[0] = p.size() - fit[0];

	fit[1] = 0;
	for(int i=int(p.size())-1; i>=0; i--) {
		if(p[i] == 0) {
			fit[1]++;
		} else {
			break;
		}
	}
	fit[1] = p.size() - fit[1];
	return true;
}

/*!
 * \brief constructor
 */
f1_problem::f1_problem() :
	f1_dim(0)
{
}

/*!
 * \brief destructor
 */
f1_problem::~f1_problem()
{
}

/*!
 * \brief initialize the number of dimensions in the problem
 */
void f1_problem::initialize()
{
	kvparse::parameter_value(keywords::DIMENSIONS, f1_dim, true);
}

/*!
 * \brief return the number of dimensions in the problem
 */
unsigned int f1_problem::dimensions() const
{
	return f1_dim;
}

/*!
 * \brief f1 is a single-objective problem
 */
unsigned int f1_problem::objectives() const
{
	return 1;
}

/*!
 * \brief bound the legal parameter values
 */
pair<double,double> f1_problem::parameter_range(unsigned int index) const
{
	return make_pair<double,double>(-5.12,5.11);
}

/*!
 * \brief evaluate the fitness of a candidate solution
 */
bool f1_problem::evaluate(const vector<double>& p, vector<double>& fit) const
{
	fit[0] = 0.0;
	for(unsigned int i=0; i<p.size(); i++) {
		fit[0] += p[i]*p[i];
	}
	return true;
}

/*!
 * \brief constructor
 */
f2_problem::f2_problem() :
	f2_dim(2)
{
}

/*!
 * \brief destructor
 */
f2_problem::~f2_problem()
{
}

/*!
 * \brief initialize the number of dimensions in the problem
 */
void f2_problem::initialize()
{
}

/*!
 * \brief return the number of dimensions in the problem
 */
unsigned int f2_problem::dimensions() const
{
	return f2_dim;
}

/*!
 * \brief f2 is a single-objective problem
 */
unsigned int f2_problem::objectives() const
{
	return 1;
}

/*!
 * \brief bound the legal parameter values
 */
pair<double,double> f2_problem::parameter_range(unsigned int index) const
{
	return make_pair<double,double>(-2.048, 2.047);
}

/*!
 * \brief evaluate the fitness of a candidate solution
 */
bool f2_problem::evaluate(const vector<double>& p, vector<double>& fit) const
{
	double x1=p[0];
	double x2=p[1];
	double sq_x1=x1*x1;
	double diff_x1=1.0-x1;

	fit[0] = 100.0 * ((sq_x1-x2) * (sq_x1-x2)) + (diff_x1 * diff_x1);
	return true;
}

/*!
 * \brief constructor
 */
f3_problem::f3_problem() :
	f3_dim(5)
{
}

/*!
 * \brief destructor
 */
f3_problem::~f3_problem()
{
}

/*!
 * \brief initialize the number of dimensions in the problem
 */
void f3_problem::initialize()
{
}

/*!
 * \brief return the number of dimensions in the problem
 */
unsigned int f3_problem::dimensions() const
{
	return f3_dim;
}

/*!
 * \brief f3 is a single-objective problem
 */
unsigned int f3_problem::objectives() const
{
	return 1;
}

/*!
 * \brief bound the legal parameter values
 */
pair<double,double> f3_problem::parameter_range(unsigned int index) const
{
	return make_pair<double,double>(-5.12,5.11);
}

/*!
 * \brief evaluate the fitness of a candidate solution
 */
bool f3_problem::evaluate(const vector<double>& p, vector<double>& fit) const
{
	double sum = 0.0;
	for(unsigned int i=0; i<f3_dim; ++i) {
		double part = (int)p[i];
		if(part > p[i]) {
			part -= 1;
		}
		sum += part;
	}
	fit[0] = 30 + sum;
	return true;
}

/*!
 * \brief constructor
 */
f4_problem::f4_problem() :
	f4_dim(30)
{
}

/*!
 * \brief destructor
 */
f4_problem::~f4_problem()
{
}

/*!
 * \brief initialize the number of dimensions in the problem
 */
void f4_problem::initialize()
{
}

/*!
 * \brief return the number of dimensions in the problem
 */
unsigned int f4_problem::dimensions() const
{
	return f4_dim;
}

/*!
 * \brief f4 is a single-objective problem
 */
unsigned int f4_problem::objectives() const
{
	return 1;
}

/*!
 * \brief bound the legal parameter values
 */
pair<double,double> f4_problem::parameter_range(unsigned int index) const
{
	return make_pair<double,double>(-1.28, 1.27);
}

/*!
 * \brief evaluate the fitness of a candidate solution
 */
bool f4_problem::evaluate(const vector<double>& p, vector<double>& fit) const
{
	double sum = 0.0;
	mtrandom rng;
	for(unsigned int i=0; i<f4_dim; ++i) {
		double quad=p[i] *p[i] * p[i] * p[i];
		quad *= (i+1);
		sum += quad + rng.gaussian();
	}
	fit[0] = sum;
	return true;
}

/*!
 * \brief constructor
 */
f5_problem::f5_problem() :
	f5_dim(2)
{
}

/*!
 * \brief destructor
 */
f5_problem::~f5_problem()
{
}

/*!
 * \brief initialize the number of dimensions in the problem
 */
void f5_problem::initialize()
{
}

/*!
 * \brief return the number of dimensions in the problem
 */
unsigned int f5_problem::dimensions() const
{
	return f5_dim;
}

/*!
 * \brief f5 is a single-objective problem
 */
unsigned int f5_problem::objectives() const
{
	return 1;
}

/*!
 * \brief bound the legal parameter values
 */
pair<double,double> f5_problem::parameter_range(unsigned int index) const
{
	return make_pair<double,double>(-65.536, 65.535);
}

/*!
 * \brief evaluate the fitness of a candidate solution
 */
bool f5_problem::evaluate(const vector<double>& p, vector<double>& fit) const
{
	static int a[2][25] = {
		{
			-32, -16,   0,  16,  32, -32, -16,   0,  16,  32, -32, -16,
			0,  16,  32, -32, -16,   0,  16,  32, -32, -16,   0,  16, 32
		},
		{
			-32, -32, -32, -32, -32, -16, -16, -16, -16, -16,  16,  16, 16,
			16,  16,  32,  32,  32,  32,  32
		}
	};

	double lowtot, prod, total = 0.002;

	for(unsigned int j=0; j<25; j+=1) {
		lowtot=1.0 + (double)j;
		for(unsigned int i=0; i<2; i+=1) {
			prod = 1.0;
			for(int power=0; power<6; power+=1) {
				prod *= p[i] - a[i][j];
			}
			lowtot += prod;
		}
		total += 1.0/lowtot;
	}

	fit[0] = 500.0 - (1.0/total);
	return true;
}

/*!
 * \brief constructor
 */
kur_problem::kur_problem()
{
}

/*!
 * \brief destructor
 */
kur_problem::~kur_problem()
{
}

/*!
 * \brief kursawe's function is defined to be 3-dimensional
 */
unsigned int kur_problem::dimensions() const
{
	return 3;
}

/*!
 * \brief kursawe's function is bi-objective
 */
unsigned int kur_problem::objectives() const
{
	return 2;
}

/*!
 * \brief bound the legal parameter values
 */
pair<double,double> kur_problem::parameter_range(unsigned int index) const
{
	return make_pair(-5.0,5.0);
}

/*!
 * \brief evaluate the fitness of a candidate solution
 */
bool kur_problem::evaluate(const vector<double>& p, vector<double>& fit) const
{
	fit[0] = fit[1] = 0.0;
	for(unsigned int i=0; i<p.size()-1; i++) {
		fit[0] += -10.0 * exp(-0.2*sqrt(p[i]*p[i] + p[i+1]*p[i+1]));
		fit[1] += pow(fabs(p[i]), 0.8) + 5.0*sin(pow(p[i],3.0));
	}
	fit[1] += pow(fabs(p[p.size()-1]),0.8) + 5.0*sin(pow(p[p.size()-1],3.0));
	return true;
}

/*!
 * \brief constructor
 */
sch_problem::sch_problem()
{
}

/*!
 * \brief destructor
 */
sch_problem::~sch_problem()
{
}

/*!
 * \brief defined to be 1-dimensional
 */
unsigned int sch_problem::dimensions() const
{
	return 1;
}

/*!
 * \brief defined as bi-objective function
 */
unsigned int sch_problem::objectives() const
{
	return 2;
}

/*!
 * \brief bound the legal parameter values
 */
pair<double,double> sch_problem::parameter_range(unsigned int index) const
{
	return make_pair(-1000.0,1000.0);
}

/*!
 * \brief evaluate the fitness of a candidate solution
 */
bool sch_problem::evaluate(const vector<double>& p, vector<double>& fit) const
{
	fit[0] = p[0] * p[0];
	fit[1] = (p[0]-2) * (p[0]-2);
	return true;
}

/*!
 * \brief constructor
 */
zdt1_problem::zdt1_problem()
{
}

/*!
 * \brief destructor
 */
zdt1_problem::~zdt1_problem()
{
}

/*!
 * \brief defined as a 30-dimensional problem
 */
unsigned int zdt1_problem::dimensions() const
{
	return 30;
}

/*!
 * \brief defined as bi-objective problem
 */
unsigned int zdt1_problem::objectives() const
{
	return 2;
}

/*!
 * \brief bound the legal parameter values
 */
pair<double,double> zdt1_problem::parameter_range(unsigned int index) const
{
	return make_pair(0.0, 1.0);
}

/*!
 * \brief evaluate the fitness of a candidate solution
 */
bool zdt1_problem::evaluate(const vector<double>& p, vector<double>& fit) const
{
	fit[0] = p[0];
	double gx = 0.0;
	for(unsigned int i=1; i<p.size(); i++) {
		gx += p[i];
	}
	gx = 1 + 9 * gx / (p.size()-1);
	fit[1] = gx * (1 - sqrt(p[0]/gx));
	return true;
}

/*!
 * \brief constructor
 */
zdt2_problem::zdt2_problem()
{
}

/*!
 * \brief destructor
 */
zdt2_problem::~zdt2_problem()
{
}

/*!
 * \brief defined as a 30-dimensional problem
 */
unsigned int zdt2_problem::dimensions() const
{
	return 30;
}

/*!
 * \brief defined as a bi-objective problem
 */
unsigned int zdt2_problem::objectives() const
{
	return 2;
}

/*!
 * \brief bound the legal parameter values
 */
pair<double,double> zdt2_problem::parameter_range(unsigned int index) const
{
	return make_pair(0.0, 1.0);
}

/*!
 * \brief evaluate the fitness of a candidate solution
 */
bool zdt2_problem::evaluate(const vector<double>& p, vector<double>& fit) const
{
	fit[0] = p[0];
	double gx = 0.0;
	for(unsigned int i=1; i<p.size(); i++) {
		gx += p[i];
	}
	gx = 1 + 9 * gx / (p.size()-1);
	fit[1] = gx * (1 - pow(p[0]/gx, 2.0));
	return true;
}

/*!
 * \brief constructor
 */
zdt3_problem::zdt3_problem()
{
}

/*!
 * \brief destructor
 */
zdt3_problem::~zdt3_problem()
{
}

/*!
 * \brief defined as a 30-dimensional problem
 */
unsigned int zdt3_problem::dimensions() const
{
	return 30;
}

/*!
 * \brief defined as a bi-objective problem
 */
unsigned int zdt3_problem::objectives() const
{
	return 2;
}

/*!
 * \brief bound the legal parameter values
 */
pair<double,double> zdt3_problem::parameter_range(unsigned int index) const
{
	return make_pair(0.0, 1.0);
}

/*!
 * \brief evaluate the fitness of a candidate solution
 */
bool zdt3_problem::evaluate(const vector<double>& p, vector<double>& fit) const
{
	fit[0] = p[0];
	double gx = 0.0;
	for(unsigned int i=1; i<p.size(); i++) {
		gx += p[i];
	}
	gx = 1 + 9 * gx / (p.size()-1);
	fit[1] = gx * (1 - sqrt(p[0]/gx) - p[0]/gx * sin(10*M_PI*p[0]));
	return true;
}

/*!
 * \brief constructor
 */
zdt4_problem::zdt4_problem()
{
}

/*!
 * \brief destructor
 */
zdt4_problem::~zdt4_problem()
{
}

/*!
 * \brief defined as a 10-dimensional problem
 */
unsigned int zdt4_problem::dimensions() const
{
	return 10;
}

/*!
 * \brief defined as a bi-objective problem
 */
unsigned int zdt4_problem::objectives() const
{
	return 2;
}

/*!
 * \brief bound the legal parameter values
 */
pair<double,double> zdt4_problem::parameter_range(unsigned int index) const
{
	if(index == 0) {
		return make_pair(0.0, 1.0);
	} else {
		return make_pair(-5.0, 5.0);
	}
}

/*!
 * \brief evaluate the fitness of a candidate solution
 */
bool zdt4_problem::evaluate(const vector<double>& p, vector<double>& fit) const
{
	double f1 = p[0];
	double g = 91.0;
	for (unsigned i=1; i<p.size(); i++) {
		g += p[i] * p[i] - 10.0 * cos(4.0 * M_PI * p[i]);
	}
	double h = 1.0 - sqrt(f1/g);
	double f2 = g * h;
	fit[0] = f1;
	fit[1] = f2;
	return true;
}

/*!
 * \brief constructor
 */
zdt6_problem::zdt6_problem()
{
}

/*!
 * \brief destructor
 */
zdt6_problem::~zdt6_problem()
{
}

/*!
 * \brief defined as a 10-dimensional problem
 */
unsigned int zdt6_problem::dimensions() const
{
	return 10;
}

/*!
 * \brief defined as a bi-objective problem
 */
unsigned int zdt6_problem::objectives() const
{
	return 2;
}

/*!
 * \brief bound the legal parameter values
 */
pair<double,double> zdt6_problem::parameter_range(unsigned int index) const
{
	return make_pair(0.0, 1.0);
}

/*!
 * \brief evaluate the fitness of a candidate solution
 */
bool zdt6_problem::evaluate(const vector<double>& p, vector<double>& fit) const
{
	fit[0] = 1.0 - exp(-4 * p[0]) * pow(sin(6 * M_PI * p[0]), 6.0);
	double gx = 0.0;
	for(unsigned int i=1; i<p.size(); i++) {
		gx += p[i];
	}
	gx = 1 + 9 * pow(gx / (p.size()-1), 0.25);
	fit[1] = gx * (1 - pow(fit[0]/gx, 2.0));
	return true;
}

/*!
 * \brief constructor
 */
dtlz1_problem::dtlz1_problem() :
	dtlz1_dim(0),
	dtlz1_obj(0)
{
}

/*!
 * \brief destructor
 */
dtlz1_problem::~dtlz1_problem()
{
}

/*!
 * \brief return the dimensionality of the problem instance
 */
unsigned int dtlz1_problem::dimensions() const
{
	return dtlz1_dim;
}

/*!
 * \brief return the number of objectives in the problem instance
 */
unsigned int dtlz1_problem::objectives() const
{
	return dtlz1_obj;
}

/*!
 * \brief initialize the problem instance
 */
void dtlz1_problem::initialize()
{
	dtlz1_obj = 3;
	dtlz1_dim = 7;
	kvparse::parameter_value(keywords::DIMENSIONS, dtlz1_dim, false);
	kvparse::parameter_value(keywords::OBJECTIVES, dtlz1_obj, false);
}

/*!
 * \brief bound the legal parameter values
 */
pair<double,double> dtlz1_problem::parameter_range(unsigned int index) const
{
	return make_pair(0.0, 1.0);
}

/*!
 * \brief evaluate the fitness of a candidate solution
 */
bool dtlz1_problem::evaluate(const vector<double>& p, vector<double>& fit) const
{
	double g = 0.0;
	unsigned int n = dtlz1_dim;
	unsigned int k = n - dtlz1_obj + 1;

	for(unsigned int i=n-k+1; i<=n; i++) {
		g += pow(p[i-1]-0.5,2) - cos(20.0 * M_PI * (p[i-1]-0.5));
	}
	g = 100 * (k+g);

	for(unsigned int i=1; i<=dtlz1_obj; i++) {
		double f = 0.5 * (1+g);
		for(unsigned int j=dtlz1_obj-i; j>=1; j--) {
			f *= p[j-1];
		}
		if(i>1) {
			f *= 1 - p[dtlz1_obj-i+1-1];
		}
		fit[i-1] = f;
	}
	return true;
}

/**
 * @brief constructor
 */
rana_problem::rana_problem()
{
}

/**
 * @brief destructor
 */
rana_problem::~rana_problem()
{
}

/**
 * @brief return the number of variables
 */
unsigned int rana_problem::dimensions() const
{
	return 10;
}

/**
 * @brief return the number of objectives
 */
unsigned int rana_problem::objectives() const
{
	return 1;
}

/**
 * @brief return the valid range of values for the input parameter
 */
pair<double,double> rana_problem::parameter_range(unsigned int index) const
{
	return make_pair<double,double>(-512.0, 511.0);
}

/**
 * @brief do the evaluation
 */
bool rana_problem::evaluate(const vector<double>& p, vector<double>& fit) const
{
	const double RANA_WEIGHTS[]= {0.3489, 0.1848, 0.3790, 0.4386, 0.9542, 0.1430, 0.7849, 0.3689, 0.9767, 0.8163};
	const double RANA_SUM_WEIGHTS=5.3953;
	double sum=0.0;
	for(unsigned int i=0; i<p.size(); i++) {
		double p1 = p[i];
		double p2 = p[(i+1)%p.size()];

		sum+=RANA_WEIGHTS[i] * (p1 * sin(sqrt(fabs(p2+1.0-p1)))*
		                        cos(sqrt(fabs(p1+p2+1.0))) +
		                        (p2+1.0)*
		                        cos(sqrt(fabs(p2+1.0-p1)))*
		                        sin(sqrt(fabs(p1+p2+1.0))));
	}
	fit[0] = sum / RANA_SUM_WEIGHTS;
	return true;
}

/*!
 * \brief create a permutation problem
 */
permutation_problem* permutation_problem_factory::construct()
{
	string prob;
	kvparse::parameter_value(keywords::PROBLEM, prob, true);

	if(prob == "qap" || prob == "mqap") {
		permutation_problem* p = new qap_problem;
		p->initialize();
		return p;
	} else {
		cerr << "invalid problem specified: " << prob << endl;
		exit(1);
		return 0;
	}
}

/*!
 * \brief create an integer problem
 */
integer_problem* integer_problem_factory::construct()
{
	string prob;
	kvparse::parameter_value(keywords::PROBLEM, prob, true);

	if(false) {
		// no valid integer problems defined yet
	} else {
		cerr << "invalid problem specified: " << prob << endl;
		exit(1);
		return 0;
	}
}

/*!
 * \brief create a generalized assignment problem instance
 */
gap_problem* gap_problem_factory::construct()
{
	string prob;
	kvparse::parameter_value(keywords::PROBLEM, prob, true);

	if(prob == "gap") {
		gap_problem* p = new gap_problem;
		p->initialize();
		return p;
	} else {
		cerr << "invalid problem specified: " << prob << endl;
		exit(1);
		return 0;
	}
}

/*!
 * \brief create a generalized sailor assignment problem instance
 */
gsap_problem* gsap_problem_factory::construct()
{
	string prob;
	kvparse::parameter_value(keywords::PROBLEM, prob, true);

	if(prob == "gsap") {
		gsap_problem* p = new gsap_problem;
		p->initialize();
		return p;
	} else {
		cerr << "invalid problem specified: " << prob << endl;
		exit(1);
		return 0;
	}
}

/*!
 * \brief create a binary problem
 */
bit_string_problem* bit_string_problem_factory::construct()
{
	string prob;
	kvparse::parameter_value(keywords::PROBLEM, prob, true);

	if(prob == "onemax") {
		bit_string_problem* p = new onemax_problem;
		p->initialize();
		return p;
	} else if(prob == "lotz") {
		bit_string_problem* p = new lotz_problem;
		p->initialize();
		return p;
	} else if(prob=="knapsack") {
		knapsack_problem* p=new knapsack_problem;
		p->initialize();
		return p;
	} else {
		cerr << "invalid problem specified: " << prob << endl;
		exit(1);
		return 0;
	}
}

/*!
 * \brief create a numeric problem
 */
numeric_problem* numeric_problem_factory::construct()
{
	string prob;
	kvparse::parameter_value(keywords::PROBLEM, prob, true);

	if(prob=="f1") {
		numeric_problem* p = new f1_problem;
		p->initialize();
		return p;
	}
	if(prob=="f2") {
		numeric_problem* p = new f2_problem;
		p->initialize();
		return p;
	}
	if(prob=="f3") {
		numeric_problem* p = new f3_problem;
		p->initialize();
		return p;
	}
	if(prob=="f4") {
		numeric_problem* p = new f4_problem;
		p->initialize();
		return p;
	}
	if(prob=="f5") {
		numeric_problem* p = new f5_problem;
		p->initialize();
		return p;
	} else if(prob=="kur") {
		numeric_problem* p = new kur_problem;
		return p;
	} else if(prob == "sch") {
		numeric_problem* p = new sch_problem;
		return p;
	} else if(prob == "zdt1") {
		numeric_problem* p = new zdt1_problem;
		return p;
	} else if(prob == "zdt2") {
		numeric_problem* p = new zdt2_problem;
		return p;
	} else if(prob == "zdt3") {
		numeric_problem* p = new zdt3_problem;
		return p;
	} else if(prob == "zdt4") {
		numeric_problem* p = new zdt4_problem;
		return p;
	} else if(prob == "zdt6") {
		numeric_problem* p = new zdt6_problem;
		return p;
	} else if(prob == "dtlz1") {
		numeric_problem* p = new dtlz1_problem;
		p->initialize();
		return p;
	} else if(prob == "rana") {
		return new rana_problem;
	} else {
		cerr << "invalid problem specified; " << prob << endl;
		exit(1);
		return 0;
	}
}

