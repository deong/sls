/*!
 * \file problems.h
 *
 * defines the objective functions to be optimized by the search
 * algorithms
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _PROBLEMS_H_
#define _PROBLEMS_H_

#include <vector>

using namespace std;

/*!
 * \class problem
 *
 * base class for all problems
 */
class problem
{
protected:
	// private constructor -- only friend can create
	// ensures proper initialization
	problem();

	// disable copying of functional class
	problem(const problem& that);
	problem& operator=(const problem& that);

public:
	virtual ~problem();

	// provide a default initialization method for problems that do not
	// require custom initialization
	virtual void initialize();

	// these methods must be provided by all subclasses
	virtual unsigned int dimensions() const = 0;
	virtual unsigned int objectives() const = 0;
};

/*!
 * \class numeric_problem
 *
 * problems defined by arithmetic expressions
 */
class numeric_problem : public problem
{
private:
	numeric_problem(const numeric_problem& that);
	numeric_problem& operator=(const numeric_problem& that);

protected:
	numeric_problem();

public:
	virtual ~numeric_problem();
	virtual pair<double,double> parameter_range(unsigned int index) const = 0;
	virtual bool evaluate(const vector<double>& p, vector<double>& fit) const = 0;
};

/*!
 * \class bit_string_problem
 *
 * problems operating directly on binary strings
 */
class bit_string_problem : public problem
{
private:
	bit_string_problem(const bit_string_problem& that);
	bit_string_problem& operator=(const bit_string_problem& that);

protected:
	bit_string_problem();

public:
	virtual ~bit_string_problem();
	virtual bool evaluate(const vector<int>& p, vector<int>& fit) const = 0;
};

/*!
 * \class permutation_problem
 */
class permutation_problem : public problem
{
public:
	typedef long FitnessType;

private:
	permutation_problem(const permutation_problem& that);
	permutation_problem& operator=(const permutation_problem& that);

protected:
	permutation_problem();

public:
	virtual ~permutation_problem();
	virtual bool evaluate(const vector<int>& p, vector<FitnessType>& fit) const = 0;
};

/*!
 * \class integer_problem
 */
class integer_problem : public problem
{
private:
	integer_problem(const integer_problem& that);
	integer_problem& operator=(const integer_problem& that);

protected:
	integer_problem();

public:
	virtual ~integer_problem();
	virtual void legal_values(unsigned int index, vector<int>& vals) const = 0;
	virtual bool evaluate(const vector<int>& p, vector<int>& fit) const = 0;
};

/*!
 * \class qap_problem
 */
class qap_problem : public permutation_problem
{
private:
	int*** mqap_flow;
	int**  mqap_dist;
	int    mqap_n;
	int    mqap_obj;
	FitnessType*** mqap_delta;

public:
	qap_problem();
	virtual ~qap_problem();

	// methods specific to the mqap fast update procedures
	virtual void compute_delta(vector<int>& p, int i, int j);
	virtual void compute_delta_part(vector<int>& p, int i, int j, int r, int s);
	virtual void update(vector<int>& p, int i, int j);
	virtual void get_delta(int i, int j, vector<FitnessType>& del);

	// general problem methods
	virtual void initialize();
	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual bool evaluate(const vector<int>& p, vector<FitnessType>& fit) const;
};

/*!
 * \class gap_problem
 */
class gap_problem : public integer_problem
{
public:
	unsigned int agents;
	unsigned int tasks;
	unsigned int mgap_obj;
	int*   capacity;
	int**  resources;
	int*** cost;

public:
	gap_problem();
	virtual ~gap_problem();

	virtual void initialize();
	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual void legal_values(unsigned int index, vector<int>& vals) const;
	virtual bool evaluate(const vector<int>& p, vector<int>& fit) const;
};

/*!
 * \class gsap_problem
 * \brief generalized sailor assignment problem
 *
 * Basic problem description is that we have N agents who must perform all
 * required work, partitioned into M<N task classes, over some period of time.
 * The time periods are broken down two levels deep: on the first level we
 * have a number of distinct days.  Then each day is further subdivided into
 * distinct non-overlapping time slots.  Each time slot requires a specified
 * set of tasks to be performed.
 *
 * Each agent is associated a specific number of work units required to perform
 * each task, and incurs a given cost for that task.  The goal is to assign
 * one agent to every task required in every time slot without causing any
 * agent to exceed a set capacity constraint and while minimizing the total
 * cost incurred.
 *
 * Note that there may be several different and conflicting notions of cost,
 * leading to a multiobjective formulation of the problem.  In addition, it
 * is desired that the number of sailors required to cover the assigned work
 * be minimized.  We may also wish to define other objectives to, for example,
 * reward solutions in which agents are assigned tasks in contiguous time
 * slots.
 */
class gsap_problem : public integer_problem
{
public:
	struct task_element {
		int day;
		int time_slot;
		int task;
		task_element(int d, int s, int t) : day(d), time_slot(s), task(t) {}
	};

public:
	gsap_problem();
	virtual ~gsap_problem();
	virtual void initialize();
	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual void legal_values(unsigned int index, vector<int>& vals) const;
	virtual bool evaluate(const vector<int>& p, vector<int>& fit) const;
	const vector<task_element>& get_elements() const;
	const vector<vector<int> >& get_agent_task_map() const;
	unsigned int agents() const;
	unsigned int tasks() const;
	const vector<unsigned int>& time_slot_boundaries() const;
	const vector<unsigned int>& capacities() const;
	const vector<vector<unsigned int> >& resources() const;

protected:
	unsigned int gsap_obj;
	unsigned int gsap_dim;
	unsigned int n_agents;      //!< total number of available agents
	unsigned int n_tasks;       //!< number of task classes
	unsigned int n_days;        //!< number of days to schedule
	unsigned int n_time_slots;  //!< number of time slots in each day
	vector<task_element> m_elements; //!< particular elements of the schedule
	vector<vector<int> > m_agents_for_task; //!< which agents can do each task
	vector<vector<vector<unsigned int> > > m_tasks;      //!< tasks required for each day/time slot
	vector<unsigned int> m_capacities;   //!< work units allowed for each agent
	vector<vector<unsigned int> > m_resources;   //!< work units required for each agent/task
	vector<vector<vector<unsigned int> > > m_costs;      //!< costs incurred by each agent/task
	vector<unsigned int> tsb;     //!< time slot boundaries (first index of each slot)
};

/*!
 * \class knapsack_problem
 * \brief multiobjective knapsack problem
 */
class knapsack_problem : public bit_string_problem
{
private:
	unsigned int ks_obj;
	unsigned int ks_dim;
	unsigned int ks_capacity;
	vector<unsigned int> ks_weights;
	vector<vector<unsigned int> > ks_values;

public:
	knapsack_problem();
	virtual ~knapsack_problem();
	virtual void initialize();
	virtual unsigned int objectives() const;
	virtual unsigned int dimensions() const;
	virtual bool evaluate(const vector<int>& p, vector<int>& fit) const;
};

/*!
 * \class onemax_problem
 */
class onemax_problem : public bit_string_problem
{
private:
	unsigned int onemax_dim;

public:
	onemax_problem();
	virtual ~onemax_problem();

	virtual void initialize();
	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual bool evaluate(const vector<int>& p, vector<int>& fit) const;
};

/*!
 * \class lotz_problem
 *
 * leading-ones-trailing-zeros
 */
class lotz_problem : public bit_string_problem
{
private:
	unsigned int lotz_dim;

public:
	lotz_problem();
	virtual ~lotz_problem();

	virtual void initialize();
	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual bool evaluate(const vector<int>& p, vector<int>& fit) const;
};

/*!
 * \class f1_problem
 */
class f1_problem : public numeric_problem
{
private:
	unsigned int f1_dim;

public:
	f1_problem();
	virtual ~f1_problem();

	virtual void initialize();
	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual pair<double,double> parameter_range(unsigned int index) const;
	virtual bool evaluate(const vector<double>& p, vector<double>& fit) const;
};

/*!
 * \class f2_problem
 */
class f2_problem : public numeric_problem
{
private:
	unsigned int f2_dim;

public:
	f2_problem();
	virtual ~f2_problem();

	virtual void initialize();
	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual pair<double,double> parameter_range(unsigned int index) const;
	virtual bool evaluate(const vector<double>& p, vector<double>& fit) const;
};

/*!
 * \class f3_problem
 */
class f3_problem : public numeric_problem
{
private:
	unsigned int f3_dim;

public:
	f3_problem();
	virtual ~f3_problem();

	virtual void initialize();
	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual pair<double,double> parameter_range(unsigned int index) const;
	virtual bool evaluate(const vector<double>& p, vector<double>& fit) const;
};

/*!
 * \class f4_problem
 */
class f4_problem : public numeric_problem
{
private:
	unsigned int f4_dim;

public:
	f4_problem();
	virtual ~f4_problem();

	virtual void initialize();
	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual pair<double,double> parameter_range(unsigned int index) const;
	virtual bool evaluate(const vector<double>& p, vector<double>& fit) const;
};

/*!
 * \class f5_problem
 */
class f5_problem : public numeric_problem
{
private:
	unsigned int f5_dim;

public:
	f5_problem();
	virtual ~f5_problem();

	virtual void initialize();
	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual pair<double,double> parameter_range(unsigned int index) const;
	virtual bool evaluate(const vector<double>& p, vector<double>& fit) const;
};

/*!
 * \class kur_problem
 */
class kur_problem : public numeric_problem
{
public:
	kur_problem();
	virtual ~kur_problem();

	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual pair<double,double> parameter_range(unsigned int index) const;
	virtual bool evaluate(const vector<double>& p, vector<double>& fit) const;
};

/*!
 * \class sch_problem
 */
class sch_problem : public numeric_problem
{
public:
	sch_problem();
	virtual ~sch_problem();

	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual pair<double,double> parameter_range(unsigned int index) const;
	virtual bool evaluate(const vector<double>& p, vector<double>& fit) const;
};

/*!
 * \class zdt1_problem
 */
class zdt1_problem : public numeric_problem
{
public:
	zdt1_problem();
	virtual ~zdt1_problem();

	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual pair<double,double> parameter_range(unsigned int index) const;
	virtual bool evaluate(const vector<double>& p, vector<double>& fit) const;
};

/*!
 * \class zdt2_problem
 */
class zdt2_problem : public numeric_problem
{
public:
	zdt2_problem();
	virtual ~zdt2_problem();

	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual pair<double,double> parameter_range(unsigned int index) const;
	virtual bool evaluate(const vector<double>& p, vector<double>& fit) const;
};

/*!
 * \class zdt3_problem
 */
class zdt3_problem : public numeric_problem
{
public:
	zdt3_problem();
	virtual ~zdt3_problem();

	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual pair<double,double> parameter_range(unsigned int index) const;
	virtual bool evaluate(const vector<double>& p, vector<double>& fit) const;
};

/*!
 * \class zdt4_problem
 */
class zdt4_problem : public numeric_problem
{
public:
	zdt4_problem();
	virtual ~zdt4_problem();

	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual pair<double,double> parameter_range(unsigned int index) const;
	virtual bool evaluate(const vector<double>& p, vector<double>& fit) const;
};

/*!
 * \class zdt6_problem
 */
class zdt6_problem : public numeric_problem
{
public:
	zdt6_problem();
	virtual ~zdt6_problem();

	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual pair<double,double> parameter_range(unsigned int index) const;
	virtual bool evaluate(const vector<double>& p, vector<double>& fit) const;
};

/*!
 * \class dtlz1_problem
 */
class dtlz1_problem : public numeric_problem
{
private:
	unsigned int dtlz1_dim;
	unsigned int dtlz1_obj;

public:
	dtlz1_problem();
	virtual ~dtlz1_problem();

	virtual void initialize();
	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual pair<double,double> parameter_range(unsigned int index) const;
	virtual bool evaluate(const vector<double>& p, vector<double>& fit) const;
};

/**
 * @class rana_problem
 * @brief
 */
class rana_problem : public numeric_problem
{
public:
	rana_problem();
	virtual ~rana_problem();

	virtual unsigned int dimensions() const;
	virtual unsigned int objectives() const;
	virtual pair<double,double> parameter_range(unsigned int index) const;
	virtual bool evaluate(const vector<double>& p, vector<double>& fit) const;
};

/*!
 * \class numeric_problem_factory
 */
class numeric_problem_factory
{
public:
	static numeric_problem* construct();
};

/*!
 * \class bit_string_problem_factory
 */
class bit_string_problem_factory
{
public:
	static bit_string_problem* construct();
};

/*!
 * \class permutation_problem_factory
 */
class permutation_problem_factory
{
public:
	static permutation_problem* construct();
};

/*!
 * \class integer_problem_factory
 */
class integer_problem_factory
{
public:
	static integer_problem* construct();
};

/*!
 * \class gap_problem_factory
 */
class gap_problem_factory
{
public:
	static gap_problem* construct();
};

/*!
 * \class gap_problem_factory
 */
class gsap_problem_factory
{
public:
	static gsap_problem* construct();
};

#endif
