/*!
 * \file mtrandom.h
 *
 * implementation of the mersenne twister.  initial incorporation
 * into the C++ framework done by Jean-Paul Watson.  Subsequently
 * modified to include additional features by Deon Garrett.
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _MTRANDOM_H_
#define _MTRANDOM_H_

#include <cmath>
#include <ctime>
#include <iostream>
#include <vector>

using namespace std;

/*!
 * \class mtrandom
 */
class mtrandom
{
private:
	static bool m_first_time;
	static unsigned long m_seed;
	static int m_iset;
	static double m_fac;
	static double m_r;
	static double m_v1;
	static double m_v2;
	static double m_gset;
	static long m_iy;
	static long m_ir[98];

public:
	// returns the random number seed - a static method
	// is provided so instantiation of the class is not
	// required.
	static unsigned long seed();

	// initialize the random number seed from the configuration file,
	// if it was specified - if unspecified, this method is a no-op.
	static bool initialize();
	static void seed_random(unsigned long n);

public:
	mtrandom();
	unsigned long get_seed();
	void   reseed();
	double random();
	int    random(int ub);
	int    random(int ub, int lb);
	double random(double ub, double lb);
	double bias_random(int range, double bias);
	double gaussian();
	double gaussian(double mean, double stddev);
	vector<int> permutation(int n);
	void   shuffle(vector<int>& v);
};

#endif
