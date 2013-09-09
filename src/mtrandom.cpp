/*!
 * \file mtrandom.cpp
 *
 * implementation of the mersenne twister.  initial incorporation
 * into the C++ framework done by Jean-Paul Watson.  Subsequently
 * modified to include additional features by Deon Garrett.
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "mtrandom.h"
#include "configuration.h"
#include "keywords.h"

using namespace std;

// allocate space for the class static's.
bool mtrandom::m_first_time=true;
unsigned long mtrandom::m_seed;
int mtrandom::m_iset;
double mtrandom::m_fac;
double mtrandom::m_r;
double mtrandom::m_v1;
double mtrandom::m_v2;
double mtrandom::m_gset;
long mtrandom::m_iy;
long mtrandom::m_ir[98];

/* Period Parameters */
#define N          624
#define M          397
#define MATRIX_A   0x9908b0df
#define UPPER_MASK 0x80000000
#define LOWER_MASK 0x7fffffff

/* Tempering parameters */
#define TEMPERING_MASK_B      0x9d2c5680
#define TEMPERING_MASK_C      0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned long mt[N];     /* the array for the state vector  */
static int           mti = N+1; /* mti==N+1 means mt[N] is not initialized */

/*!
 * \brief return the seed used for the generator
 */
unsigned long mtrandom::seed()
{
	// creating a dummy instance of the class forces the
	// random number seed to be initialized, if it wasn't already.
	mtrandom rn;

	(void) rn.get_seed(); // a dummy to get rid of the 'un-used variable' compiler warning

	return m_seed;
}

/*!
 * \brief initialize the seed for the generator
 */
bool mtrandom::initialize()
{
	unsigned int seed;
	if(configuration::keyword_exists(keywords::RANDOM_SEED)) {
		if(configuration::unsigned_integer_parameter(keywords::RANDOM_SEED,seed)) {
			cerr << "Failed to read value for keyword " << keywords::RANDOM_SEED << endl;
			return false;
		}
		m_seed = static_cast<unsigned long>(seed);
		m_first_time = false;
	} else {
		// try to read from /dev/random; if that fails, use the system time
		int fd=open("/dev/random",O_RDONLY|O_NONBLOCK);
		if((fd!=-1) && (read(fd,&m_seed,sizeof(unsigned long))==sizeof(unsigned long))) {
			m_seed=static_cast<unsigned long>(m_seed);
		} else {
			cerr << "failed to read random bytes from /dev/random...falling back to system clock" << endl;
			m_seed = static_cast<unsigned long>(time(NULL));
		}
		m_first_time=false;
	}

	return true;
}

/*!
 * \brief constructor
 *
 * if no instance has been created previously, seeds the generator using the
 * system clock
 */
mtrandom::mtrandom()
{
	if(m_first_time == true) {
		// try to read from /dev/random; if that fails, use the system time
		int fd=open("/dev/random",O_RDONLY|O_NONBLOCK);
		if((fd!=-1) && (read(fd,&m_seed,sizeof(unsigned long))==sizeof(unsigned long))) {
			m_seed=static_cast<unsigned long>(m_seed);
		} else {
			cerr << "failed to read random bytes from /dev/random...falling back to system clock" << endl;
			m_seed = static_cast<unsigned long>(time(NULL));
		}

		m_first_time = false;
		//m_seed = static_cast<unsigned long>(time(NULL));
		//m_seed = ((m_seed>0) ? -m_seed : m_seed);
		seed_random(m_seed);
		m_iset=0;
	}
}

/*!
 * \brief seed the generator with a user defined seed
 */
void mtrandom::seed_random(unsigned long n)
{
	int i;
	unsigned long seedValue = static_cast<unsigned long>(n);
	for(i=0; i<N; i++) {
		mt[i] = seedValue & 0xffff0000;
		seedValue = 69069 * seedValue + 1;
		mt[i] |= (seedValue & 0xffff0000) >> 16;
		seedValue = 69069 * seedValue + 1;
	}
	mti = N;
}

/*!
 * \brief forces the next invocation to reseed from the clock
 */
void mtrandom::reseed()
{
	m_first_time = true;
}

/*!
 * \brief return the seed for the generator
 */
unsigned long mtrandom::get_seed()
{
	return m_seed;
}

/*!
 * \brief return a random number in the range [0,1)
 */
double mtrandom::random()
{
	unsigned long y;
	static unsigned long mag01[2] = {0x0, MATRIX_A};
	/* mag01[x] = x * MATRIX_A  for x=0,1 */

	if(mti >= N) { /* generate N words at one time */
		int kk;

		if(mti == N+1) {      /* if sgenrand() has not been called, */
			seed_random(4357); /* a default initial seed is used   */
		}

		for(kk=0; kk<N-M; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
			mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		for(; kk<N-1; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
			mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

		mti = 0;
	}

	y = mt[mti++];
	y ^= TEMPERING_SHIFT_U(y);
	y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
	y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
	y ^= TEMPERING_SHIFT_L(y);

	return ((double) y * 2.3283064365386963e-10); /* reals: [0,1)-interval */
	/* return y; */ /* for integer generation */
}

/*!
 * \brief return a random integer in the range [0,ub)
 */
int mtrandom::random(int ub)
{
	return random(0, ub);
}

/*!
 * \brief return a random integer in the range [lb,ub)
 */
int mtrandom::random(int lb, int ub)
{
	assert(ub>=lb);
	return (int(floor(this->random() * double(ub-lb))) + lb);
}

/*!
 * \brief return a random double in the range [lb,ub)
 */
double mtrandom::random(double lb, double ub)
{
	assert(ub>=lb);
	return (double(this->random() * (ub-lb)) + lb);
}

/*!
 * \brief return a biased random number in the range [0,range)
 */
double mtrandom::bias_random(int range, double bias)
{
	if(bias > 1.0)
		return double(double(range) *
		              (bias - sqrt(bias * bias - 4.0 * (bias-1) * random())) /
		              2.0/(bias-1.0));
	else {
		return(double (range) * random());
	}
}

/*!
 * \brief return a normally distributed number with mu=0, sigma=1
 */
double mtrandom::gaussian()
{
	if(m_iset == 0) {
		do {
			m_v1 = 2.0 * this->random() - 1.0;
			m_v2 = 2.0 * this->random() - 1.0;
			m_r  = m_v1 * m_v1 + m_v2 * m_v2;
		} while(m_r >= 1.0);

		m_fac  = sqrt(-2.0 * log(m_r) / m_r);
		m_gset = m_v1 * m_fac;
		m_iset = 1;
		return(m_v2 * m_fac);
	} else {
		m_iset = 0;
		return(m_gset);
	}
}

/*!
 * \brief return a normally distributed number with specified mu, sigma
 */
double mtrandom::gaussian(double mean, double stddev)
{
	if(m_iset == 0) {
		do {
			m_v1 = 2.0 * this->random() - 1.0;
			m_v2 = 2.0 * this->random() - 1.0;
			m_r  = m_v1 * m_v1 + m_v2 * m_v2;
		} while(m_r >= 1.0);

		m_fac  = sqrt(-2.0 * log(m_r) / m_r);
		m_gset = m_v1 * m_fac;
		m_iset = 1;
		return((m_v2*m_fac) * stddev + mean);
	} else {
		m_iset = 0;
		return(m_gset * stddev + mean);
	}
}

/*!
 * \brief return a random permutation of the numbers [0,n)
 */
vector<int> mtrandom::permutation(int n)
{
	vector<int> temp(n);
	vector<int> result(n);

	for(int i=0; i<n; i++) {
		temp[i] = i;
	}

	int ub = n - 1;
	for(int i=0; i<n; i++) {
		int next   = random(0,ub);
		result[i]  = temp[next];
		temp[next] = temp[ub];
		ub--;
	}

	return result;
}

/*!
 * \brief randomly permute the elements of the given vector
 */
void mtrandom::shuffle(vector<int>& v)
{
	int ub = (int)v.size()-1;

	for(unsigned int i=0; i<v.size(); i++) {
		int next = random(0,ub--);
		swap(v[i],v[next]);
	}
}

