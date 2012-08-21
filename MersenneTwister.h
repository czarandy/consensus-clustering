#pragma once

#include <ctime>
using namespace std;

/* This is the Mersenne Twister PRNG class, adapted from the code at:
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html
 */
class MersenneTwister
{
	/* Period parameters */  
	static const int N = 624;
	static const int M = 397;
	static const unsigned long MATRIX_A = 0x9908b0dfUL; /* constant vector a */  
	static const unsigned long UPPER_MASK = 0x80000000UL; /* most significant w-r bits */
	static const unsigned long LOWER_MASK = 0x7fffffffUL; /* least significant r bits */

	unsigned long* mt; /* the array for the state vector  */
	int mti; /* mti == N + 1 means mt[N] is not initialized */

public:
	// For generating real numbers, this gives the type of interval, i.e [0, 1] vs. (0, 1) vs. (0, 1] vs. [0, 1)
	enum IntervalType { Open, Closed, RightOpen };

	MersenneTwister();
	MersenneTwister(unsigned long seed);
	MersenneTwister(unsigned long init_key[], int key_length);
	~MersenneTwister();

	// These give the random numbers
	unsigned long RandomUL();
	int RandomInt();
	int RandomInt(int a, int b); 
	long double RandomReal(IntervalType = RightOpen);
};

