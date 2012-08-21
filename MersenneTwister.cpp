#include "MersenneTwister.h"

MersenneTwister::MersenneTwister() : mti(N + 1), mt(new unsigned long[N])
{
	unsigned long seed = static_cast<unsigned long>(time(NULL));

	/* initializes mt[N] with a seed */

    mt[0]= seed & 0xffffffffUL;
    for (mti = 1; mti < N; ++mti)
	{
        mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

MersenneTwister::MersenneTwister(unsigned long seed) : mti(N + 1), mt(new unsigned long[N])
{
	/* initializes mt[N] with a seed */

    mt[0]= seed & 0xffffffffUL;
    for (mti = 1; mti < N; ++mti)
	{
        mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
MersenneTwister::MersenneTwister(unsigned long init_key[], int key_length)
{
    int i = 1, j = 0, k = (N > key_length ? N : key_length);

    MersenneTwister(19650218UL);

    for (; k; --k) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        ++i; 
		++j;
        if (i >= N)
		{ 
			mt[0] = mt[N - 1]; 
			i = 1; 
		}
        if (j >= key_length) 
			j = 0;
    }
    for (k = N - 1; k; --k) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        ++i;
        if (i >= N) 
		{ 
			mt[0] = mt[N - 1]; 
			i = 1; 
		}
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

// Frees up memory used by value array
MersenneTwister::~MersenneTwister()
{
	delete [] mt;
}

// generates a random number on [0,0xffffffff]-interval 
unsigned long MersenneTwister::RandomUL()
{
	unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            MersenneTwister(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
int MersenneTwister::RandomInt()
{
    return static_cast<int>(RandomUL() >> 1);
}

// Generates a real number on [0, 1] or [0, 1) or (0, 1)
long double MersenneTwister::RandomReal(IntervalType it)
{
	switch (it)
	{
	case Closed:
		return RandomUL() * (1.0 / 4294967295.0);
	case Open:
		return ((static_cast<double>(RandomUL())) + 0.5) * (1.0 / 4294967296.0); 
	case RightOpen:	 
	default: // Default to RightOpen
		return RandomUL() * (1.0 / 4294967296.0);
	}
}

// Gives a random integer in the range [a, b]
int MersenneTwister::RandomInt(int a, int b)
{
	return static_cast<int>((RandomReal(RightOpen) * static_cast<double>(b - a + 1)) + static_cast<double>(a));
}

