// Contains misc. utility functions/classes

#include "Utility.h"

unsigned int BlockSizeChooseTwo(unsigned int partialSum, const vector<int>& block)
{
	unsigned int N = block.size();
	return partialSum + ChooseTwo(N);
}

long double Lambert_W_Initial(long double x) {
	long double lx1;
	if (x <= 500.0)
	{
		lx1 = log(x + 1.0);
		return 0.665 * (1 + 0.0195 * lx1) * lx1 + 0.04;
	}
	return log(x - 4.0) - (1.0 - 1.0/log(x)) * log(log(x));
}

long double Lambert_W(long double x, long double precision = 1e-15)
{
	long double w = Lambert_W_Initial(x);
	int iters = 100;
	while (iters--)
	{
		long double wTimesExpW = w * exp(w);
		long double wPlusOneTimesExpW = (w + 1) * exp(w);
		if (precision > abs((x - wTimesExpW) / wPlusOneTimesExpW))
			break;
		w = w - (wTimesExpW - x) / (wPlusOneTimesExpW - (w + 2) * (wTimesExpW - x) / (2 * w + 2));
	}
	return w;
}

// Returns the probability the max block is of size x in a set partition of size n
long double MaxBlockProbability(unsigned int n, unsigned int x)
{
	// Some special cases
	if (x > n || x == 0) return 0.0;
	if (x == n || x == 1) return 1.0 / BellNumber<double>(n);

	static long double E = exp(1.0);
	static long double PI = acos(-1.0);

	long double R = Lambert_W(n);
	long double lsqrtr = log(sqrt(R));
	long double delta = E * R - log(sqrt(2 * PI * E * R)) - log(E - 1);
	long double deltaFract = delta - floor(delta);

	long double X = x - floor(delta) + 0.5;
	
	return exp(-1.0 * exp(-1.0 * X + deltaFract)) - exp(-1.0 * exp(-1.0 * (X - 1) + deltaFract));
}

// Returns complimentary error function
long double erfc(long double x)
{
	long double t,z,ans;
	z=fabs(x);
	t=1.0/(1.0+0.5*z);
	ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
	t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
	t*(-0.82215223+t*0.17087277)))))))));
	return x >= 0.0 ? ans : 2.0-ans;
}

// Returns standard normal cdf
long double NormalCDF(long double x)
{
	double erf = 1.0 - erfc(x / sqrt(2.0));
	return (1.0 + erf) / 2.0;
}

long double BlockDistributionProbability(map<unsigned int, unsigned int>& m, unsigned int n)
{
	long double p = Factorial<long double, unsigned int>(n);
	long double denom = BellNumber<long double>(n);

	for (map<unsigned int, unsigned int>::iterator i = m.begin(); i != m.end(); ++i)
	{
		denom *= pow(Factorial<long double, unsigned int>(i->first), static_cast<long double>(i->second)) * Factorial<long double, unsigned int>(i->second);
	}
	
	return p / denom;
}

long double mean(vector<long double>& v)
{
	long double sum = accumulate(v.begin(), v.end(), 0.0);
	return sum / v.size();
}

long double variance(vector<long double>& v)
{
	long double sum = 0.0, mu = mean(v);
	for (int i = 0; i < v.size(); ++i)
		sum += pow(static_cast<long double>(v[i] - mu), 2.0L);

	return sum / v.size();
}

long double skewness(vector<long double>& v)
{
	long double sum = 0.0, mu = mean(v);
	for (int i = 0; i < v.size(); ++i)
		sum += pow(v[i] - mu, 3.0L);

	return sum / v.size();
}
