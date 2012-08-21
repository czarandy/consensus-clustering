// Contains misc. utility functions/classes
#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <cmath>
#include <map>
#include <utility>
#include <numeric>
#include <functional>
using namespace std;

#include "SetPartition.h"

long double MaxBlockProbability(unsigned int n, unsigned int x);
long double Lambert_W_Initial(long double x);
long double Lambert_W(long double x, long double precision);
long double erfc(long double x);
long double NormalCDF(long double x);
long double BlockDistributionProbability(map<unsigned int, unsigned int>& m, unsigned int n);
long double mean(vector<long double>& v);
long double variance(vector<long double>& v);
long double skewness(vector<long double>& v);

// Compares two elements as equal if they have the same block in two set partitions
class InducedElementComparer
{
	const SetPartition& a;
	const SetPartition& b;

public:
	InducedElementComparer(const SetPartition& A, const SetPartition& B) : a(A), b(B) {}
	bool operator()(int x, int y)
	{
		if (a.elements[x] == a.elements[y])
			return b.elements[x] < b.elements[y];
		return a.elements[x] < a.elements[y];
	}
};

// This should be in <numeric>, but Visual C++ is dumb
template <class ForwardIterator, class T>
void iota(ForwardIterator first, ForwardIterator last, T value)
{
	T n = value;
	for (ForwardIterator i = first; i != last; ++i, ++value)
		*i = value;
}

unsigned int BlockSizeChooseTwo(unsigned int partialSum, const vector<int>& block);

template<class T> // e.g., unsigned int or unsigned long int
T ChooseTwo(T N)
{
	return N * (N - 1) / 2;
}

template<class Ret, class Par> // e.g., unsigned int or unsigned long int
Ret Factorial(Par N)
{
	if (N <= 1)
		return 1;

	Ret prod = N;
	while (--N)
		prod *= N;

	return prod;
}



template<class T> // e.g., unsigned int or unsigned long int
T BellNumber(unsigned int N)
{
	static long double E_INV = exp(-1.0); // 1/e

	long double sum = 0.0;
	for (long double m = 0; m < 2 * N; ++m)
		sum += pow(m, static_cast<long double>(N)) / Factorial<T, T>(m);
	
	return static_cast<T>(0.5 + sum * E_INV);
}

