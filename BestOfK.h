#pragma once

// Best of K algorithm using Rand distance metric

#include "CCHeuristic.h"
#include "MajorityRule.h"
#include "Matrix.h"
#include <cmath>
using std::exp;

class BestOfK : public CCHeuristic
{
protected:
	const SetPartition& FindBestOfK(const SetPartitionVector& SPV);

public:
	virtual SetPartition Run(const SetPartitionVector& SPV);
};

class BestOneElementMove : public BestOfK
{
protected:
	void GenerateK(const SetPartitionVector& SPV, Matrix<int>& K);
	void GenerateM(const SetPartitionVector& SPV, const SetPartition& sp, const
      Matrix<int>& K, Matrix<int>& M);

public:
	virtual SetPartition Run(const SetPartitionVector& SPV, bool AL = false);
  virtual SetPartition Run(const SetPartitionVector& SPV, SetPartition sp);
};

class SimulatedAnnealingOEM : public BestOneElementMove
{
protected:
	const double Alpha;
	const int IterLength;
	const double FinalTemp;
	const int SAConst;

	double ReduceTemp(double t)
	{
		return t * Alpha;
	}
	double InitTemp()
	{
		return 1.0;
	}
	int DeltaCost(const SetPartition& P, const Matrix<int>& K, int e, int blockfrom, int blockto);

public:
	SimulatedAnnealingOEM() : Alpha(0.999), IterLength(250), FinalTemp(0.001),
                            SAConst(1000) {}
	SimulatedAnnealingOEM(double _Alpha, int _IterLength, double _FinalTemp, int _SAConst) : Alpha(_Alpha), IterLength(_IterLength), FinalTemp(_FinalTemp), SAConst(_SAConst) {}
	virtual SetPartition Run(const SetPartitionVector& SPV, bool AL = false);
  virtual SetPartition Run(const SetPartitionVector& SPV, SetPartition sp);
};
