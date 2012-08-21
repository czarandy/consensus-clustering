#pragma once

// Best of K algorithm using Adjusted Rand metric

#include "CCHeuristic.h"
#include <cmath>
using std::exp;

class AdjBestOfK : public CCHeuristic
{
protected:
	const SetPartition& FindBestOfK(const SetPartitionVector& SPV);

public:
	virtual SetPartition Run(const SetPartitionVector& SPV);
};

class AdjBestOneElementMove : public AdjBestOfK
{
protected:

public:
	virtual SetPartition Run(const SetPartitionVector& SPV);
};
