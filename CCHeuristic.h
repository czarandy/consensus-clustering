#pragma once

// Abstract base class for all consensus clustering heuristics/algorithms

#include "SetPartition.h"
#include "SetPartitionVector.h"

// Consensus Clustering heuristic
class CCHeuristic
{
public:
	virtual SetPartition Run(const SetPartitionVector& SPV)=0;
};
