#pragma once
#include "CCHeuristic.h"

// Simple optimal algorithm that tries all possible set partitions

class CCOptimal : public CCHeuristic
{
public:
	virtual SetPartition Run(const SetPartitionVector& SPV);
};
