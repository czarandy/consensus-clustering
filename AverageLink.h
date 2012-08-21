#pragma once

// Clusters set partitions together using average link clustering
// with probabilistic distance metric and a "noise" cutoff

#include "CCHeuristic.h"
#include "SetPartitionVector.h"

class AverageLink : public CCHeuristic
{
  double noiseCutoff;

public:
  AverageLink() : noiseCutoff(0.25) {}
  AverageLink(double n) : noiseCutoff(n) {}
	virtual SetPartition Run(const SetPartitionVector& SPV);
};
