#pragma once

#include "CCHeuristic.h"
#include "SetPartitionVector.h"

class CCAverageLink : public CCHeuristic
{
private:
  double cutoff;
public:
  CCAverageLink(double _cutoff) : cutoff(_cutoff) {}
	virtual SetPartition Run(const SetPartitionVector& SPV);
};
