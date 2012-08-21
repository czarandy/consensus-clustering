#pragma once
#include "CCHeuristic.h"
#include "SetPartitionVector.h"

// Heuristic using majority rule to create a consensus clustering


class MajorityRule : public CCHeuristic
{
  private:
    double cutoff;
public:
    MajorityRule(double _cutoff) : cutoff(_cutoff) {}
	virtual SetPartition Run(const SetPartitionVector& SPV);
};
