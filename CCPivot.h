#pragma once
#include <set>
#include "CCHeuristic.h"
#include "SetPartitionVector.h"


class CCPivot : public CCHeuristic
{
  private:
    double cutoff;
    
    int ChooseRandom(const set<int>& isClustered, const SetPartition& sp)
      const;

  
public:
    CCPivot(double _cutoff) : cutoff(_cutoff) {}
    SetPartition _Run(const SetPartitionVector& SPV);
	virtual SetPartition Run(const SetPartitionVector& SPV);
};
