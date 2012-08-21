#pragma once

#include <vector>
#include <utility>
#include <fstream>
using namespace std;

#include "SetPartition.h"

class SetPartitionVector : public vector<SetPartition>
{
	static const unsigned int MAX_N = 14; // Max n for all partitions

public:
	void FillRandom(unsigned int n, unsigned int k);
	void FillAll(unsigned int n);
	void FillAllOneElementMoves(const SetPartition& sp);
  void FillNoisy(unsigned int n, unsigned int k, unsigned int noise);

	unsigned int SumOfRandDistance(const SetPartition& sp) const;
	double SumOfNormalizedRandDistance(const SetPartition& sp) const;
	double SumOfNormalRandDistance(const SetPartition& sp) const;
	double SumOfAdjustedRandDistance(const SetPartition& sp) const;
  double AverageSOD(const SetPartition& sp) const
  {
    return SumOfNormalizedRandDistance(sp) / size();
  }
	unsigned int MinRandDistance(const SetPartition& sp) const;
	unsigned int MaxRandDistance(const SetPartition& sp) const;
	double MinAdjustedRandDistance(const SetPartition& sp) const;
	double MaxAdjustedRandDistance(const SetPartition& sp) const;

	bool MajorityClustering(int a, int b) const;

  double SumOfDistanceMatrix() const;
  double SubsetSumOfDistanceMatrix(const vector<int>& subset) const;
  double EmpiricalNoiseProb(double noise = 0.25) const;
  double SubsetEmpiricalNoiseProb(const vector<int>& subset, double noise =
      0.25) const;

  double CCLowerBound() const;

	void RandomizeAll();

  void LoadFile(const char* fileName);
  // load file with labels
  void LoadFile(const char* fileName, vector<string>& elLabels);

};
