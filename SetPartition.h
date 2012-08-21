#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <string>
#include <ctime>
#include <map>
#include <fstream>
#include <utility>
#include <sstream>
#include <numeric>
#include <functional>
#include <string>
using namespace std;

#include "MersenneTwister.h"

class SetPartition
{
private:
	vector<int> elements;
	vector< vector<int> > blocks;
	static MersenneTwister mt;
  string label;

	int RandomEquivalenceRelation(); // Generates the random set partition
	void SetUpBlocks(unsigned int blockCount); // Sets up blocks

public:
	SetPartition(const SetPartition& sp) : elements(sp.elements),
                                         blocks(sp.blocks), label(sp.label) {}
	SetPartition(unsigned int n) : elements(n + 1) {}
  SetPartition(unsigned int n, const string& lbl) : elements(n + 1), label(lbl) {}
  SetPartition(unsigned int n, const char* lbl) : elements(n + 1), label(lbl) {}
	~SetPartition()
	{
	}

  void LoadData(string data);
  static SetPartition LoadClusterFile(const char* filename);

	unsigned int GetN() const { return elements.size() - 1; }

  const string& GetLabel() const { return label; }
  void SetLabel(const string& s) { label = s; }

	// Generate the partition
	void Randomize();
	void Induce(const SetPartition& a, const SetPartition& b);
	unsigned int RandDistanceTo(const SetPartition& sp) const;
	double NormalizedRandDistanceTo(const SetPartition& sp) const;
	double AdjustedRandDistanceTo(const SetPartition& sp) const;
	double NormalRandDistanceTo(const SetPartition& sp) const;
	void MakeFirst(); // Each element is seperate
	void MakeLast(); // All elements are together
	void Next(); // Next set partition

  void MergeBlocks(int a, int b);
	// These are faster than checking for equality
	bool IsFirst();
	bool IsLast();

	void PrintDistributionOfBlockSizes();

	int NumberOfBlocks() const { return blocks.size() - 1; }
	int BlockSize(int b) const { if (b < blocks.size()) return blocks[b].size(); return 0; }
	int BlockElement(int b, int i) const { return blocks[b][i]; }
	int WhichBlock(int element) const { return elements[element]; }

	int RandomBlockElement(int b) { return blocks[b][mt.RandomInt(0, blocks[b].size() - 1)]; }

	bool CoClustered(int a, int b) const;

  const vector<int>& GetBlock(int b) { return blocks[b]; }

  bool Verify() const;

	// Make a one element move
	bool OneElementMove(int element, int block);
	void RandomOneElementMove(); 
  void RandomOneElementMoves(int n);
	void RandomTransition(int& element, int& block);

  void RandomElementSwaps(int N);
  void RandomElementTripPerms(int N);

	unsigned int LargestBlockSize();
	long double ProbabilityOfBlockConfiguration(); // Returns the probability of this block configuration

	long double RandomReal() const { return mt.RandomReal(); }
  int RandomInt(int a, int b) const { return mt.RandomInt(a, b); }

	// Some static helper methods
	static SetPartition Random(unsigned int n)
	{
		SetPartition SP(n);
		SP.Randomize();
		return SP;
	}
	static SetPartition Induced(const SetPartition& a, const SetPartition& b)
	{
		SetPartition SP(a.elements.size() - 1);
		SP.Induce(a, b);
		return SP;
	}
	static SetPartition First(unsigned int n)
	{
		SetPartition SP(n);
		SP.MakeFirst();
		return SP;
	}
	static SetPartition Last(unsigned int n)
	{
		SetPartition SP(n);
		SP.MakeLast();
		return SP;
	}

	bool operator==(const SetPartition& sp);	
	bool operator!=(const SetPartition& sp)
	{
		return !(*this == sp);
	}
	SetPartition& operator++() { Next(); return *this; }
	SetPartition& operator+=(int n) 
	{ 
		while (n--) Next(); 
		return *this;
	}
	SetPartition operator+(int n)
	{
		SetPartition SP(*this);
		SP += n;
		return SP;
	}

	friend ostream& operator<<(ostream&, const SetPartition&);
	friend class InducedElementComparer;
	friend class SetPartitionVector;
};

