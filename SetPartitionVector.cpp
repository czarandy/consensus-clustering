#include "SetPartitionVector.h"
#include "Utility.h"

void SetPartitionVector::FillRandom(unsigned int n, unsigned int k)
{
	this->clear();
	while (k--)
	{
		this->push_back(SetPartition(n));
		this->back().Randomize();
	}
}

void SetPartitionVector::FillNoisy(unsigned int n, unsigned int k, unsigned int
    noise)
{
  this->clear();
  SetPartition sp(n);
  sp.Randomize();

  unsigned int moves = (n * noise) / 100;

  while (k--)
  {
    this->push_back(SetPartition(sp));

    for (unsigned int i = 0; i < moves; ++i)
      this->back().RandomOneElementMove();
  }
}

void SetPartitionVector::RandomizeAll()
{
	for (int i = 0; i < this->size(); ++i)
		(*this)[i].Randomize(); }

unsigned int SetPartitionVector::SumOfRandDistance(const SetPartition& sp) const
{
	unsigned int sum = 0;
	for (vector<SetPartition>::const_iterator i = this->begin(); i != this->end(); ++i)
		sum += i->RandDistanceTo(sp);
	return sum;
}

double SetPartitionVector::SumOfNormalizedRandDistance(const SetPartition& sp) const
{
	double sum = 0;
	for (vector<SetPartition>::const_iterator i = this->begin(); i != this->end(); ++i)
		sum += i->NormalizedRandDistanceTo(sp);
	return sum;
}

double SetPartitionVector::SumOfNormalRandDistance(const SetPartition& sp) const
{
	double sum = 0;
	for (vector<SetPartition>::const_iterator i = this->begin(); i != this->end(); ++i)
		sum += i->NormalRandDistanceTo(sp);
	return sum;
}

double SetPartitionVector::SumOfAdjustedRandDistance(const SetPartition& sp) const
{
	double sum = 0;
	for (vector<SetPartition>::const_iterator i = this->begin(); i != this->end(); ++i)
		sum += i->AdjustedRandDistanceTo(sp);
	return sum;
}

unsigned int SetPartitionVector::MinRandDistance(const SetPartition& sp) const
{
	unsigned int min = ChooseTwo(sp.elements.size());
	for (vector<SetPartition>::const_iterator i = this->begin(); i != this->end(); ++i)
	{
		unsigned int d = i->RandDistanceTo(sp);
		if (d < min)
			min = d;
	}
	return min;
}

unsigned int SetPartitionVector::MaxRandDistance(const SetPartition& sp) const
{
	unsigned int max = 0;
	for (vector<SetPartition>::const_iterator i = this->begin(); i != this->end(); ++i)
	{
		unsigned int d = i->RandDistanceTo(sp);
		d = i->RandDistanceTo(sp);
		if (d > max)
			max = d;
	}
	return max;
}

double SetPartitionVector::MinAdjustedRandDistance(const SetPartition& sp) const
{
	double min = ChooseTwo(sp.elements.size());
	for (vector<SetPartition>::const_iterator i = this->begin(); i != this->end(); ++i)
	{
		double d = i->AdjustedRandDistanceTo(sp);
		if (d < min)
			min = d;
	}
	return min;
}

double SetPartitionVector::MaxAdjustedRandDistance(const SetPartition& sp) const
{
	double max = 0;
	for (vector<SetPartition>::const_iterator i = this->begin(); i != this->end(); ++i)
	{
		double d = i->AdjustedRandDistanceTo(sp);
		d = i->AdjustedRandDistanceTo(sp);
		if (d > max)
			max = d;
	}
	return max;
}

// Fill with ALL set partitions of size n (note: this is a lot)
void SetPartitionVector::FillAll(unsigned int n)
{
	if (n > MAX_N || n == 0)
		return;

	this->resize(BellNumber<unsigned int>(n), n);
	iota(this->begin(), this->end(), SetPartition::First(n));
}

// Fill will ALL partitions reachable by OEMs (note: this can be a lot, up to roughly n^2)
void SetPartitionVector::FillAllOneElementMoves(const SetPartition& sp)
{
	this->clear();

	for (int element = 1; element < sp.elements.size(); ++element)
	{
		for (int block = 1; block <= sp.blocks.size(); ++block)
		{
			if (sp.elements[element] == block)
				continue; // No move will be done

			if (sp.blocks[sp.elements[element]].size() == 1 && block == sp.blocks.size())
				continue; // No move will be done

			this->push_back(sp);
			this->back().OneElementMove(element, block);
		}
	}
}

// Returns true if a and b are clustered in at least half of the set partitions
bool SetPartitionVector::MajorityClustering(int a, int b) const
{
	int count = 0;
	for (int i = 0; i < this->size(); ++i)
	{
		if ((*this)[i].CoClustered(a, b))
			++count;
	}

	return count >= (this->size() / 2);
}

double SetPartitionVector::SumOfDistanceMatrix() const
{
  double D = 0.0;
  /*
  for (int i = 0; i < size(); ++i)
  {
    for (int j = 0; j < size(); ++j)
    {
      D += (*this)[i].AdjustedRandDistanceTo((*this)[j]);
    }
  }
  */
  // Sample O(k) values
  for (int i = 0; i < size(); ++i)
  {
    int j = (*this)[0].mt.RandomInt(1, size()) - 1;
    D += (*this)[i].AdjustedRandDistanceTo((*this)[j]);
  }
  return D;
}

double SetPartitionVector::SubsetSumOfDistanceMatrix(const vector<int>& subset) const
{
  double D = 0.0;
  /*
  for (int i = 0; i < size(); ++i)
  {
    for (int j = 0; j < size(); ++j)
    {
      D += (*this)[i].AdjustedRandDistanceTo((*this)[j]);
    }
  }
  */
  // Sample O(k) values
  for (int i = 0; i < subset.size(); ++i)
  {
    int j = subset[(*this)[0].mt.RandomInt(0, subset.size() - 1)] - 1;
    D += (*this)[subset[i] - 1].AdjustedRandDistanceTo((*this)[j]);
  }
  return D;
}

void SetPartitionVector::LoadFile(const char* fileName)
{
  int N, k;
  ifstream fin(fileName);
  string s;

  fin >> N >> k;
  getline(fin, s);
  this->clear();
  while (k--)
  {
    getline(fin, s);
    this->push_back(SetPartition(N));
    this->back().LoadData(s);
  }
}
void SetPartitionVector::LoadFile(const char* fileName, vector<string>&
    elLabels)
{
  int N, k;
  ifstream fin(fileName);
  string s, s2;

  fin >> N >> k;
  getline(fin, s);

  for (int i = 0; i < N; ++i)
  {
    getline(fin, s);
    elLabels.push_back(s);
  }
  
  this->clear();
  while (k--)
  {
    getline(fin, s2);
    getline(fin, s);
    this->push_back(SetPartition(N, s2));
    this->back().LoadData(s);
  }
}

double SetPartitionVector::EmpiricalNoiseProb(double noise) const
{
  // for each element of the vector, generate m noisy samples
  // then get the sum of distance matrix
  int numMoves = (int)(noise * (*this)[0].GetN());
  vector<double> values;

  const int samples = 100;

  for (int i = 0; i < samples; ++i)
  {
    SetPartitionVector SPV;
    int which = SPV[0].mt.RandomInt(0, size() - 1);
    for (int j = 0; j < size(); ++j)
    {
      SPV.push_back((*this)[which]);
      //SPV.back().RandomOneElementMoves(numMoves);
      for (int k = 0; k < numMoves; ++k)
        SPV.back().RandomOneElementMove();
    }
    values.push_back(SPV.SumOfDistanceMatrix());
  }

  sort(values.begin(), values.end());

  double myDist = SumOfDistanceMatrix();
  
  int i;
  for (i = 0; i < values.size(); ++i)
  {
    if (values[i] > myDist) break;
  }
  return ((double) i / (double) values.size());
}

double SetPartitionVector::SubsetEmpiricalNoiseProb(const vector<int>& subset, double noise) const
{
  // for each element of the vector, generate m noisy samples
  // then get the sum of distance matrix
  int numMoves = (int)(noise * (*this)[0].GetN());
  vector<double> values;

  const int samples = 100;

  for (int i = 0; i < samples; ++i)
  {
    SetPartitionVector SPV;
    int which = subset[SPV[0].mt.RandomInt(0, subset.size() - 1)] - 1;
    for (int j = 0; j < subset.size(); ++j)
    {
      SPV.push_back((*this)[which]);
      //SPV.back().RandomOneElementMoves(numMoves);
      for (int k = 0; k < numMoves; ++k)
        SPV.back().RandomOneElementMove();
    }
    values.push_back(SPV.SumOfDistanceMatrix());
  }

  sort(values.begin(), values.end());

  double myDist = SubsetSumOfDistanceMatrix(subset);
  
  int i;
  for (i = 0; i < values.size(); ++i)
  {
    if (values[i] > myDist) break;
  }
  return ((double) i / (double) values.size());
}

double SetPartitionVector::CCLowerBound() const
{
  int n = (*this)[0].GetN();
  int k = size();

  int ret = 0;
  for (int i = 1; i <= n; ++i)
  {
    for (int j = i + 1; j <= n; ++j)
    {
      int count = 0;
      for (int m = 0; m < k; ++m)
        if ((*this)[m].CoClustered(i, j))
          ++count;

      ret += min(count, k - count);
    }
  }
  return ((double)ret / (double)k) /(double) ChooseTwo(n);
}
