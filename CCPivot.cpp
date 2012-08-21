#include "CCPivot.h"
#include "Matrix.h"
#include "SetPartition.h"
#include "MersenneTwister.h"
#include <cstdlib>
#include <ctime>

int CCPivot::ChooseRandom(const set<int>& isClustered, const SetPartition&
    sp) const
{
  int n = sp.GetN();

  if (isClustered.size() == 0)
    return -1;

  double count = 1;
  for (set<int>::const_iterator itr = isClustered.begin(); itr !=
      isClustered.end(); ++itr)
  {
    if (sp.RandomReal() <= (count / (double)isClustered.size()))
      return *itr;

    ++count;
  }
}

SetPartition CCPivot::Run(const SetPartitionVector& SPV)
{
  srand((unsigned)time(NULL));
  SetPartition sp(_Run(SPV));
  int mindist = SPV.SumOfRandDistance(sp);

  for (int i = 0; i < 49; ++i)
  {
    SetPartition sp2(_Run(SPV));
    int d = SPV.SumOfRandDistance(sp2);
    if (d < mindist)
    {
      mindist = d;
      sp = sp2;
    }
  }
  return sp;
}

SetPartition CCPivot::_Run(const SetPartitionVector &SPV)
{
	unsigned int k = SPV.size();
	unsigned int n = SPV[0].GetN();
	SetPartition consensus(n);
	consensus.MakeLast();

  Matrix<double> mat(n + 1, n + 1);
  mat.FillWithValue(0.0);
	for (int p = 0; p < SPV.size(); ++p)
	{
		for (int i = 1; i <= n; ++i)
		{
			for (int j = 1; j <= n; ++j)
			{
				if (SPV[p].CoClustered(i, j))
					mat(i, j) += 1.0 / (double)k;
			}
		}
	}

  set<int> isClustered;
  for (int i = 1; i <= n; ++i)
    isClustered.insert(i);

  while (true)
  {
    int which = ChooseRandom(isClustered, SPV[0]);
    if (which == -1)
      break;
    isClustered.erase(which);

    // cluster everyone with which
    for (int i = 1; i <= n; ++i)
    {
      if (isClustered.find(i) == isClustered.end())
        continue;

      if (mat(which, i) > cutoff)
      {
        isClustered.erase(i);
				consensus.MergeBlocks(consensus.WhichBlock(i), consensus.WhichBlock(which));
      }
    }
  }

  return consensus;
}
