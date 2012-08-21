#include "MajorityRule.h"
#include "Matrix.h"
#include "SetPartition.h"

SetPartition MajorityRule::Run(const SetPartitionVector &SPV)
{
	unsigned int k = SPV.size();
	unsigned int n = SPV[0].GetN();
	SetPartition majority(n);
	majority.MakeLast();

  Matrix<double> mat(n + 1, n + 1);
  mat.FillWithValue(0.0);

	for (int p = 0; p < k; ++p)
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

	for (int i = 1; i <= n; ++i)
		{
			for (int j = i + 1; j <= n; ++j)
			{
				if (mat(i, j) > cutoff && !majority.CoClustered(i, j))
				{
						majority.MergeBlocks(majority.WhichBlock(i), majority.WhichBlock(j));
				}
			}
		}

	return majority;
}
