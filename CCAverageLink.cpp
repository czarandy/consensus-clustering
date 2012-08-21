#include "CCAverageLink.h"
#include "SetPartition.h"
#include "Matrix.h"
#include <cmath>

SetPartition CCAverageLink::Run(const SetPartitionVector &SPV)
{
	unsigned int k = SPV.size();
	unsigned int n = SPV[0].GetN();
	SetPartition cluster(n);
	cluster.MakeLast();

  Matrix<double> mat(n + 1);
  mat.FillWithValue(0.0);
  for (int p = 0; p < k; ++p)
  {
    for (int i = 1; i <= n; ++i)
    {
      for (int j = i + 1; j <= n; ++j)
      {
        if (!SPV[p].CoClustered(i, j))
        {
          mat(i, j) += 1.0 / (double)k;
          mat(j, i) += 1.0 / (double)k;
        }
      }
    }
  }

  while (true)
  {
    // find new min distance
    double minVal = mat(1, 2);
    int minI = 1;
    int minJ = 2;
    for (int i = 1; i <= n; ++i)
      for (int j = i + 1; j <= n; ++j)
      {
        if (cluster.CoClustered(i, j))
          continue;
        
        if (mat(i, j) < minVal)
        {
          minVal = mat(i, j);
          minI = i;
          minJ = j;
        }
      }

    if (minVal >= 0.5) // too far
      break;

    // recalculate distances
    double iSize = cluster.BlockSize(cluster.WhichBlock(minI));
    double jSize = cluster.BlockSize(cluster.WhichBlock(minJ));
    for (int k = 1; k <= n; ++k)
    {
      if (k == minI || k == minJ)
        continue;

      mat(minI, k) = mat(k, minI) = mat(minJ, k) =
        mat(k, minJ) = ((mat(minI, k) *
            iSize) + (mat(minJ, k) * jSize)) / (iSize + jSize);
    }

    // merge
    cluster.MergeBlocks(cluster.WhichBlock(minI),
        cluster.WhichBlock(minJ));

    // quit if all together
    if (cluster.BlockSize(1) == cluster.GetN())
      break;
  }

  return cluster;
}
