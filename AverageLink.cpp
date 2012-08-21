#include "AverageLink.h"
#include "SetPartition.h"
#include "Matrix.h"
#include <cmath>

SetPartition AverageLink::Run(const SetPartitionVector &SPV)
{
	unsigned int k = SPV.size();
	unsigned int n = SPV[0].GetN();
	SetPartition cluster(k);
	cluster.MakeLast();

  double minVal = 1000.0;
  int minI, minJ;

  Matrix<double> mat(k);
	for (int i = 0; i < k; ++i)
    for (int j = 0; j < k; ++j)
    {
      if (i == j)
      {
        mat(i, j) = -1.0;
        continue;
      }
		  mat(i, j) = mat(j, i) = SPV[i].AdjustedRandDistanceTo(SPV[j]);
      
      if (mat(i, j) < minVal)
      {
        minVal = mat(i, j);
        minI = i;
        minJ = j;
      }
    }

  int iter = 0;
  while (true)
  {
    ++iter;

    // merge
    SetPartition newCluster(cluster);
    newCluster.MergeBlocks(newCluster.WhichBlock(minI + 1),
        newCluster.WhichBlock(minJ + 1));

    // should we quit?
    if (newCluster.BlockSize(newCluster.WhichBlock(minI + 1)) > 2)
    {
      double p;
      if ((p = SPV.SubsetEmpiricalNoiseProb(newCluster.GetBlock(
         newCluster.WhichBlock(minI + 1)), noiseCutoff)) > 0.95)
        { 
          return cluster;
        }
    }
    
    // or quit if everything is clustered
    if (newCluster.BlockSize(1) == newCluster.GetN())
      return newCluster;


    minVal = 1000.0;
    // recalculate the distances & find new min
    for (int j = 0; j < k; ++j)
    {
      if (j == minI || j == minJ)
        continue;
      double iSize = cluster.BlockSize(cluster.WhichBlock(minI + 1));
      double jSize = cluster.BlockSize(cluster.WhichBlock(minJ + 1));
      mat(minI, j) = mat(j, minI) = mat(minJ, j) =
        mat(j, minJ) = ((mat(minI, j) *
            iSize) + (mat(minJ, j) * jSize)) / (iSize + jSize);
    }

    // We are ok so do the real merge
    cluster.MergeBlocks(cluster.WhichBlock(minI + 1),
        cluster.WhichBlock(minJ + 1));

    for (int i = 0; i < k; ++i)
      for (int j = 0; j < k; ++j)
      {
        if (i == j || cluster.WhichBlock(i + 1) == cluster.WhichBlock(j + 1))
          continue;
        
        if (mat(i, j) < minVal)
        {
          minVal = mat(i, j);
          minI = i;
          minJ = j;
        }
      }

  }
}
