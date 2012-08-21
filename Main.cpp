#include <iostream>
#include <ctime>
#include <string>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include <sstream>
using namespace std;

#include "MersenneTwister.h"
#include "RefinedClustering.h"
#include "SetPartition.h"
#include "SetPartitionVector.h"
#include "Utility.h"
#include "BestOfK.h"
#include "AdjBestOfK.h"
#include "MajorityRule.h"
#include "AverageLink.h"
#include "Matrix.h"
#include "CCPivot.h"
#include "CCAverageLink.h"

int main(int argc, char** argv)
{
  // Example of refined consensus clustering
  // Sample input is "mushroom" data
  // Using a cutoff of 0.25 is typical--we consider it the maximum point where
  // random takes over
  RefinedClustering::Run("mushroom.in", 0.25);
  
  // Example of generating two random data sets and finding a consensus
  // clustering for each of them
  SetPartitionVector RandomData, NoisyData;

  // the arguments are n and k, respectively
  RandomData.FillRandom(100, 50);
  NoisyData.FillNoisy(100, 50, 10); // last argument is noise (in %)
  
  // We'll use two CC algorithms, majority rule and cc pivot
  // the arguments are the cutoffs
  CCPivot Pivot(0.5);
  MajorityRule MR(0.5);

  // Let's see which algorithm runs better on each data set
  double dPivotRandom = RandomData.AverageSOD(Pivot.Run(RandomData));
  double dPivotNoisy = NoisyData.AverageSOD(Pivot.Run(NoisyData));
  double dMajRandom = RandomData.AverageSOD(MR.Run(RandomData));
  double dMajNoisy = NoisyData.AverageSOD(MR.Run(NoisyData));

  if (dPivotRandom < dMajRandom)
    cout << "CCPivot wins at random data." << endl;
  else
    cout << "Majority Rule wins at random data." << endl;

  // Compare the best algorithm to the "lower bound"
  double minNoisy = min(dMajNoisy, dPivotNoisy);

  cout << "We got " << minNoisy << " whereas lower bound is " <<
    NoisyData.CCLowerBound() << endl;

  return 0;
}
