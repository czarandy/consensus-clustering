#include "CCOptimal.h"

SetPartition CCOptimal::Run(const SetPartitionVector& SPV)
{
	if (SPV.size() < 1)
		return SetPartition::First(3);

	int n = SPV[0].GetN();
	int min = n * n * SPV.size(); // large enough
	SetPartition minSP(n);
	SetPartition sp = SetPartition::First(n);
	do
	{
		int dist = SPV.SumOfRandDistance(sp);
		if (dist < min)
		{
			min = dist;
			minSP = sp;
		}
		++sp;
	} while (!sp.IsFirst());

	return minSP;
}
