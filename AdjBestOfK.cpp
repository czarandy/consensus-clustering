#include "AdjBestOfK.h"
#include <cassert>

const SetPartition& AdjBestOfK::FindBestOfK(const SetPartitionVector& SPV)
{
	unsigned int min = SPV.SumOfAdjustedRandDistance(SPV[0]);
	unsigned int minSP = 0;

	for (int i = 1; i < SPV.size(); ++i)
	{
		unsigned int d = SPV.SumOfAdjustedRandDistance(SPV[i]);
		if (d < min)
		{
			min = d;
			minSP = i;
		}
	}
	return SPV[minSP];
}

SetPartition AdjBestOfK::Run(const SetPartitionVector& SPV)
{
	return FindBestOfK(SPV);
}


SetPartition AdjBestOneElementMove::Run(const SetPartitionVector& SPV)
{
	SetPartition candidate(FindBestOfK(SPV));
	

	unsigned int n = candidate.GetN();

	// Main loop
	while (true)
	{
		// Find best move
		double S = SPV.SumOfAdjustedRandDistance(candidate);
		double bestDeltaS = 0.0;
		int bestElement = 1, bestBlock = 1;
		
		for (int el = 1; el <= n; ++el)
		{
			for (int block = 1; block <= candidate.NumberOfBlocks(); ++block)
			{
				int blockFrom = candidate.WhichBlock(el);
				if (blockFrom == block)
					continue;

				SetPartition nextMove(candidate);

				nextMove.OneElementMove(el, block);
				double deltaS = SPV.SumOfAdjustedRandDistance(nextMove) - S;
				if (deltaS < bestDeltaS)
				{
					bestDeltaS = deltaS;
					bestElement = el;
					bestBlock = block;
				}
			}
		}

		if (bestDeltaS >= 0.0)
			break;

		candidate.OneElementMove(bestElement, bestBlock);
	}

	return candidate;
}
