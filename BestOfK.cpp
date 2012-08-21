#include "BestOfK.h"
#include "CCAverageLink.h"
#include "Matrix.h"

const SetPartition& BestOfK::FindBestOfK(const SetPartitionVector& SPV)
{
	unsigned int min = SPV.SumOfRandDistance(SPV[0]);
	unsigned int minSP = 0;

	for (int i = 1; i < SPV.size(); ++i)
	{
		unsigned int d = SPV.SumOfRandDistance(SPV[i]);
		if (d < min)
		{
			min = d;
			minSP = i;
		}
	}
	return SPV[minSP];
}

SetPartition BestOfK::Run(const SetPartitionVector& SPV)
{
	return FindBestOfK(SPV);
}

void BestOneElementMove::GenerateK(const SetPartitionVector& SPV, Matrix<int>& K)
{
  int n = K.GetRows();
	unsigned int k = SPV.size();

	for (int i = 1; i < n; ++i)
	{
		for (int j = 1; j < n; ++j)
		{
			K(i, j) = k;
			for (int p = 0; p < k; ++p)
			{
				if (SPV[p].CoClustered(i, j))
					K(i, j) -= 2;
			}
		}
	}
}

void BestOneElementMove::GenerateM(const SetPartitionVector& SPV, const
    SetPartition& sp, const Matrix<int>& K, Matrix<int>& M)
{
	unsigned int b = sp.NumberOfBlocks();
  int n = M.GetRows();

  M.FillWithValue(0);

	for (int i = 1; i < n; ++i)
	{
		for (int j = 1; j <= b; ++j)
		{
			int blockSize = sp.BlockSize(j);
			for (int l = 0; l < blockSize; ++l)
			{
				int el = sp.BlockElement(j, l);
				if (el != i)
					M(i, j) += K(i, el);
			}
		}
	}
}

SetPartition BestOneElementMove::Run(const SetPartitionVector& SPV, bool AL)
{

  if (AL)
  {
    CCAverageLink AL(0.5);;
    return Run(SPV, AL.Run(SPV));
  }
  else
    return Run(SPV, FindBestOfK(SPV));
}

SetPartition BestOneElementMove::Run(const SetPartitionVector& SPV,
    SetPartition sp)
{

  SetPartition candidate(sp);
	unsigned int n = candidate.GetN();

  Matrix<int> K(n + 1, n + 1);
  Matrix<int> M(n + 1, n + 1);
	GenerateK(SPV, K);
	GenerateM(SPV, candidate, K, M);

	unsigned int b = candidate.NumberOfBlocks();

  vector<int> mb(n + 1);
  vector<int> mv(n + 1);

	// Setup mb and mv
	for (int i = 1; i <= n; ++i)
	{
    int m = M.GetMinRowElementColumn(i, 1, b + 1);
		mv[i] = M(i, m);
		mb[i] = m;
	}

	// Main loop
	while (true)
	{
		// Find BOEM
		int el = 1;
		int delta = M(1, candidate.WhichBlock(1)) - mv[1];

		for (int x = 1; x <= n; ++x)
		{
			int newDelta = M(x, candidate.WhichBlock(x)) - mv[x];
			if (newDelta > delta)
			{
				delta = newDelta;
				el = x;
			}
		}


		int blockfrom = candidate.WhichBlock(el);
		int blockto =  mb[el];

		if (delta <= 0 || blockfrom == blockto)
			break;

		for (int i = 1; i <= n; ++i) 
		{
			M(i, blockfrom) -= K(i, el);
			M(i, blockto) += K(i, el);
      
			if (M(i, blockfrom) < mv[i]) {
				mv[i] = M(i, blockfrom);
				mb[i] = blockfrom;
			}
			if (M(i, blockto) < mv[i]) {
				mv[i] = M(i, blockto);
				mb[i] = blockto;
			}
		}

		if (candidate.BlockSize(blockfrom) == 1)
		{
			for (int i = 1; i <= n; ++i)
			{
				M(blockfrom, i) = M(b, i);
				M(b, i) = 0;
			}
		}

		candidate.OneElementMove(el, blockto);
		b = candidate.NumberOfBlocks();
	}

	return candidate;
}


// SIMULATED ANNEALING

int SimulatedAnnealingOEM::DeltaCost(const SetPartition& P, const Matrix<int>& K, int e, 
    int blockfrom, int blockto)
{
	int suma = 0, sumb = 0;

	if (blockfrom == blockto) return 0;

	// for each element e in blockfrom
	// calculate suma(K[e][])
	for (int i = 0; i < P.BlockSize(blockfrom); ++i)
		suma += K(e, P.BlockElement(blockfrom, i));

	for (int i = 0; i < P.BlockSize(blockto); ++i)
		sumb += K(e, P.BlockElement(blockto, i));

	return suma - sumb;
}

SetPartition SimulatedAnnealingOEM::Run(const SetPartitionVector& SPV, bool AL)
{
  if (AL)
  {
    CCAverageLink AL(0.5);
    return Run(SPV, AL.Run(SPV));
  }
  else
    return Run(SPV, FindBestOfK(SPV));
}

SetPartition SimulatedAnnealingOEM::Run(const SetPartitionVector& SPV,
    SetPartition sp)
{
  unsigned int n = SPV[0].GetN();
  Matrix<int> K(n + 1, n + 1);
	GenerateK(SPV, K);

  SetPartition P(sp);

  SetPartition Best(P);
  long int mindS = 1000000000;
	long int dS = 0; // dS (delta sum of distance) after one move
	int element, blockfrom, blockto; // element to move, and blocks to move it from/to

	double T = InitTemp(); // temperature
	long double  simanncond; // condition and the random number for the sim ann decision 
	long int numcond = 0;
	long int SumdS = 0;
	long double Const;


	// determine the constant
	// such that in the beginning
	// > 80% of the bad moves are accepted
	for (int i = 0; i < 15000; ++i)
	{
		P.RandomTransition(element, blockto);
		blockfrom = P.WhichBlock(element);
		dS = DeltaCost(P, K, element, blockfrom, blockto);
		if (dS < 0)
		{
		  SumdS -= dS;
		  ++numcond;
		}
	}

  Const = (long double ) SumdS / (long double ) numcond;

  long int dsSum = 0;
  do
  {
    for (int i = 1; i <= IterLength; ++i)
	  {
		do {
			P.RandomTransition(element, blockto);
			blockfrom = P.WhichBlock(element);
		} while ((blockto == blockfrom) || ((P.BlockSize(blockfrom) + P.BlockSize(blockto)) == 1 ));

        dS = DeltaCost(P, K, element, blockfrom, blockto);

    dsSum += dS;
    if (dsSum < mindS)
    {
      Best = P;
      dsSum = mindS;
    }

		simanncond = exp(dS / (Const * T));

		if ( dS > 0 ) 
		{
			P.OneElementMove(element, blockto);
		}
        else 
		{
			if (simanncond > P.RandomReal()) 
			{
				P.OneElementMove(element, blockto);
			}
		}

    }
    T = ReduceTemp (T); 
  }
  while (T > FinalTemp);

  return Best;
}
