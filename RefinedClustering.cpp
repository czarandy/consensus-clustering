#include <iostream>
#include <ctime>
#include <string>
#include <sstream>
#include <cstdlib>
#include <fstream>
#include <sstream>
using namespace std;

#include "RefinedClustering.h"
#include "SetPartition.h"
#include "SetPartitionVector.h"
#include "BestOfK.h"
#include "AverageLink.h"

void RefinedClustering::Run(const char* inputfile, double noiseCutoff)
{
	SetPartitionVector SPV;
  vector<string> elLabels;

  SPV.LoadFile(inputfile, elLabels);
  
  string output(inputfile);
  output += ".results";
  ofstream fout(output.c_str());

  AverageLink AL(noiseCutoff);
  SimulatedAnnealingOEM SA(0.999, 500, 0.0001, 100);

  cout << "Input: " << inputfile << endl << "Output: " << output << endl;
  cout << "N: " << SPV[0].GetN() << '\t'  << "K: " << SPV.size() << endl;
  cout << endl;
  fout << "Input: " << inputfile << endl << "Output: " << output << endl;
  fout << "N: " << SPV[0].GetN() << '\t'  << "K: " << SPV.size() << endl;
  fout << endl;

  SetPartition SP(AL.Run(SPV));

  cout << "Found " << SP.NumberOfBlocks() << " clusters." << endl;
  fout << "Found " << SP.NumberOfBlocks() << " clusters." << endl;
  fout << endl << "Clustered SPs" << endl;

  for (int i = 0; i < SP.NumberOfBlocks(); ++i)
  {
    fout << "C" << (i + 1) << ": ";
    for (int j = 0; j < SP.BlockSize(i + 1); ++j)
    {
      fout << SPV[SP.BlockElement(i + 1, j) - 1].GetLabel() << ' ';
    }
    fout << endl;
  }
  fout << endl;
  fout << "Consenses (only for clusters with > 2 SP)" << endl;
    
  for (int i = 0; i < SP.NumberOfBlocks(); ++i)
  {
    if (SP.BlockSize(i + 1) < 3)
      continue;

    fout << "C" << (i + 1) << ": # SP: " << SP.BlockSize(i + 1);
    fout << ", " << "# Blocks: ";
    
    SetPartitionVector Subset;
    for (int j = 0; j < SP.BlockSize(i + 1); ++j)
    {
      Subset.push_back(SPV[SP.BlockElement(i + 1, j) - 1]);
    }
    SetPartition CC(SA.Run(Subset));

    fout << CC.NumberOfBlocks() << ", " << "Avg SOD: "
      << Subset.AverageSOD(CC) << ", Lower Bound: " << Subset.CCLowerBound() <<
      endl;

    for (int b = 0; b < CC.NumberOfBlocks(); ++b)
    {
      fout << '\t';
      for (int j = 0; j < CC.BlockSize(b + 1); ++j)
      {
        fout << elLabels[CC.BlockElement(b + 1, j) - 1] << ' ';
      }
      fout << endl;
    }
  }
}
