#include "SetPartition.h"
#include "Utility.h"

#include <queue>
#include <cstdlib>
using namespace std;

MersenneTwister SetPartition::mt;

// Prints a block of a set partition in the form {a, b, c}
namespace std
{
	ostream& operator<<(ostream& os, const vector<int>& v)
	{
		if (v.size() == 0)
			return os;

		os << '{';
		ostream_iterator<int> osi(cout, ", ");
		copy(v.begin(), v.end() - 1, osi);
		os << v.back();
		return os << '}';
	}
}

// Prints a set partition in the standard form: {{a, b, c}, {d, e, f}, {g, h}}
ostream& operator<<(ostream& os, const SetPartition& sp)
{
	os << '{';
	ostream_iterator< vector<int> > osi(cout, ", ");
	copy(sp.blocks.begin() + 1, sp.blocks.end() - 1, osi);
	return os << sp.blocks.back() << '}';
}

// Compares two blocks by looking at the first elements
bool operator<(const vector<int>& lhs, const vector<int>& rhs)
{
	if (rhs.size() == 0)
		return false;
	else if (lhs.size() == 0)
		return true;
	return lhs[0] < rhs[0];
}

void SetPartition::SetUpBlocks(unsigned int blockCount)
{
	blocks.clear();
	blocks.resize(blockCount + 1);
	
	for (unsigned int i = 1; i < elements.size(); ++i)
		blocks[elements[i]].push_back(i);

}

// Randomizes partition
void SetPartition::Randomize()
{	
	SetUpBlocks(RandomEquivalenceRelation());
}

// Sets this partitions to be that induced by the two given partitions
// where 1 is clustered with 2 iff it is in both a and b
void SetPartition::Induce(const SetPartition& a, const SetPartition& b)
{
	if (a.elements.size() != b.elements.size() || a.elements.size() != elements.size())
		return;

	// Set up sorted element array
	unsigned int* sortedElements = new unsigned int[elements.size()];
	iota(sortedElements, sortedElements + elements.size(), 0);

	sort(sortedElements + 1, sortedElements + elements.size(), InducedElementComparer(a, b));

	unsigned int block = 1;
	for (unsigned int i = 2; i < elements.size(); ++i)
	{
		unsigned int prevEl = sortedElements[i - 1];
		unsigned int el = sortedElements[i];
		if (a.elements[prevEl] == a.elements[el] && b.elements[prevEl] == b.elements[el])
			elements[el] = elements[prevEl] = block;
		else
			elements[el] = ++block;
	}
	SetUpBlocks(block);

	delete [] sortedElements;
}

// Computes the absolute Rand distance between two set partitions
unsigned int SetPartition::RandDistanceTo(const SetPartition& sp) const
{
	
	if (elements.size() != sp.elements.size())
		return 0;

	unsigned int b = accumulate(blocks.begin() + 1, blocks.end(), 0, BlockSizeChooseTwo);
	unsigned int c = accumulate(sp.blocks.begin() + 1, sp.blocks.end(), 0, BlockSizeChooseTwo);

	SetPartition induced(this->GetN());
	induced.Induce(*this, sp);
	unsigned int a = accumulate(induced.blocks.begin() + 1, induced.blocks.end(), 0, BlockSizeChooseTwo);

	return b + c - 2 * a;
}

double SetPartition::NormalizedRandDistanceTo(const SetPartition& sp) const
{
	return static_cast<double>(RandDistanceTo(sp)) / static_cast<double>(ChooseTwo(GetN()));
}

double SetPartition::AdjustedRandDistanceTo(const SetPartition& sp) const
{
	if (elements.size() != sp.elements.size())
		return 0;

	SetPartition induced = SetPartition::Induced(*this, sp);

	unsigned int NChooseTwo = ChooseTwo(GetN());

	unsigned int a = accumulate(induced.blocks.begin() + 1, induced.blocks.end(), 0, BlockSizeChooseTwo);
	unsigned int b = accumulate(blocks.begin() + 1, blocks.end(), 0, BlockSizeChooseTwo) - a;
	unsigned int c = accumulate(sp.blocks.begin() + 1, sp.blocks.end(), 0, BlockSizeChooseTwo) - a;
	unsigned int d = NChooseTwo - a - b - c;

	double nc = static_cast<double>((a + b) * (a + c)  + (c + d) * (b + d)) / static_cast<double>(NChooseTwo);

	if (NChooseTwo == nc)
		return 0.0;
	return 1.0 - ((((a + d - nc) / (NChooseTwo - nc)) + 1.0) / 2.0);
}

double SetPartition::NormalRandDistanceTo(const SetPartition& sp) const
{
	unsigned int N = GetN();
	double B = BellNumber<double>(N - 1) / BellNumber<double>(N);
  //B = log(N - 1) / (N - 1);
	double p = 2 * B * (1 - B);
	double EX = p * ChooseTwo(N);
	double VX = EX * (1 - p);
	double d = static_cast<double>(RandDistanceTo(sp));
	return (d - EX) / sqrt(VX);
}

// Sets up a random equivalence relation on elements, returns numbers of blocks
int SetPartition::RandomEquivalenceRelation() 
{
	unsigned int n = elements.size() - 1;
	long double* b = new long double[elements.size() + 1];
	
	b[1] = 1.0;

	long double sum1;
	for (int l = 1; l <= n - 1; ++l)
	{
		sum1 = 1 / static_cast<long double>(l);
		for (int k = 1; k <= l - 1; ++k)
			sum1 = (sum1 + b[k]) / static_cast<long double>(l - k);
		b[l + 1] = (sum1 + b[l]) / static_cast<long double>(l + 1);
	}

	unsigned int m = n;
	unsigned int npart = 0;

	while (true)
	{
		long double z = mt.RandomReal() * b[m] * static_cast<long double>(m);
		int k = 0;
		++npart;

		while (0 <= z)
		{
			elements[m] = npart;
			--m;

			if (m == 0)
				break;

			z = z - b[m];
			++k;
			z = z * k;
		}

		if (m == 0)
			break;
	}

	for (int m = 1; m < n; ++m)
		swap(elements[mt.RandomInt(m, n)], elements[m]); // Randomly permutes elements array

	return npart;
}

void SetPartition::MakeFirst()
{
	fill(elements.begin(), elements.end(), 1);
	SetUpBlocks(1);
}

void SetPartition::MakeLast()
{
	iota(elements.begin(), elements.end(), 0);
	SetUpBlocks(elements.size() - 1);
}

bool SetPartition::operator==(const SetPartition& sp)
{
	elements[0] = sp.elements[0];
	return elements == sp.elements;
}

void SetPartition::Next()
{ 
	vector<int>::iterator i, first = elements.begin() + 1;
	for (i = elements.end() - 1; i > first; --i)
	{
		if (*i <= *max_element(first, i))
		{
			++*i;
			fill(i + 1, elements.end(), 1);
			break;
		}
	}

	if (i == first)
		MakeFirst(); // We hit the end
	else
		SetUpBlocks(*max_element(first, elements.end()));
}

// Returns true iff a move ocurred
bool SetPartition::OneElementMove(int element, int block)
{
	if (block == elements[element])
		return false; // Nothing to be done

	vector< vector<int> > v(blocks);

	// Remove element from old block
	blocks[elements[element]].erase(remove_if(blocks[elements[element]].begin(), blocks[elements[element]].end(), 
		bind2nd(equal_to<int>(), element)), blocks[elements[element]].end());

	if (block != blocks.size())
		blocks[block].push_back(element);
	else
	{
		blocks.push_back(vector<int>());
		blocks.back().push_back(element);
	}

	swap(elements[element], block); // block is now the old block

	// Check for empty block
	if (blocks[block].size() == 0)
	{
		// Move the empty block to the end, and update references
		swap(blocks[block], blocks.back());
		for (int i = 1; i < elements.size(); ++i)
			if (elements[i] == blocks.size() - 1)
				elements[i] = block;
		blocks.pop_back();
	}

	return true;
}

void SetPartition::RandomOneElementMove()
{
	int element = mt.RandomInt(1, elements.size() - 1);
	int block = elements[element];

	while (block == elements[element])
		block = mt.RandomInt(1, blocks.size());

	OneElementMove(element, block);
}

void SetPartition::RandomOneElementMoves(int n)
{
  vector<int> sizes(blocks.size() * 2);
  int blockCount = blocks.size() - 1;
  for (int i = 1; i < blocks.size(); ++i)
    sizes[i] = blocks[i].size();
  while (n--)
  {
    int element = mt.RandomInt(1, elements.size() - 1);
    int block = elements[element];

    while (block == elements[element])
      block = mt.RandomInt(1, blockCount);

    --sizes[elements[element]];
    if (block == blockCount && sizes[elements[element]] == 0)
    {
      ++sizes[elements[element]];
      continue; // this move is pointless
    }
    else if (block == blockCount && sizes[elements[element]] > 0)
    {
      sizes[blockCount] = 1;
      ++blockCount; 
    }
    else if (block != blockCount && sizes[elements[element]] == 0)
    {
      --blockCount;
      for (int i = 1; i <= elements.size() - 1; ++i)
      {
        if (elements[i] == blockCount) // move from last block into empty block
          elements[i] = elements[element];
      }
      sizes[elements[element]] = sizes[blockCount];
      sizes[blockCount] = 0;
    }
    else
      ++sizes[block];
    elements[element] = block;
  }
  SetUpBlocks(blockCount);
}

void SetPartition::RandomTransition(int& element, int& block)
{
	element = mt.RandomInt(1, elements.size() - 1);
	block = mt.RandomInt(1, blocks.size());
}

bool SetPartition::IsFirst()
{
	return blocks.size() == 2; 
}

bool SetPartition::IsLast()
{
	return blocks.size() == elements.size();
}

long double SetPartition::ProbabilityOfBlockConfiguration()
{
	int n = elements.size() - 1;

	map<unsigned int, unsigned int> m;
	for (int i = 1; i < blocks.size(); ++i)
		++m[blocks[i].size()];

	return BlockDistributionProbability(m, n);

	vector<unsigned int> blockSizes(blocks.size() - 1);
	for (int i = 1; i < blocks.size(); ++i)
		blockSizes[i - 1] = blocks[i].size();

	return MaxBlockProbability(n, *max_element(blockSizes.begin(), blockSizes.end()));

	sort(blockSizes.begin(), blockSizes.end());

	double product = 1.0;

	for (int i = 1; i < blocks.size(); ++i)
	{
		product *= MaxBlockProbability(n, blockSizes[i - 1]);
		n -= blockSizes[i - 1];
	}
	return product;
}

unsigned int SetPartition::LargestBlockSize()
{
	unsigned int max = 0;
	for (int i = 1; i < blocks.size(); ++i)
		if (blocks[i].size() > max)
			max = blocks[i].size();

	return max;
}

void SetPartition::PrintDistributionOfBlockSizes()
{
	vector<int> blockSizes(blocks.size() - 1);
	for (int i = 1; i < blocks.size(); ++i)
		blockSizes[i - 1] = blocks[i].size();

	sort(blockSizes.begin(), blockSizes.end(), greater<int>());
	cout << blockSizes;
}

bool SetPartition::CoClustered(int a, int b) const
{
	int n = elements.size() - 1;
	if (a <= n && b <= n && a >= 1 && b >= 1)
		return elements[a] == elements[b];
	return false;
}

void SetPartition::MergeBlocks(int a, int b)
{
	// try this way, but it might be inefficient
	if (BlockSize(b) < BlockSize(a)) 
		swap(a, b);

  queue<int> q;
  for (int i = 0; i < BlockSize(a); ++i)
    q.push(blocks[a][i]);

  while (!q.empty())
  {
    OneElementMove(q.front(), b);
    q.pop();
  }
  return;
}

void SetPartition::RandomElementSwaps(int N)
{
  if (GetN() < 2) return ;
  for (int i = 0; i < N; ++i)
  {
	  int element = mt.RandomInt(1, elements.size() - 1);
    int element2 = element;
	  while (element2 == element)
      element2 = mt.RandomInt(1, elements.size() - 1);

    swap(elements[element], elements[element2]);
  }
  SetUpBlocks(NumberOfBlocks());
}

void SetPartition::RandomElementTripPerms(int N)
{
  if (GetN() < 2) return ;
  for (int i = 0; i < N; ++i)
  {
	  int element = mt.RandomInt(1, elements.size() - 1);
    int element2 = element;
	  while (element2 == element)
      element2 = mt.RandomInt(1, elements.size() - 1);
    int element3 = element;
    while (element3 == element || element3 == element2)
      element3 = mt.RandomInt(1, elements.size() - 1); 

    switch (mt.RandomInt(1, 6))
    {
      case 1: swap(elements[element], elements[element2]); break;
      case 2: swap(elements[element], elements[element3]); break;
      case 3: swap(elements[element2], elements[element3]); break;
      case 4: swap(elements[element], elements[element3]); 
              swap(elements[element2], elements[element3]); break;
      case 5: swap(elements[element], elements[element2]);
              swap(elements[element2], elements[element3]); break;
    }
  } 
  
  SetUpBlocks(NumberOfBlocks());
}

void SetPartition::LoadData(string data)
{
  fill(elements.begin(), elements.end(), -1);
  int curBlock = 1;
  int n;
  char c;
  istringstream istr(data);
  while (istr >> n)
  {
    elements[n] = curBlock;
    istr >> c;
    if (c == ';')
      ++curBlock;
    else if (c == '.')
      break;
  }
  SetUpBlocks(curBlock);

  // test this partition
  for (int i = 1; i <= GetN(); ++i)
  {
    if (elements[i] == -1)
    {
      cout << "element " << i << " has -1!" << endl;
      exit(1);
    }
    
  }
}

bool SetPartition::Verify() const
{
  int n = GetN();
  int b = NumberOfBlocks();

  for (int i = 1; i <= n; ++i)
    if (WhichBlock(i) < 1 || WhichBlock(i) > b)
      return false;

  return true;
}

SetPartition SetPartition::LoadClusterFile(const char* filename)
{
  ifstream fin(filename);
  vector<int> v;
  int i;
  while (fin >> i) 
    v.push_back(i);

  int blockCount = *max_element(v.begin(), v.end());

  SetPartition sp(v.size());
  copy(v.begin(), v.end(), sp.elements.begin() + 1);
  sp.SetUpBlocks(blockCount);
  return sp;
}
