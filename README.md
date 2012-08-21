Set Partition & Consensus Clustering Library

Code written as part of research described in the paper [Consensus Clustering Algorithms: Comparison & Refinement](http://www.cs.ucdavis.edu/~filkov/papers/consensuseval.pdf) (Andrey Goder & Vladimir Filkov).

I haven't updated it since it was originally written in 2007, apart from ensuring it compiles. No guarantees.

Included are files that implement the following:
  1. A SetPartition class.
  2. A SetPartitionVector class that stores an arbitrary amount of SetPartitions and can generate them from a file, by choosing randomly, or by creating "noisy partitions".
  3. A Timer class allowing for timing of performance of algorithms
  4. The Mersenne Twister pseudo-random number generator, for generating random set partitions and for the randomized algorithms
  5. A Matrix class, used in some of the algorithms.
  6. A suite of consensus clustering algorithms (described in the paper):
    - BestOfK
    - BestOneElementMove
    - SimulatedAnnealingOEM
    - MajorityRule
    - CCPivot
    - CCAverageLink
  7. A refined consensus clustering algorithm using AverageLink clustering that breaks up a set of clustering into smaller groups that will have better consenses.

The Main.cpp file includes examples showing some of these uses.
