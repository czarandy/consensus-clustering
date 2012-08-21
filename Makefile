CC = g++
CFLAGS = -O3 -ffast-math -fomit-frame-pointer -msse -mmmx
DIR = .

a.out : AdjBestOfK.o AverageLink.o BestOfK.o CCOptimal.o Main.o MajorityRule.o MersenneTwister.o SetPartition.o SetPartitionVector.o Utility.o CCPivot.o CCAverageLink.o RefinedClustering.o
	$(CC) $(CFLAGS) -o $(DIR)/sp.out $(DIR)/*.o

AdjBestOfK.o : AdjBestOfK.cpp AdjBestOfK.h CCHeuristic.h
	$(CC) $(CFLAGS) -c -o $(DIR)/AdjBestOfK.o AdjBestOfK.cpp

AverageLink.o : AverageLink.h AverageLink.cpp
	$(CC) $(CFLAGS) -c -o $(DIR)/AverageLink.o AverageLink.cpp

BestOfK.o : BestOfK.cpp BestOfK.h MajorityRule.h CCHeuristic.h
	$(CC) $(CFLAGS) -c -o $(DIR)/BestOfK.o BestOfK.cpp

CCOptimal.o : CCOptimal.cpp CCOptimal.h CCHeuristic.h
	$(CC) $(CFLAGS) -c -o $(DIR)/CCOptimal.o CCOptimal.cpp

CCAverageLink.o : CCAverageLink.cpp CCAverageLink.h
	$(CC) $(CFLAGS) -c -o $(DIR)/CCAverageLink.o CCAverageLink.cpp

Main.o : Main.cpp MersenneTwister.h SetPartition.h SetPartitionVector.h Utility.h BestOfK.h AdjBestOfK.h MajorityRule.h
	$(CC) $(CFLAGS) -c -o $(DIR)/Main.o Main.cpp

MajorityRule.o : MajorityRule.h MajorityRule.cpp SetPartition.h CCHeuristic.h SetPartitionVector.h
	$(CC) $(CFLAGS) -c -o $(DIR)/MajorityRule.o MajorityRule.cpp

MersenneTwister.o : MersenneTwister.cpp MersenneTwister.h
	$(CC) $(CFLAGS) -c -o $(DIR)/MersenneTwister.o MersenneTwister.cpp

SetPartition.o : SetPartition.h SetPartition.cpp Utility.h MersenneTwister.h
	$(CC) $(CFLAGS) -c -o $(DIR)/SetPartition.o SetPartition.cpp

SetPartitionVector.o : SetPartitionVector.cpp SetPartitionVector.h SetPartition.h Utility.h
	$(CC) $(CFLAGS) -c -o $(DIR)/SetPartitionVector.o SetPartitionVector.cpp

Utility.o : Utility.h Utility.cpp
	$(CC) $(CFLAGS) -c -o $(DIR)/Utility.o Utility.cpp

CCPivot.o : CCPivot.h CCPivot.cpp
	$(CC) $(CFLAGS) -c -o $(DIR)/CCPivot.o CCPivot.cpp

RefinedClustering.o : RefinedClustering.h RefinedClustering.cpp
	$(CC) $(CFLAGS) -c -o $(DIR)/RefinedClustering.o RefinedClustering.cpp

clean :
	rm $(DIR)/*.o
	rm $(DIR)/sp.out
