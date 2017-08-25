ALL = test1
objects=LUSolve.o LE_SymSprsMat.o
test1: LUSolve.cpp LE_SymSprsMat.cpp
	g++ -std=c++11 -O3 -fopenmp -o $@ $^

clean:
	rm -f $(ALL) *.o *.map log

