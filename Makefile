ALL = test1
test1: LUSolve.cpp LE_SymSprsMat.cpp
	g++ -std=c++11 -O3 -fopenmp -o $@ $^
debug:LUSolve.cpp LE_SymSprsMat.cpp
	g++ -std=c++11 -g -fopenmp -o $@ $^
clean:
	rm -f $(ALL) *.o *.map log

