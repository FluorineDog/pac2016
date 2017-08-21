ALL = test1
objects=LUSolve.o LE_SymSprsMat.o
test1:$(objects)
	g++ -O2 -fopenmp -o $@ $^
LUSolve.o:LUSolve.cpp LE_SymSprsMat.cpp
	g++ -O2 -fopenmp -c $^

clean:
	rm -f $(ALL) *.o *.map log

