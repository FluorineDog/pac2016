ALL = test1
CC=icc
objects=LUSolve.o LE_SymSprsMat.o
LIB=-lpthread #-xmic-avx512
CXXFLAGS=-std=c++11 -O2 -xmic-avx512 -qopt-report
test1:$(objects)
	$(CC) -o $@ $^ $(LIB) $(CXXFLAGS)
LUSolve.o:LUSolve.cpp LE_SymSprsMat.cpp
	$(CC) -c $^  $(LIB) $(CXXFLAGS)

.PHONY: clean

clean:
	rm -f $(ALL) *.o *.map log


