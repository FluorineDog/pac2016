ALL = test1
objects=LUSolve.o LE_SymSprsMat.o
test1:$(objects)
	icc -o $@ $^
LUSolve.o:LUSolve.cpp LE_SymSprsMat.cpp
	icc -c $^

clean:
	rm -f $(ALL) *.o *.map log

