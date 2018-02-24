
CC := c++

CFLAGS := -std=c++11 -lm -larmadillo -lblas -llapack -O3 

PROG := eig.exe

OBJS := eig.o eig_solver.o

all: $(PROG) test

$(PROG): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

eig.o: eig.cpp eig_solver.hpp 
	$(CC) -c $< $(CFLAGS)

eig_solver.o: eig_solver.cpp eig_solver.hpp 
	$(CC) -c $< $(CFLAGS)
    
test: UnitTest.o eig_solver.o
	$(CC) -o $@.exe $^ $(CFLAGS)

UnitTest.o: UnitTest.cpp eig_solver.o catch.hpp 
	$(CC) -c $< $(CFLAGS)

clean:
	-rm -f *.exe *.o
