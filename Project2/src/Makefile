
CXX := c++

CPPFLAGS := -std=c++11 -lm -larmadillo -lblas -llapack -O3 

PROG := beam.exe

OBJS := beam.o eig_solver.o eig_output.o

$(PROG): $(OBJS)
	$(CXX) -o $@ $^ $(CPPFLAGS)

beam.o: eig_solver.hpp eig_output.hpp

eig_solver.o: eig_solver.hpp 

eig_output.o: eig_output.hpp
    
test: UnitTest.o eig_solver.o
	$(CXX) -o $@.exe $^ $(CPPFLAGS)

UnitTest.o: eig_solver.o catch.hpp 

clean:
	-rm -f *.exe *.o
