
CXX := c++

CPPFLAGS := -std=c++11 -lm -larmadillo -lblas -llapack -O3 

PROG1 := beam.exe

OBJS1 := beam.o eig_solver.o eig_output.o

PROG2 := ho.exe

OBJS2 := ho.o eig_solver.o eig_output.o

PROG3 := ws.exe

OBJS3 := ws.o eig_solver.o eig_output.o

all: $(PROG1) $(PROG2) 

$(PROG1): $(OBJS1)
	$(CXX) -o $@ $^ $(CPPFLAGS)
 
$(PROG2): $(OBJS2)
	$(CXX) -o $@ $^ $(CPPFLAGS)
    
$(PROG3): $(OBJS3)
	$(CXX) -o $@ $^ $(CPPFLAGS)    
    
beam.o: eig_solver.hpp eig_output.hpp

ho.o: eig_solver.hpp eig_output.hpp

ws.o: eig_solver.hpp eig_output.hpp

eig_solver.o: eig_solver.hpp 

eig_output.o: eig_output.hpp
    
test: UnitTest.o eig_solver.o
	$(CXX) -o $@.exe $^ $(CPPFLAGS)

UnitTest.o: eig_solver.o catch.hpp 

clean:
	-rm -f *.exe *.o
