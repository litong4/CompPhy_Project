
CXX := c++

CPPFLAGS := -std=c++11 -lm -O3 

PROG1 := solar_noclass.exe 

PROG2 := solar.exe

OBJS1 := solar_noclass.o 

OBJS2 := solar.o planet.o 

all: $(PROG1) $(PROG2) 

$(PROG1): $(OBJS1)
	$(CXX) -o $@ $^ $(CPPFLAGS) 
    
$(PROG2): $(OBJS2)
	$(CXX) -o $@ $^ $(CPPFLAGS) 
    
planet.o: planet.hpp

solar.o: planet.hpp planet.cpp

clean:
	-rm -f *.exe *.o 
