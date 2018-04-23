
CXX := c++

CPPFLAGS := -std=c++11 -lm -O3 

PROG1 := ising.exe

OBJS1 := ising.o

all: $(PROG1)

$(PROG1): $(OBJS1)
	$(CXX) -o $@ $^ $(CPPFLAGS) 

clean:
	-rm -f *.exe *.o 
