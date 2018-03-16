
CXX := c++

CPPFLAGS := -std=c++11 -lm -O3 

PROG1 := solar.exe 

OBJS1 := solar.o 

$(PROG1): $(OBJS1)
	$(CXX) -o $@ $^ $(CPPFLAGS) 

clean:
	-rm -f *.exe *.o
