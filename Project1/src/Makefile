
CC := c++

CFLAGS := -std=c++11 -lm -larmadillo -lblas -llapack -O3

diff1d.exe: diff1d.cpp
	$(CC) -o diff1d.exe ./diff1d.cpp $(CFLAGS)

clean:
	-rm *.exe
