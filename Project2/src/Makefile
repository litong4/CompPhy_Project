
CC := c++

CFLAGS := -std=c++11 -lm -larmadillo -lblas -llapack -O3

eig.exe: eig.cpp
	$(CC) -o eig.exe ./eig.cpp $(CFLAGS)

clean:
	-rm *.exe
