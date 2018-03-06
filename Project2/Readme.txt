This is the 2nd project of PHY 480 (MSU). 

Eigenvalue decomposition of a real symmetric matrix, using Jacobi's method. 

The program is compiled and run on a Linux subsystem (Ubuntu 16.04.3 LTS, Linux version 4.4.0-43-Microsoft, GCC version 5.4.0) on Windows 10 Version 1709 (OS Build 16299.214), with Armadillo 6.500.5. 

Compile: 
go to folder /src and make (or "make test" to get unit tests). 

Run: 
go to folder /src and type
"./test.exe" for unit tests. 
"./beam.exe n filename"
"./ho.exe n filename"
where "n" should be replaced by the number of grid points and "filename" is the prefix of names of output files. 
In beam.exe, we will diagonalize a tri-diagonal symmetric matrix, with d=2/h^2 and a=-1/h^2.
In ho.exe, the quantum dot problem with one and two electrons trapped inside a harmonic oscillator potential is calculated. For two-electron case, cases with wr=0.01, 0.5, 1, 5 are computed. 

Benchmark: 
go to folder /benchmark and run generate_benchmark*.sh
