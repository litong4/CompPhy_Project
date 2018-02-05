This is the first project of PHY 480 (MSU). 

Use three numerical methods to solve a one-dimensional Poisson with Dirichlet boundary conditions. 

The program is compiled and run on a Linux subsystem (Ubuntu 16.04.3 LTS, Linux version 4.4.0-43-Microsoft, GCC version 5.4.0) on Windows 10 Version 1709 (OS Build 16299.214), with Armadillo 6.500.5. 

Compile: 
go to folder /src and make. 
No optimization flag is used. 

Run: 
go to folder /src and type
"./diff1d.exe n filename"
where "n" should be replaced by the number of grid points and "filename" is the prefix of names of output files. 
The program will generate output files "filename_gen.txt" from a general tri-diagonal matrix solver, "filename_spe.txt" from a special solver for this problem, and "filename_arma.txt" from Armadillo's LU decomposition. 
For n>5000 armadillo will not be used. 
For n>10000 only maximum error and time will be output. 

Benchmark: 
go to folder /benchmark and run generate_benchmark.sh 
