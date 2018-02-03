#!/bin/bash
# 2/2/2018
#
# This script generates all the benchmark results under this folder. 
# The number of points n=10, 100, 1000, ..., 10^8. 
# The corresponding file names are n10e1, n10e2, n10e3, ..., n10e8. 
#
# "_gen", "_spe" and "_arma" are added to the end of the file names 
# for general tridiagonal matrix solver, special solver (for this project)
# and armadillo LU decomposition. 
#
# For n>5000 armadillo will not be used. 
# For n>10000 only maximum error and time will be output. 
#

execu="diff1d.exe"

cd ../src/
make
cp ${execu} ../benchmark/
cd ../benchmark/

for power in 1 2 3 4 5 6 7 8 
do
    num=1
    for ((loop=0; loop<${power}; loop++)); do
        num=${num}0
    done
    ./${execu} ${num} n10e${power}
done
rm ${execu}

