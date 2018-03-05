#!/bin/bash
# 3/5/2018
#
# This script generates benchmark results under this folder. 
# The maximum rmax=4
# The number of points n=10,20,40,80,100,200
# The corresponding filenames are r4n10, r4n20, ...
#
# For ho.exe: 
# "one" and "two" in filenames indicate one- and two-electron case, respectively. 
# for two-electron case, \omega_r is also shown in the filename as "wr..."
#
# "_arma" and "_jacobi" are added to the end of the file names 
# for armadillo and jacobi diagonalization, respectively. 
#

execu="ho.exe"

cd ../src/
make
cp ${execu} ../benchmark/
cd ../benchmark/

rmax=4
for num in 10 20 40 80 100 200
do
    ./${execu} ${num} ${rmax} r${rmax}n${num}
done
rm ${execu}
