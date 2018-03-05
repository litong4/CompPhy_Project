#!/bin/bash
# 3/5/2018
#
# This script generates benchmark results under this folder. 
# The maximum rmax=1,2,3,4,5
# Number of points n is selected to make grid spacing be always 0.01. 
# The corresponding file names are r1n100, r2n200, ...
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

for rmax in 1 2 3 4 5
do
    num=${rmax}00
    ./${execu} ${num} ${rmax} r${rmax}n${num}
done
rm ${execu}
