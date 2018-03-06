#!/bin/bash
# 3/5/2018
#
# This script generates benchmark results under this folder. 
# The maximum rmax=4,6,8,10,15,20,30,40,50
# The number of points n=200 300 400 500 1000
# The corresponding filenames are r4n500, r4n1000, ...
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

for rmax in 4 6 8 10 15 20 30 40 50 60 
do
    for num in 200 300 400 500 1000
    do
        ./${execu} ${num} ${rmax} r${rmax}n${num}
    done
done
rm ${execu}
