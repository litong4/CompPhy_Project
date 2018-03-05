#!/bin/bash
# 3/5/2018
#
# This script generates benchmark results under this folder. 
# Number of points n=10,20,40,80
# Filenames are beam_n10, beam_n20, ...

# For ho.exe: 
# "one" and "two" in filenames indicate one- and two-electron case, respectively. 
# for two-electron case, \omega_r is also shown in the filename as "wr..."
#
# "_arma" and "_jacobi" are added to the end of the file names 
# for armadillo and jacobi diagonalization, respectively. 
#

execu="beam.exe"

cd ../src/
make
cp ${execu} ../benchmark/
cd ../benchmark/

for num in 10 20 40 80 
do
    ./${execu} ${num} beam_n${num}
done
rm ${execu}
