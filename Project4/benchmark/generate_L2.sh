#!/bin/bash
# 5/1/2018
#
# This script generates benchmark results for 2 times 2 lattice. 
#

execu="ising.exe"
size=2
temp=1.0
init=1
detail=0
filename=L2

cd ../src/
make
cp ${execu} ../benchmark/
cd ../benchmark/

for mc in 10 100 1000 10000 100000
do 
    ./${execu} ${size} ${temp} ${mc} ${init} ${detail} ${filename}_${mc} 
done 

rm ${execu}
