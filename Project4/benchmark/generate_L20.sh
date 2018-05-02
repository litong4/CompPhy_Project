#!/bin/bash
# 5/1/2018
#
# This script generates benchmark results for 20 times 20 lattice. 
#

execu="ising.exe"
size=20
init=1
detail=1
filename=L20
mc=1000

cd ../src/
make
cp ${execu} ../benchmark/
cd ../benchmark/

for init in 0 1
do
    for temp in 1 2.4 
    do
        ./${execu} ${size} ${temp} ${mc} ${init} ${detail} ${filename}_init${init}_temp${temp} 
    done
done 

rm ${execu}
