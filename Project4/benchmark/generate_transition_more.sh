#!/bin/bash
# 5/1/2018
#
# This script generates benchmark results for the study of phase transition. 
#

execu="ising.exe"
init=1
detail=0
filename=tran
mc=1000000

cd ../src/
make
cp ${execu} ../benchmark/
cd ../benchmark/

for size in 100 140 
do
    for temp1 in 2.1 2.2 2.3 
    do
        for temp2 in "" 2 4 6 8 
        do
        temp=${temp1}${temp2}
        ./${execu} ${size} ${temp} ${mc} ${init} ${detail} ${filename}_size${size}_temp${temp} 
        done
    done
    temp=2.4
    ./${execu} ${size} ${temp} ${mc} ${init} ${detail} ${filename}_size${size}_temp${temp} 
done 

rm ${execu}
