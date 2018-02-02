#!/bin/bash

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

