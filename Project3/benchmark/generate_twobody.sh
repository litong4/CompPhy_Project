#!/bin/bash
# 4/2/2018
#
# This script generates benchmark results for Earth-Sun system. 
#

execu="solar.exe"
time=10
stepsize=0.001
filename=twobody

cd ../src/
make
cp ${execu} ../benchmark/
cd ../benchmark/

for input in 1 2 3
do
    ./${execu} ${time} ${stepsize} ${filename}${input}
done
rm ${execu}
