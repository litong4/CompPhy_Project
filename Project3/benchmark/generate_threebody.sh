#!/bin/bash
# 4/2/2018
#
# This script generates benchmark results for Earth-Jupiter-Sun system. 
#

execu="solar.exe"
time=100
stepsize=0.001
filename=threebody

cd ../src/
make
cp ${execu} ../benchmark/
cd ../benchmark/

for input in 1 2 3 4 
do
    ./${execu} ${time} ${stepsize} ${filename}${input}
done
rm ${execu}
