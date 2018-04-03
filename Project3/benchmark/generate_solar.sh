#!/bin/bash
# 4/2/2018
#
# This script generates benchmark results for the whole solar system 
#

execu="solar.exe"
time=200
stepsize=0.001
filename=solar

cd ../src/
make
cp ${execu} ../benchmark/
cd ../benchmark/

./${execu} ${time} ${stepsize} ${filename}

rm ${execu}
