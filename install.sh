#!/bin/bash

#lighter compilation
echo 'Starting Lighter installation'
cd Lighter
make

cd ..

echo 'Starting bcalm2 installation'
#bcalm compilation
cd bcalm
mkdir build;  cd build;  cmake ..;  make -j 8

cd ../..

echo 'Starting REINDEER installation'
cd REINDEER
make

cd ..

Rscript install.R
