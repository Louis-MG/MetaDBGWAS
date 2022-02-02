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

R install.packages("ape")
R install.packages("phangorn")
R install.packages("https://raw.githubusercontent.com/sgearle/bugwas/master/build/bugwas_1.0.tar.gz", repos=NULL, type="source")

