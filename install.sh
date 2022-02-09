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

mkdir DBGWAS
cd DBGWAS
curl  https://www.dropbox.com/s/s9oojqfl1kgi4l5/DBGWAS-0.5.4-Linux-precompiled.tar.gz?dl=1 -L --output ./dbgwas.Linux.precompiled.tar.gz
tar -zxvf dbgwas.Linux.precompiled.tar.gz
rm dbgwas.Linux.precompiled.tar.gz
cd ..

echo "Installation complete ! Do not forget to install the R libraries."
