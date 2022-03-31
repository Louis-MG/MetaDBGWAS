#!/bin/bash

GREEN='\132[1;32m' # green color
NC='\033[0m' # No Color
GREEN_BACKGROUND='\e[42m' #self explanatory right?

#lighter compilation
echo "${GREEN}Starting Lighter installation${NC}"
cd Lighter
make

cd ..

echo "${GREEN}Starting bcalm2 installation${NC}"
#bcalm compilation
cd bcalm
mkdir build;  cd build;  cmake ..;  make -j 4

cd ../..

echo "${GREEN}Starting REINDEER installation${NC}"
cd REINDEER
make

cd ..

echo "${GREEN}Starting DBGWAS installation${NC}"
mkdir DBGWAS
cd DBGWAS
curl  https://www.dropbox.com/s/s9oojqfl1kgi4l5/DBGWAS-0.5.4-Linux-precompiled.tar.gz?dl=1 -L --output ./dbgwas.Linux.precompiled.tar.gz
tar -zxvf dbgwas.Linux.precompiled.tar.gz
rm dbgwas.Linux.precompiled.tar.gz
cd ..

echo "${GREEN}Compiling MetaDBGWAS in src/${NC}"
cmake .
make -j 4

cd ..

echo "${GREEN}Installation complete !${NC} ${GREEN_BACKGROUND}Do not forget to install the R and C++ libraries.${NC}"
