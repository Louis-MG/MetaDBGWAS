#!/bin/sh

GREEN='\e[0;32m' # green color
NC='\e[0m' # No Color
GREEN_BACKGROUND='\e[42m' #self explanatory right?

#R lib
echo "${GREEN}Starting R libraries installation${NC}"
Rscript -e "install.packages(c('ape', 'phangorn'))"
Rscript -e "install.packages('https://raw.githubusercontent.com/sgearle/bugwas/master/build/bugwas_1.0.tar.gz', repos=NULL, type='source')"

#C++ lib
echo "${GREEN}Starting C++ libraries installation${NC}"
sudo apt install -y libgatbcore-dev libhdf5-dev libboost-all-dev libpstreams-dev zlib1g-dev g++ cmake git r-base-core

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
sed -i "51i#include <limits>" ./blight/robin_hood.h #temporary fix for REINDEER compilation
make

cd ..

echo "${GREEN}Compiling MetaDBGWAS in tool/src/${NC}"
cmake .
make -j 4

cd ..

chmod +x Metadbgwas/metadbgwas.sh
chmod +x Metadbgwas//tools/gemma/gemma.0.93b
chmod +x Metadbgwas/tools/phantomjs/phantomjs
echo "${GREEN}Installation complete !${NC} ${GREEN_BACKGROUND}Do not forget to install the R libraries !${NC}"
