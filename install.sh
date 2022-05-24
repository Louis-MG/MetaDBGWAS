#!/bin/bash

GREEN='\e[0;32m' # green color
NC='\e[0m' # No Color
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

echo "${GREEN}Compiling MetaDBGWAS in tool/src/${NC}"
cmake .
make -j 4

cd ..

chmod +x metadbgwas.sh
echo "${GREEN}Installation complete !${NC} ${GREEN_BACKGROUND}Do not forget to install the R libraries !${NC}"
