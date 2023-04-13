#!/bin/bash

GREEN='\e[0;32m' # green color
NC='\e[0m' # No Color
GREEN_BACKGROUND='\e[42m' #self explanatory right?

#R lib
echo -e "${GREEN}Starting R libraries installation${NC}"
Rscript -e "install.packages(c('ape', 'phangorn'))"
Rscript -e "install.packages('https://raw.githubusercontent.com/sgearle/bugwas/master/build/bugwas_1.0.tar.gz', repos=NULL, type='source')"

#C++ lib
echo -e "${GREEN}Starting C++ libraries installation${NC}"
sudo apt install -y libgatbcore-dev libhdf5-dev libboost-all-dev libpstreams-dev zlib1g-dev g++ cmake git r-base-core build-essentials

#lighter compilation
echo -e "${GREEN}Starting Lighter installation${NC}"
cd Lighter
make

cd ..

echo -e "${GREEN}Starting bcalm2 installation${NC}"
#bcalm compilation
cd bcalm
mkdir build;  cd build;  cmake ..;  make -j 4

cd ../..

echo -e "${GREEN}Starting Bifrost installation${NC}"
cd Bifrost
cd bifrost && mkdir build && cd build
cmake ..
make

cd ../..

echo -e "${GREEN}Compiling MetaDBGWAS in tool/src/${NC}"
cmake .
make -j 4

cd ..

chmod +x ./metadbgwas.sh
chmod +x ./tools/gemma/gemma.0.93b
chmod +x ./tools/phantomjs/phantomjs
echo -e "${GREEN}Installation complete !${NC}"
