#!/bin/bash

#lighter compilation
cd Lighter
make

cd ..

#bcalm compilation
cd bcalm
mkdir build;  cd build;  cmake ..;  make -j 8

cd ..

