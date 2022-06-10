apt-get update
apt install -y libgatbcore-dev libhdf5-dev libboost-all-dev libpstreams-dev zlib1g-dev g++ cmake git r-base-core
Rscript -e "install.packages(c('ape', 'phangorn'))"
Rscript -e "install.packages('https://raw.githubusercontent.com/sgearle/bugwas/master/build/bugwas_1.0.tar.gz', repos=NULL, type='source')"
git clone --recursive https://github.com/Louis-MG/Metadbgwas.git
cd Metadbgwas
chmod +x metadbgwas.sh
chmod +x ./tools/phantomjs/phantomjs
chmod +x ./tools/gemma/gemma.0.93b
sed -i "51i#include <limits>" ./REINDEER/blight/robin_hood.h #temporary fix for REINDEER compilation
sh install.sh
