cmake_minimum_required(VERSION 3.16) #verison minimum de cmake demandee, en genral la plus basse qui marche et par automatisme celle qui est utilisee
project(MetaDBGWAS) #nom du projet a utiliser pour le reste du script

set(CMAKE_CXX_STANDARD 11) # version C++

find_package(Boost REQUIRED COMPONENTS filesystem)

find_library(GATB_LIBRARY gatbcore)
find_package(PkgConfig)
pkg_search_module(hdf5 REQUIRED IMPORTED_TARGET hdf5)

add_executable(MetaDBGWAS build_dbg.cpp global.cpp Utils.cpp main.cpp) #fichier qui une fois compile sera l'executable, faut tous les mettre pour avoir l'executable
target_link_libraries(MetaDBGWAS PUBLIC Boost::filesystem PkgConfig::hdf5 ${GATB_LIBRARY})
