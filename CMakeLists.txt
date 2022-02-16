cmake_minimum_required(VERSION 3.16) #verison minimum de cmake demandee, en genral la plus basse qui marche et par automatisme celle qui est utilisee
project(MetaDBGWAS) #nom du projet a utiliser pour le reste du script

set(CMAKE_CXX_STANDARD 11) # version C++

find_package(Boost REQUIRED COMPONENTS filesystem)

# Use manual find_package script until Gatb supports cmake config files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
find_package(GatbCore REQUIRED)

add_executable(MetaDBGWAS build_dbg.cpp global.cpp Utils.cpp) #fichier qui une fois compile sera l'executable, faut tous les mettre pour avoir l'executable
target_link_libraries(MetaDBGWAS PUBLIC Boost::filesystem GatbCore)
