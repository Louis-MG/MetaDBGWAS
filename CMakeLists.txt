cmake_minimum_required(VERSION 3.10) #verison minimum de cmake demandee, en genral la plus basse qui marche et par automatisme celle qui est utilisee
project(MetaDBGWAS) #nom du projet a utiliser pour le reste du script

set(CMAKE_CXX_STANDARD 20) # version C++

find_package(Boost REQUIRED COMPONENTS filesystem serialization)

# Use manual find_package script until Gatb supports cmake config files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
find_package(GatbCore REQUIRED)

add_subdirectory(src)
