add_executable(MetaDBGWAS build_dbg.cpp global.cpp Utils.cpp main.cpp) #fichier qui une fois compile sera l'executable, faut tous les mettre pour avoir l'executable
target_link_libraries(MetaDBGWAS PUBLIC Boost::filesystem GatbCore)
