add_executable(MetaDBGWAS Blast.cpp PhenoCounter.cpp input_files_gen.cpp Utils.cpp generate_output.cpp global.cpp main.cpp statistical_test.cpp whereami.c)
target_link_libraries(MetaDBGWAS PUBLIC Boost::serialization Boost::filesystem Boost::graph GatbCore)
