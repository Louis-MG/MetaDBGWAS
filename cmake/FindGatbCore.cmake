# Try to find an install of Gatb-core and generate a CMake IMPORTED target from it ("GatbCore").
# This only exists because Gatb-core does not yet have cmake config files.
# When Gatb-core does have cmake config files, a simple find_package will be sufficient and this can be removed.

include(FindPackageHandleStandardArgs)

# Find header path and library path
find_path(GATB_CORE_INCLUDE_DIR gatb/gatb_core.hpp)
find_library(GATB_CORE_LIBRARY gatbcore)

# Standard checks and message
find_package_handle_standard_args(
	GatbCore
	DEFAULT_MSG
	GATB_CORE_LIBRARY GATB_CORE_INCLUDE_DIR 
)
mark_as_advanced(GATB_CORE_INCLUDE_DIR GATB_CORE_LIBRARY) # For Cmake GUI

# Generate target. Settings are declared as INTERFACE as we declare an external library, we do not build one ourselves.
add_library(GatbCore UNKNOWN IMPORTED)
set_property(TARGET GatbCore PROPERTY IMPORTED_LOCATION "${GATB_CORE_LIBRARY}")
target_include_directories(GatbCore INTERFACE "${GATB_CORE_INCLUDE_DIR}")

# Gatb_core links against hdf5 which has no cmake files but has working pkg-config
find_package(PkgConfig)
pkg_search_module(hdf5 REQUIRED IMPORTED_TARGET hdf5)
target_link_libraries(GatbCore INTERFACE PkgConfig::hdf5)

