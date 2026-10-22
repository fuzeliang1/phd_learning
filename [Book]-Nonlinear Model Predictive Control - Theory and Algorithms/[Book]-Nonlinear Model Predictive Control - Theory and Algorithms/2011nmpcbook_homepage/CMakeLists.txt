PROJECT ( application_2011_nmpcbook )

set ( LIB_NAME Models )
set ( LIB_VERSION 1.0.0 )

cmake_minimum_required ( VERSION 2.6 )

#############################################################################
set ( CMAKE_BUILD_TYPE Release )
message ( STATUS "Build Type: ${CMAKE_BUILD_TYPE} | C-Flags: ${CMAKE_C_FLAGS_RELEASE} | CXX-Flags: ${CMAKE_CXX_FLAGS_RELEASE} \n" )

#############################################################################
# include some macros
set ( CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules" ${CMAKE_MODULE_PATH} )

include ( MacroInstallHeaderWithDirectory )
include ( MacroListOperations )
include ( MacroOutOfSourceBuild )
include ( MacroPackageFind )

macro_ensure_out_of_source_build ( "${PROJECT_NAME} requires an out of source build. In-source builds are not allowed. CMake would overwrite the makefiles distributed with this library. Please create a directory and run cmake from there, passing the path to this source directory as the last argument. This process created the file `CMakeCache.txt' and the directory `CMakeFiles'. Please delete them." )

#############################################################################
include_directories ( ${LIB_NAME} )
add_subdirectory ( ${LIB_NAME} )


#############################################################################
include_directories ( Examples )
add_subdirectory ( Examples )
