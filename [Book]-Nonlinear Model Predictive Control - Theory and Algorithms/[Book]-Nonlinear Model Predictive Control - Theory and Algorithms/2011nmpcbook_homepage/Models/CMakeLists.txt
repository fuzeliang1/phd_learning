# for #include "config.h and #include "libexport.h"
include_directories ( ${CMAKE_BINARY_DIR}/src )
# create config.h and libexport.h based on related *.h.cmake-files
configure_file ( ${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h @ONLY )
configure_file ( ${CMAKE_CURRENT_SOURCE_DIR}/libexport.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/libexport.h @ONLY )

#############################################################################
# Include and link libs
message ( STATUS "Searching for required libraries" )
PackageFindVersionCheck ( yane "REQUIRED" 1.0.0 yane_INCLUDE_DIRS yane_LIBRARIES yane_VERSION )
PackageFindVersionCheck ( sqpfortran "REQUIRED" 1.0.0 sqpfortran_INCLUDE_DIRS sqpfortran_LIBRARIES sqpfortran_VERSION )

#############################################################################
# create file-lists for shared library
message ( STATUS "Generating list of library source files" )
file ( GLOB lib_SRC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp )
list ( SORT lib_SRC )
message ( STATUS "List of target source files = ${lib_SRC}" )

message ( STATUS "Generating makefile for library ${LIB_NAME}" )
add_library ( ${LIB_NAME} SHARED ${lib_SRC} )
set_target_properties ( ${LIB_NAME} PROPERTIES VERSION ${LIB_VERSION} SOVERSION ${LIB_VERSION} )
target_link_libraries ( ${LIB_NAME} ${yane_LIBRARIES} ${sqpfortran_LIBRARIES} ${examplemodels2_LIBRARIES} )
message ( STATUS "Done" )
