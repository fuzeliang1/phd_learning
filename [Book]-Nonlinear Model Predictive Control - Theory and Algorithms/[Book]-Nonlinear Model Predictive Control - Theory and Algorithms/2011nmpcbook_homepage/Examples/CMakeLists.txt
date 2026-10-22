#############################################################################
# Include and link libs
message ( STATUS "Searching for required libraries" )
PackageFindVersionCheck ( yane "REQUIRED" 1.0.0 yane_INCLUDE_DIRS yane_LIBRARIES yane_VERSION )
PackageFindVersionCheck ( sqpfortran "REQUIRED" 1.0.0 sqpfortran_INCLUDE_DIRS sqpfortran_LIBRARIES sqpfortran_VERSION )
#PackageFindVersionCheck ( sqpnagc3 "REQUIRED" 1.0.0 sqpnagc_INCLUDE_DIRS sqpnagc_LIBRARIES sqpnagc_VERSION )
#PackageFind ( nagc "REQUIRED" nagc_INCLUDE_DIRS nagc_LIBRARIES )
#PackageFindVersionCheck ( yaneipopt "REQUIRED" 1.0.0 yaneipopt_INCLUDE_DIRS yaneipopt_LIBRARIES yaneipopt_VERSION )
#PackageFind ( ipopt "REQUIRED" ipopt_INCLUDE_DIRS ipopt_LIBRARIES )
#PackageFindVersionCheck ( examplemodels2 "REQUIRED" 2.0.0 examplemodels2_INCLUDE_DIRS examplemodels2_LIBRARIES examplemodels2_VERSION )

#############################################################################
# include libraries
#include_directories ( ${yane_INCLUDE_DIRS} ${sqpfortran_INCLUDE_DIRS} ${nagc_INCLUDE_DIRS} ${sqpnagc_INCLUDE_DIRS} ${ipopt_INCLUDE_DIRS} ${yaneipopt_INCLUDE_DIRS} ${examplemodels2_INCLUDE_DIRS} )
include_directories ( ${yane_INCLUDE_DIRS} ${sqpfortran_INCLUDE_DIRS} )

#############################################################################
# Create bin
message ( STATUS "Generating list of target source files" )
file ( GLOB target_SRC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp )
message ( STATUS "Done" )
message ( STATUS "List of target source files = ${target_SRC}" )

foreach ( target ${target_SRC} )
    message ( STATUS "Processing file: ${target}" )
    string ( REGEX REPLACE ".cpp" "" target ${target} )
    add_executable ( ${target} ${target}.cpp )
    target_link_libraries ( ${target}
        ${yane_LIBRARIES}
        ${sqpfortran_LIBRARIES}
        #${sqpnagc3_LIBRARIES}
        #${nagc_LIBRARIES}
        #${yaneipopt_LIBRARIES}
        #${ipopt_LIBRARIES}
        #${examplemodels2_LIBRARIES}
        ${LIB_NAME}
    )
    message ( STATUS "Generating makefile for binary" )
    install ( TARGETS ${target} RUNTIME DESTINATION bin )
    message ( STATUS "Done" )
endforeach ( target )