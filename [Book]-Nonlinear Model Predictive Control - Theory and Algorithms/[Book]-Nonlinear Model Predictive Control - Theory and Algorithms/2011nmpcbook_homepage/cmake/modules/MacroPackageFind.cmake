# Macros for finding Packages
# Both macros includes directories to headerfiles.


# Macro without checking versions
# input-var:
#	package - Package-Name
#	option  - String containing REQUIRED and/or QUIET, use "" for no Option
# output-var:
#	package_include - contains include-path
#	package_libs    - contains library-name with path
macro( PackageFind package option package_include package_libs )

	set( find_options "" )
	if( ${option} MATCHES ".*REQUIRED.*" )
		set( find_options ${find_options} REQUIRED )
	endif( ${option} MATCHES ".*REQUIRED.*" )
	if( ${option} MATCHES ".*QUIET.*" )
		set( find_options ${find_options} QUIET )
	endif( ${option} MATCHES ".*QUIET.*" )

	find_package( ${package} ${find_options} )
	
	if( ${package}_FOUND )
		message( STATUS "${package}_INCLUDE_DIRS = ${${package}_INCLUDE_DIRS}" )
		message( STATUS "${package}_LIBRARIES = ${${package}_LIBRARIES}" )
	endif( ${package}_FOUND )

	set( ${package_include} ${${package}_INCLUDE_DIRS} )
	set( ${package_libs} ${${package}_LIBRARIES} )

	include_directories( ${${package}_VERSION} )

endmacro( PackageFind package option package_include package_libs )


# Macro without checking versions
# input-var:
#	package     - Package-Name
#	req_version - Required Package-Version (not as string)
#	option      - String containing REQUIRED and/or QUIET, use "" for no Option
# output-var:
#	package_include - contains include-path of installed package
#	package_libs    - contains library-name with path of installed packege
#	package_version - contains version of installed package
macro( PackageFindVersionCheck package option req_version package_include package_libs package_version )

	set( find_options "" )
	if( ${option} MATCHES ".*REQUIRED.*" )
		set( find_options ${find_options} REQUIRED )
	endif( ${option} MATCHES ".*REQUIRED.*" )
	if( ${option} MATCHES ".*QUIET.*" )
		set( find_options ${find_options} QUIET )
	endif( ${option} MATCHES ".*QUIET.*" )

	find_package( ${package} ${find_options} )

	if( ${package}_FOUND )
		message( STATUS "${package}_INCLUDE_DIRS = ${${package}_INCLUDE_DIRS}" )
		message( STATUS "${package}_LIBRARIES = ${${package}_LIBRARIES}" )
		message( STATUS "${package}_VERSION = ${${package}_VERSION}" )
		if( ${${package}_VERSION} VERSION_LESS ${req_version} )
			message( FATAL_ERROR "installed version not compatible, minimum required version: ${req_version}" )
		endif( ${${package}_VERSION} VERSION_LESS ${req_version} )
	endif( ${package}_FOUND )

	set( ${package_include} ${${package}_INCLUDE_DIRS} )
	set( ${package_libs} ${${package}_LIBRARIES} )
	set( ${package_version} ${${package}_VERSION} )

	include_directories( ${${package}_VERSION} )

endmacro( PackageFindVersionCheck package option req_version package_include package_libs package_version ) 

