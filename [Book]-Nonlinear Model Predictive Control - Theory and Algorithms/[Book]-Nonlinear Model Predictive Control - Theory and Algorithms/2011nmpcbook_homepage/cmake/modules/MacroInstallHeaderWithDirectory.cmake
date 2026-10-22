# Macro to install all headerfiles from a list with relative paths
# Destination: ${CMAKE_INSTALL_PREFIX}/include/${relative_path}/${headerfile}

macro( macro_install_headers_with_directory header_list )
	foreach( header ${${header_list}} )
		string( REGEX MATCH "(.*)[/\\]" DIR ${header} )
		install( FILES ${header} DESTINATION include/${DIR} )
	endforeach( header )
endmacro( macro_install_headers_with_directory )