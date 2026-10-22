# Some Macros for Listoperations

# macro to add directory-name to every entry of a list
macro( macro_add_dir_to_list list dir output_list )
	foreach( header ${${list}} )
		list( APPEND ${output_list} ${dir}/${header} )
	endforeach( header )
endmacro( macro_add_dir_to_list )

# macro to remove all components in remove_list from list
macro( macro_remove_headers_from_list list remove_list )
	foreach( header ${${remove_list}} )
		list( REMOVE_ITEM ${list} ${header} )
	endforeach( header )
	list( SORT ${list} )
endmacro( macro_remove_headers_from_list )