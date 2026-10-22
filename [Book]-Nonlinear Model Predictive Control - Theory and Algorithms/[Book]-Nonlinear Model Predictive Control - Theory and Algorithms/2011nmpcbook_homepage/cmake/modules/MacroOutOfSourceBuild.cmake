# Ensures that we do an out of source build

macro( macro_ensure_out_of_source_build msg)
	string( COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} insource )
	get_filename_component( parent_dir ${CMAKE_SOURCE_DIR} PATH )
	string( COMPARE EQUAL ${parent_dir} ${CMAKE_SOURCE_DIR} insourcedir )
	if( insource OR insourcedir )
		message( FATAL_ERROR "${MSG}" )
	endif( insource OR insourcedir )
endmacro( macro_ensure_out_of_source_build )

# MACRO(MACRO_ENSURE_OUT_OF_SOURCE_BUILD MSG)
#     STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}"
# "${CMAKE_BINARY_DIR}" insource)
#     GET_FILENAME_COMPONENT(PARENTDIR ${CMAKE_SOURCE_DIR} PATH)
#     STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}"
# "${PARENTDIR}" insourcesubdir)
#     IF(insource OR insourcesubdir)
#         MESSAGE(FATAL_ERROR "${MSG}")
#     ENDIF(insource OR insourcesubdir)
# ENDMACRO(MACRO_ENSURE_OUT_OF_SOURCE_BUILD)
