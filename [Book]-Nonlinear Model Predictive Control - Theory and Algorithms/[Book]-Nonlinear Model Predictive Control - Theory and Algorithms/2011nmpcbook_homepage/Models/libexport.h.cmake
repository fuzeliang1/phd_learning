
#ifndef _LIBEXPORT_H_
#define _LIBEXPORT_H_

/* Cmake will define ${LIB_NAME}_EXPORTS on Windows when it
configures to build a shared library with target-name ${LIB_NAME}.
If you are going to use another build system on windows or create 
the visual studio projects by hand you need to define 
MyLibrary_EXPORTS when building a DLL on windows.
*/


#if defined (_WIN32)
	#if defined( @LIB_NAME@_EXPORTS )
		#define LIB_EXPORT __declspec(dllexport)
	#else
		#define LIB_EXPORT __declspec(dllimport)
	#endif /* COMPILING_DLL */
#else /* defined (_WIN32) */
	#define LIB_EXPORT
#endif

#endif /* _LIBEXPORT_H_ */

