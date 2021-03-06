### CMake file for cyan_imageviewer 
###
###     Author : glaconbleu@codices.org
###
###############################################################################

cmake_minimum_required (VERSION 2.6)

project (cyan_transforms)
set (cyan_transforms_VERSION_MAJOR 0)
set (cyan_transforms_VERSION_MINOR 1)
set (cyan_transforms_VERSION_PATCH 1)
set (cyan_transforms_VERSION_STRING ${cyan_transforms_VERSION_MAJOR}.${cyan_transforms_VERSION_MINOR}.${cyan_transforms_VERSION_PATCH} )

 set( CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /lib /usr/lib /usr/local/lib )
 set( CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} /usr/include /usr/local/include )
set( CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake" ) 

set ( INSTALL_BIN_PATH bin )
set ( INSTALL_LIB_PATH lib )
set ( INSTALL_AR_PATH  lib )
set ( INSTALL_HEADERS_PATH include )

##set( CMAKE_INSTALL_RPATH "/usr/local/lib/")
##set( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath -Wl,/usr/local/lib") 

# Make relative paths absolute (needed later on)

foreach(p LIB BIN HEADERS AR)
        set(var INSTALL_${p}_PATH)
        if(NOT IS_ABSOLUTE "${${var}}")
                set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
        endif()
endforeach()

### Options ###

### Set include directories 
        
include_directories( "${CMAKE_CURRENT_SOURCE_DIR}/include/"
	                 "${CMAKE_CURRENT_BINARY_DIR}" )
		 if(PNG_FOUND)
		 include_directories("${PNG_INCLUDE_DIRS}")
		endif()
### File lists

file( GLOB src "src/*.c"  )

### Targets

add_executable( "cyan_transforms" ${src} ${hdr} ) 

### Linking

target_link_libraries( "cyan_transforms" "m" )
##target_link_libraries( "cyan_imageviewer" ${PNG_LIBRARIES} )
target_link_libraries( "cyan_transforms" "cyan" )

### Install 
install( TARGETS "cyan_transforms"
	    RUNTIME DESTINATION ${INSTALL_BIN_PATH}  COMPONENT main
	    LIBRARY DESTINATION ${INSTALL_LIB_PATH}  COMPONENT main
        ARCHIVE DESTINATION ${INSTALL_AR_PATH}   COMPONENT main )

