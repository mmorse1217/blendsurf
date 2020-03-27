cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
project(Blendsurf)
set(Blendsurf_VERSION 0.1)

FUNCTION(PREPEND var )
   SET(listVar "")
   FOREACH(f ${${var}})
       LIST(APPEND listVar "${CMAKE_CURRENT_SOURCE_DIR}/${f}")
   ENDFOREACH(f)
   SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")

find_package(OpenMP REQUIRED)
find_package(LAPACK REQUIRED)

# include source code and headers
add_subdirectory(src)
add_subdirectory(include)

# add the library as a target with src and include as dependencies
add_library(blendsurf STATIC ${Blendsurf_SRC} ${Blendsurf_INC})
add_library(Blendsurf::blendsurf ALIAS blendsurf)
set_property(TARGET blendsurf PROPERTY EXPORT_NAME Blendsurf::blendsurf)

set(Blendsurf_INCLUDE_DEST "include/blendsurf-${Blendsurf_VERSION}")
set(Blendsurf_LIB_DEST "lib/blendsurf-${Blendsurf_VERSION}")


#target_include_directories(blendsurf PUBLIC ${Blendsurf_SOURCE_DIR}/include )

target_include_directories(blendsurf PUBLIC
   $<BUILD_INTERFACE:${Blendsurf_SOURCE_DIR}/include> # for headers when building
   $<BUILD_INTERFACE:${Blendsurf_BINARY_DIR}/include> # for headers when building
   #$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # for config_impl.hpp when building
   $<INSTALL_INTERFACE:include> # for client in install mode
   #$<INSTALL_INTERFACE:${Blendsurf_LIB_DEST}> # for config_impl.hpp in install mode)
   )
target_link_libraries(blendsurf ${LAPACK_LIBRARIES})
target_link_libraries(blendsurf stdc++)

# compile the renderer
if (COMPILE_RENDERER)
    add_subdirectory(vis)
    # TODO add install support for renderer...?
endif()

# Install blendsurf in CMAKE_INSTALL_PREFIX (defaults to /usr/local on linux)
install(
    # targets to install
    TARGETS blendsurf 
    # name of the CMake "export group" containing the targets we want to install
    EXPORT blendsurfTargets
    # Dynamic, static library and include destination locations after running
    # "make install"
    LIBRARY DESTINATION ${Blendsurf_LIB_DEST}
    ARCHIVE DESTINATION ${Blendsurf_LIB_DEST} 
    INCLUDES DESTINATION ${Blendsurf_INCLUDE_DEST}
    )
install(
    # export group we want to save (matches name defined above containing the
    # install targets)
    EXPORT blendsurfTargets
    # CMake file in which to store the export group info
    FILE  blendsurfTargets.cmake
    # namespace prepends all targets in the export group
    NAMESPACE Blendsurf::
    # where to place the resulting file (here we're putting it with the library)
    DESTINATION ${Blendsurf_LIB_DEST}
    )

# Copies the files in ${Blendsurf_INC} to ${Blendsurf_INCLUDE_DEST}
install(FILES ${Blendsurf_INC} DESTINATION ${Blendsurf_INCLUDE_DEST})

# Write a "version file" in case some wants to only load a particular version of
# blendsurf @BUG this fails at the moment.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    blendsurfConfigVersion.cmake
    VERSION ${Blendsurf_VERSION}
    COMPATIBILITY AnyNewerVersion
    )

# Copies the resulting CMake config files to the installed library directory
install(FILES 
    "cmake/blendsurf-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/blendsurfConfigVersion.cmake"
    DESTINATION ${Blendsurf_LIB_DEST})

