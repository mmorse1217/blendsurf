cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
project(Blendsurf)

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


target_include_directories(blendsurf PUBLIC ${Blendsurf_SOURCE_DIR}/include )
target_link_libraries(blendsurf ${LAPACK_LIBRARIES})
target_link_libraries(blendsurf stdc++)

# compile the renderer
if (COMPILE_RENDERER)
    add_subdirectory(vis)
endif()

