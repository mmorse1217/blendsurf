cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
project(BlendSurf)

find_package(OpenMP REQUIRED)
find_package(lapack REQUIRED)

FUNCTION(PREPEND var )
   SET(listVar "")
   FOREACH(f ${${var}})
       LIST(APPEND listVar "${CMAKE_CURRENT_SOURCE_DIR}/${f}")
   ENDFOREACH(f)
   SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")

add_subdirectory(src)
add_subdirectory(include)

add_library(blend STATIC ${Blendsurf_SRC})
target_include_directories(blend PUBLIC ${BlendSurf_SOURCE_DIR}/include ${LAPACK_LIBRARIES})

