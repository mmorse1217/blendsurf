#get_filename_component(SELF_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
#include(${SELF_DIR}/blendsurf.cmake)

include(CMakeFindDependencyMacro)
# Capturing values from configure (optional)
#set(my-config-var @my-config-var@)

# Same syntax as find_package
find_dependency(OpenMP REQUIRED)
find_dependency(LAPACK REQUIRED)

# Any extra setup

# Add the targets file
include("${CMAKE_CURRENT_LIST_DIR}/blendsurfTargets.cmake")
