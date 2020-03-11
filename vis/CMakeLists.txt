# This is just template code to pass the source code variables to the parent
# scope, in case the renderer needs to work in the future.
set(Blendsurf_Visualizer_SRC 
    vis/viewer3d.cpp  
    vis/bdsurfobj.cpp 
    vis/reflmap.cpp 
    vis/ccsurfobj.cpp  
    vis/ballviewer.cpp 
    vis/viewer.cpp 
    vis/psOpenGL.cpp  
    vis/camera.cpp  
    vis/arcball.cpp 
    vis/texture.cpp 
    )
PREPEND(Blendsurf_SRC)
set(Blendsurf_Visualizer_SRC ${Blendsurf_Visualizer_SRC}  PARENT_SCOPE)

#find_library(GLUT_LIB glut REQUIRED)
#find_package(png REQUIRED)
#include_directories($ENV{OPENGL_DIR})
#include_directories(${LAPACK_LIBRARIES})
#find_library(OPENGL_LIB OpenGL REQUIRED)
#add_executable(bdrender bdrender.cpp)
#add_executable(vistt0_sampleout vistt0_sampleout.cpp)
#target_link_libraries(vistt0_sampleout blend ${LAPACK_LIBRARIES} png ${GLUT_LIB}
#    ${OPENGL_LIB} )
#${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )
