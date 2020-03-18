# list of required source files
set(Blendsurf_Visualizer_SRC 
    viewer3d.cpp  
    bdsurfobj.cpp 
    reflmap.cpp 
    ccsurfobj.cpp  
    ballviewer.cpp 
    viewer.cpp 
    psOpenGL.cpp  
    camera.cpp  
    arcball.cpp 
    texture.cpp 
    )

PREPEND(Blendsurf_Visualizer_SRC)
set(Blendsurf_Visualizer_SRC ${Blendsurf_Visualizer_SRC}  PARENT_SCOPE)


# find packages  for rendering
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

# set up an intermediate library to link
add_library(blendsurf_render STATIC ${Blendsurf_Visualizer_SRC})
add_library(Blendsurf::render ALIAS blendsurf_render)
set_property(TARGET blendsurf_render PROPERTY EXPORT_NAME Blendsurf::render)

target_link_libraries(blendsurf_render OpenGL::GL)
target_link_libraries(blendsurf_render OpenGL::GLU)
target_link_libraries(blendsurf_render GLUT::GLUT)
target_link_libraries(blendsurf_render Blendsurf::blendsurf)

# build the executable and link
add_executable(blendsurf_renderer  vistt0_sampleout.cpp)
target_link_libraries(blendsurf_renderer Blendsurf::render)

#note: several files in visoption3d need to be present in order for renderer to
#succeed. See the README.md
