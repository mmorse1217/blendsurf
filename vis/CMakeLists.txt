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

# find packages  for rendering
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(png REQUIRED)

# set up an intermediate library to link
add_library(blendsurf_render STATIC ${Blendsurf_Visualizer_SRC})
target_link_libraries(blendsurf_render OpenGL::GL)
target_link_libraries(blendsurf_render GLUT::GLUT)
target_link_libraries(blendsurf_render blendsurf)

# build the executable and link
add_executable(blendsurf_renderer  vistt0_sampleout.cpp)
target_link_libraries(blendsurf_renderer blendsurf_render)

#note: several files in visoption3d need to be present in order for renderer to
#succeed. See the README.md
