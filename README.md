![Docker Image CI](https://github.com/mmorse1217/blendsurf/workflows/Docker%20Image%20CI/badge.svg)

Blendsurf is a C++ library to represent surfaces of arbitrary smoothness (C^\infty and C^k) written by Lexing Ying and Elif Tosun. Minor additions added by Matt Morse. Currently NOT  maintained. It is made available for use as is in other projects.

Requires:
* OpenMP
* BLAS/LAPACK
* CMake 3.1
    
The renderer requires:
* OpenGL version that supports `EXT_texture_cube_map`
    
Library compilation tested on OSX 10.10 and Cent OS 7.3 with gcc-4.8 and icc/icpc 17.0.1. The renderer has been tested with gcc-4.8 on OSX 10.10.

To compile the library:

    mkdir build
    cd build
    CC=<your-c-compiler> CXX=<your-cxx-compiler> cmake ..
    make 

To compile the renderer, use

    CC=<your-c-compiler> CXX=<your-cxx-compiler> cmake -DCOMPILE_RENDERER=True ..

Legacy rendering code is available in `vis/`. This is not maintained. It compiles and runs as expected with the current CMakeLists.txt's, but it requires `EXT_texture_cube_map` in your linked OpenGL implementation. This is deprecated and current versions of OpenGL do not support it. Moreover, the options file `vis/visoptions3d` requires `-bdsurf_submatlibfile` and `-bdsurf_bdulibfile` to point to the location of `ccsubmatall.dat` and `bdsurf_U_ONE.dat` respectively, which are downloaded from the link below. Legacy makefiles are included for future reference.

On the off chance someone gets the viewer running again:
The options file for the visualization code is `visoption3d`, with some comments explaining each value.  In the legacy viewer, the following keys are active:

-   P    = display/remove control mesh
-   R    = render surface
-   S    = surface control, switches between different views
	    reflection map, checkerboard map, higher order derivatives
-   V    = switch from one vertex to another for surface control
-   F    = wireframe
-   I    = influence boundary
-   C    = center
-   home = return to original view
-   W    = write ppm file
-   Q    = exit

Legacy test cases, precomputed data required for some options, and scripts are available at https://drive.google.com/open?id=13Gpjq5fFwwzCQ4XnaRQQOliNbOwTLaNP,  but I'm not exactly sure what most of the files are for.
