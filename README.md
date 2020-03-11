Blendsurf is a C++ library to represent surfaces of arbitrary smoothness (C^\infty and C^k) written by Lexing Ying and Elif Tosun. Minor additions added by Matt Morse. Currently NOT  maintained. It is made available for use as is in other projects.

Requires:
* OpenMP
* BLAS/LAPACK
* CMake 3.1

Compilation tested on OSX 10.10 and Cent OS 7.3 with gcc-4.8 and icc/icpc 17.0.1.

To compile the library:

    mkdir build
    cd build
    cmake ..
    make 


Legacy rendering code is available in vis/. This is not maintained and does not compile with the current CMakeLists.txt files. The original makefile and makefile.in are included but not maintained.


On the off chance someone gets the viewer running again:
The options file for the visualization code is visoption3d, with some comments explaining each value.  In the legacy viewer, the following keys are active:

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
	    
