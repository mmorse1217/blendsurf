#!/bin/bash

# code to execute during CI build
mkdir -p /blendsurf/build/ 
cd /blendsurf/build/
cmake -DCMAKE_MODULE_PATH=/usr/share/cmake-3.10/Modules/ -DCOMPILE_RENDERER=True ..  
make 
