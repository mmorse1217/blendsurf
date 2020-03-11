
set(Blendsurf_SRC
    bdsurf.cpp  
    bdulib.cpp  
    ccsubmatlib.cpp 
    ccsurf.cpp  
    ccsurfop.cpp  
    gpmesh.cpp 
    evalfromw.cpp 
    vecmatop.cpp    
    quaternion.cpp  
    svdrep.cpp 
    evalfromw.cpp 
    splinebasis.c
    utils.cpp
    )

PREPEND(Blendsurf_SRC)
set(Blendsurf_SRC ${Blendsurf_SRC}  PARENT_SCOPE)
