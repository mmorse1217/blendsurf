
set(Blendsurf_INC
    bdsurf.hpp  
    bdulib.hpp  
    blas_lapack.h
    ccsubmatlib.hpp 
    ccsurf.hpp  
    ccsurfop.hpp  
    common.hpp  
    cvect.hpp  
    evalfromw.hpp 
    gpmesh.hpp 
    hmatrix.hpp 
    mat2t.hpp 
    mat3t.hpp 
    nummatrix.hpp 
    numvector.hpp 
    quat.hpp  
    quaternion.hpp  
    svdrep.hpp 
    utils.hpp
    vec2t.hpp
    vec3t.hpp
    vecmatop.hpp
    )
PREPEND(Blendsurf_INC)
set(Blendsurf_INC ${Blendsurf_INC}  PARENT_SCOPE)
