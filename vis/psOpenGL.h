#ifndef __psOpenGL_h__
#define __psOpenGL_h__

#include "common.hpp"


#ifdef __cplusplus
extern "C" 
{
#endif
  
    void psOpenGL_open( const char* filename );
    void psOpenGL_close( void );

    extern int psOpenGL_fill_gray;
#ifdef __cplusplus
}
#endif

#ifdef PSOPENGL
#ifndef PSOPENGL_CODE



#define glEnd   psOpenGL_End

void psOpenGL_End( void );


#endif /* !PSOPENGL_CODE */
#endif /* PSOPENGL */


#endif /* __psOpenGL_h__ */
