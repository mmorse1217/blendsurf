#ifdef __APPLE__
    #include <GLUT/glut.h>
#else
    #include <GL/glut.h>
#endif


#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#ifdef __CYGWIN__
#include <GL/glext.h>
#endif


#ifndef GL_EXT_texture_cube_map
# define GL_NORMAL_MAP_EXT                   0x8511
# define GL_REFLECTION_MAP_EXT               0x8512
# define GL_TEXTURE_CUBE_MAP_EXT             0x8513
# define GL_TEXTURE_BINDING_CUBE_MAP_EXT     0x8514
# define GL_TEXTURE_CUBE_MAP_POSITIVE_X_EXT  0x8515
# define GL_TEXTURE_CUBE_MAP_NEGATIVE_X_EXT  0x8516
# define GL_TEXTURE_CUBE_MAP_POSITIVE_Y_EXT  0x8517
# define GL_TEXTURE_CUBE_MAP_NEGATIVE_Y_EXT  0x8518
# define GL_TEXTURE_CUBE_MAP_POSITIVE_Z_EXT  0x8519
# define GL_TEXTURE_CUBE_MAP_NEGATIVE_Z_EXT  0x851A
# define GL_PROXY_TEXTURE_CUBE_MAP_EXT       0x851B
# define GL_MAX_CUBE_MAP_TEXTURE_SIZE_EXT    0x851C
#endif

static GLenum faceTarget[6] = {
    GL_TEXTURE_CUBE_MAP_POSITIVE_X_EXT,
    GL_TEXTURE_CUBE_MAP_NEGATIVE_X_EXT,
    GL_TEXTURE_CUBE_MAP_POSITIVE_Y_EXT,
    GL_TEXTURE_CUBE_MAP_NEGATIVE_Y_EXT,
    GL_TEXTURE_CUBE_MAP_POSITIVE_Z_EXT,
    GL_TEXTURE_CUBE_MAP_NEGATIVE_Z_EXT
};

void checkCubeMapExtension() {

    if (!glutExtensionSupported("GL_EXT_texture_cube_map")) {
        printf("cm_demo: Your OpenGL implementation does not support EXT_texture_cube_map.\n");
        printf("cm_demo: This program requires it to run.\n");
//        exit(1);
    }
}


unsigned char* makeWhite(int dim) { 
    unsigned char* imdata; 
    imdata = (unsigned char*)malloc(dim*dim*3);
    memset(imdata, 0xff, dim*dim*3);
    return imdata;
}

unsigned char* makeStripes(int dim, int nstripes) { 
    unsigned char* imdata;
    int i;
    int stripesize;
    stripesize = dim/nstripes; 
    imdata = (unsigned char*)malloc(dim*dim*3);
    memset(imdata, 0xff, dim*dim*3);
    //int stt = dim/16; int end = dim-stt;	for(i=stt; i<end; i+=2*stripesize) {	  for(int s=0; s<stripesize; s++)		 memset(imdata+(i+s)*dim*3+stt, 0x00, (end-stt)*3);	}		 
    for(i = 0; i <dim; i+= 2*stripesize) { 	  memset(imdata+i*3*dim, 0x00, stripesize*dim*3);	}
    return imdata;
}

void
makeCubeMap(int dim, int ns)
{
    int i;
  
    unsigned char *stripes, *white; 

    stripes = makeStripes(dim, ns);
    white = makeWhite(dim); 
    for (i=0; i<2; i++) 	 gluBuild2DMipmaps(faceTarget[i], 3, dim, dim, GL_RGB, GL_UNSIGNED_BYTE, white);
    for (i=2; i<4; i++) 	 gluBuild2DMipmaps(faceTarget[i], 3, dim, dim, GL_RGB, GL_UNSIGNED_BYTE, white);
    for (i=4; i<6; i++) 	 gluBuild2DMipmaps(faceTarget[i], 3, dim, dim, GL_RGB, GL_UNSIGNED_BYTE, stripes);
  
    //gluBuild2DMipmaps(faceTarget[i], 3, dim, dim, GL_RGB, GL_UNSIGNED_BYTE, (i<4?white:stripes));
    //gluBuild2DMipmaps(faceTarget[i], 3, dim, dim, GL_RGB, GL_UNSIGNED_BYTE, stripes);
    //gluBuild2DMipmaps(faceTarget[i], 3, dim, dim, GL_RGB, GL_UNSIGNED_BYTE, (i>=4?white:stripes));
  
    glTexParameteri(GL_TEXTURE_CUBE_MAP_EXT, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP_EXT, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
 
    glEnable(GL_TEXTURE_CUBE_MAP_EXT);

    glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_REFLECTION_MAP_EXT);
    glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_REFLECTION_MAP_EXT);
    glTexGeni(GL_R, GL_TEXTURE_GEN_MODE, GL_REFLECTION_MAP_EXT);
  
    glTexParameteri(GL_TEXTURE_CUBE_MAP_EXT, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP_EXT, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glEnable(GL_TEXTURE_GEN_S);
    glEnable(GL_TEXTURE_GEN_T);
    glEnable(GL_TEXTURE_GEN_R);
}

void cubeMapOn() { 
    glEnable(GL_TEXTURE_CUBE_MAP_EXT);
    glEnable(GL_TEXTURE_GEN_S);
    glEnable(GL_TEXTURE_GEN_T);
    glEnable(GL_TEXTURE_GEN_R);
}


void cubeMapOff() { 
    glDisable(GL_TEXTURE_CUBE_MAP_EXT);
    glDisable(GL_TEXTURE_GEN_S);
    glDisable(GL_TEXTURE_GEN_T);
    glDisable(GL_TEXTURE_GEN_R);
}
