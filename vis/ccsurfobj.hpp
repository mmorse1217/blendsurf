#ifndef _CCSURFOBJ_HPP_
#define _CCSURFOBJ_HPP_

#include "ccsurf.hpp"
#include "ccsurfop.hpp"
#include "geoobject.hpp"

class CCSurfObj: public GeoObject {
public:
    enum {
        RENDER_SURF = 1,
        RENDER_FRAME = 2,
        RENDER_INFLBDRY = 4,
        RENDER_CPT = 8,
        RENDER_CCPTS = 16,
        RENDER_CPT_FACE = 32,
    };
    enum {
        SURF_NONE = 0,
        SURF_CUBEMAP = 1,
        SURF_CHECKBOARD = 2,
        SURF_GX = 3,
        SURF_GY = 4,
        SURF_GXX = 5,
        SURF_GXY = 6,
        SURF_GYY = 7,
        //...
        SURF_TTL = 8
    };
protected:
    CCSurf* _ccSurf;
    int _lvl;
    //-----------------
    int _renderctrl;
    int _surfctrl;
    int _activevert;
    GLuint _chkname; //check board texture name
    vector< CCRect<Point3> > _lmt;
    vector< CCRect<Point3> > _du;
    vector< CCRect<Point3> > _dv;
    vector< CCRect<Point3> > _nor;
    vector< CCRect<Point3> > _duu;
    vector< CCRect<Point3> > _duv;
    vector< CCRect<Point3> > _dvv;
  
public:
  CCSurfObj(CCSurf* surf, int lvl, int renderctrl = RENDER_CPT, 
	    int surfctrl = SURF_NONE): 
            GeoObject(), _ccSurf(surf), _lvl(lvl),
            _renderctrl(renderctrl), _surfctrl(surfctrl), _activevert(0) {
  }
    void render();
    void idle() {;}
    void key(unsigned char);
    void specialKey(unsigned char){};
    void mouse(int,int,int,int) {};
    void motion(int, int){};
    void setPerspectiveParam(float, float, float, float){};
    vec3f centerPoint() { Point3 c=_ccSurf->ctr(); return vec3f(c[0],c[1],c[2]); }
    vec3f minPoint() { Point3 a,b; _ccSurf->bbox(a,b); return vec3f(a[0],a[1],a[2]); }
    vec3f maxPoint() { Point3 a,b; _ccSurf->bbox(a,b); return vec3f(b[0],b[1],b[2]); }
    //others
    void setMaterial();
};

#endif
