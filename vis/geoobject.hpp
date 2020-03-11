#ifndef _GEOOBJECT_HPP_
#define _GEOOBJECT_HPP_

#include <GL/gl.h>
#include <GL/glext.h>
#include "psOpenGL.h"
#include <GL/glu.h>
#include <GL/glut.h>
#include "glcheck.hpp"

#include "common.hpp"
#include "vec3t.hpp"

typedef Vec3T<float> vec3f;

class GeoObject
{
public:  //GeoObject(const char* name) { assert(name!=NULL); strcpy(_name, name); }
    GeoObject(): _dummy(0) {;}
    virtual ~GeoObject() {;}  //const char* name() const { return _name; }
    virtual void render() = 0;
    virtual void idle() = 0;
    virtual void key(unsigned char) = 0;
    virtual void specialKey(unsigned char) = 0;
    virtual void mouse(int, int, int, int) = 0;
    virtual void motion(int, int) = 0;
    virtual vec3f centerPoint() = 0;
    virtual vec3f minPoint() = 0;
    virtual vec3f maxPoint() = 0;
    virtual void setPerspectiveParam(float, float, float, float) = 0;
  virtual int advanceRender() {return 0; }
protected:
    int _dummy;
};

#endif
