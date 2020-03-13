//#include "CCcommon.h"
//#include "compat.h"

#include <GL/gl.h>
#include "psOpenGL.h"
#include <GL/glu.h>
#ifdef __APPLE__
    #include <GLUT/glut.h>
#else
    #include <GL/glut.h>
#endif

#include "glcheck.hpp"
#include "viewer3d.hpp"
#include "geoobject.hpp"

//static void spositionCamera(Camera* camera, GeoObject* object, int* vp);

Viewer3D::Viewer3D(char* t, int w, int h) : Viewer(t, w, h) {
    _obj = NULL;
    _camera = new Camera();
    glCheck();
}

Viewer3D::~Viewer3D() { delete _camera; }

void Viewer3D::setObject(GeoObject* object) 
{
    _obj = object;
    _needPosition = true;    //int vp[4] = {0, 0, getWidth(), getHeight()};  //spositionCamera(_camera, object, vp);
}

void Viewer3D::setCamera(Camera* cam) { 
  delete _camera;  
  _camera = cam;
  _needPosition = false; 
  cerr << "setting _need to false" << endl; 
}


vec3f Viewer3D::minPoint() const
{
    assert(_obj!=NULL);
    return _obj->minPoint();  //for(int i=0; i<_objvec.size(); i++)	 _mi = min( _mi, _objvec[i]->minPoint() );
}

vec3f Viewer3D::maxPoint() const
{
    assert(_obj!=NULL);
    return _obj->maxPoint();
}

vec3f Viewer3D::centerPoint() const
{
    assert(_obj!=NULL);
    return _obj->centerPoint();
}


//----------------------------------------------------------------------
void Viewer3D::display() 
{
    //LEXING: if(_needPosition && (getObject() != 0))
    if(_needPosition && _obj!=NULL)    {
        //int vp[4] = {0, 0, getWidth(), getHeight()};
        //spositionCamera(_camera, object, vp);
        positionCamera();
        _needPosition = false;
    }
    glClearColor( 1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glCheck();  assert(0);
    _obj->render();
    glCheck();
    glutSwapBuffers();
    glCheck();
}

void Viewer3D::idle()
{
    _obj->idle();
    glCheck();
}

void Viewer3D::reshape(int w, int h) 
{
    Viewer::reshape(w, h);
  
    int vp[4] = {0, 0, w, h};
    glViewport(0, 0, w, h);
    _camera->setViewport(vp);
    _camera->computeModelview();
    _camera->computeProjection();
    _camera->loadMatrices();
    glutPostRedisplay();
    glCheck();
}

void Viewer3D::positionCamera() 
{
    int vp[4] = { 0, 0,getWidth(), getHeight() }; 

    glCheck();
    Camera* camera = getCamera();
    vec3f _mi = minPoint();
    vec3f _ma = maxPoint();

    glCheck();
    camera->reset();

    glCheck();
    camera->setViewport(vp);
    vec3f mid = 0.5f * (_mi + _ma);
  
    float f = max(_ma.x()-_mi.x(),_ma.y()-_mi.y()) / 2.0f / tan(0.5f * camera->fovy() * M_PI / 180.0f);
    vec3f trans(mid.x(), mid.y(), _ma.z() + 2 * f);
    glCheck();
  
    camera->setPerspectiveParams(camera->fovy(), 
                                 0.1f * (trans.z() - _mi.z()), 
                                 10.0f * (trans.z() - _ma.z()));
    glCheck();

    camera->translate(-trans);
    camera->computeModelview();
    camera->computeProjection();
    camera->loadMatrices();  
    glCheck();
}
