// -*- Mode: c++ -*-
// $Id$
// $Source$

/* Subdivide V2.0
   Copyright (C) 2000 Henning Biermann, Denis Zorin, NYU

   This file is part of Subdivide.

   Subdivide is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   Subdivide is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   for more details.

   You should have received a copy of the GNU General Public License
   along with Subdivide; see the file COPYING.  If not, write to the Free
   Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
   02111-1307, USA.  */

//#include "compat.h"
#include "camera.hpp"

// storage for current matrices
GLdouble Camera::_proj[16];
GLdouble Camera::_view[16];
GLint Camera::_vp[4];

// ******************** Camera memeber functions ***********************

void Camera::identityModel() {
    _ballModel.setIdentity();
}

Camera::Camera(const Camera & c) :
        _projectionMatrix(c._projectionMatrix),
        _modelviewMatrix(c._modelviewMatrix),
        _trans(c._trans),
        _ballModel(c._ballModel),
        _znear(c._znear), 
        _zfar(c._zfar),
        _fovy(c._fovy), 
        _aspect(c._aspect),
        _wincenterx(c._wincenterx),
        _wincentery(c._wincentery),
        _winscale(c._winscale)
{

    _viewport[0] = c._viewport[0];
    _viewport[1] = c._viewport[1];
    _viewport[2] = c._viewport[2];
    _viewport[3] = c._viewport[3];
}

const Camera & Camera::operator=(const Camera & c) {
    _trans = c._trans;
    _fovy = c._fovy; _znear = c._znear; _zfar = c._zfar;
    _wincenterx =  c._wincenterx; _wincentery =  c._wincentery;
    _winscale = c._winscale;
    _ballModel = c._ballModel;
    _modelviewMatrix = c._modelviewMatrix;
    _projectionMatrix = c._projectionMatrix;

    _aspect = c._aspect;
    _viewport[0] = c._viewport[0];
    _viewport[1] = c._viewport[1];
    _viewport[2] = c._viewport[2];
    _viewport[3] = c._viewport[3];
    return *this;
}

void Camera::reset() {
    _trans = Vec3T<float>(0,0,0);
    _fovy = 34.f;
    //    _fovy = 5.f;
    
    _znear = .01f;
    _zfar = 100.f;
    _winscale = 1.0f;
    _wincenterx = 0.0f;
    _wincentery = 0.0f;
    _aspect = 1.0;
    identityModel();

    computeModelview();
    computeProjection();
}

void Camera::setFovy(float fovy_degrees) {
    _fovy = fovy_degrees;
    computeProjection();
}

void Camera::setPerspectiveParams(float fovy_degrees, float znear,float
                                  zfar) {
    _fovy = fovy_degrees;  
    _znear = znear;
    _zfar = zfar;
    computeProjection();
}


// ******************** OpenGL matrices ********************

// TODO: these need not use OpenGL
void Camera::computeModelview() {
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(_trans.x(), _trans.y(), _trans.z());

    glMultMatrixd( (double*)(_ballModel));
    glGetDoublev(GL_MODELVIEW_MATRIX,(double*)(_modelviewMatrix));
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

}

// TODO: these need not use OpenGL
void Camera::computeProjection() {

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    float ysize = tan(0.5*_fovy*M_PI/180.0)*_znear;
    float xsize = _aspect*ysize;
    float left   = _wincenterx -  xsize*_winscale; 
    float right  = _wincenterx +  xsize*_winscale; 
    float bottom = _wincentery -  ysize*_winscale; 
    float top    = _wincentery +  ysize*_winscale; 

    glFrustum(left, right, bottom, top, _znear, _zfar);
    //    glOrtho(left, right, bottom, top, _znear, _zfar);

    glGetDoublev(GL_PROJECTION_MATRIX,(double*)(_projectionMatrix));
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
  
    glMatrixMode(GL_MODELVIEW);
  
    //ofstream os("cam.log");
    //os<<_projectionMatrix<<endl;
    //os.close();
}

void Camera::loadMatrices()const {
    glMatrixMode(GL_PROJECTION);
    glLoadMatrixd((const double*)(_projectionMatrix));
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixd((const double*)(_modelviewMatrix));
}
    

// ******************** project/unproject ********************

// project a point from world space to screen space
Vec3T<float> Camera::project(const Vec3T<float>&  worldPoint) const {
    GLdouble winx, winy, winz;

    gluProject(worldPoint.x(),worldPoint.y(),worldPoint.z(),
               (const double*)(_modelviewMatrix),
               (const double*)(_projectionMatrix),
               _viewport,
               &winx,&winy,&winz);
    return Vec3T<float>(winx,winy,winz);
}

// compute the world space coordinate for an image coordinate
Vec3T<float> Camera::unproject(const Vec3T<float>& imagePoint) const {
    GLdouble objx, objy, objz;

    gluUnProject(imagePoint.x(),imagePoint.y(),imagePoint.z(),
                 (const double*)(_modelviewMatrix),
                 (const double*)(_projectionMatrix)
                 ,_viewport,
                 &objx,&objy,&objz);

    return Vec3T<float>(objx,objy,objz);
}

// compute the object space coordinate for an image coordinate
Vec3T<float> Camera::unproject(const HMatrix& model, const Vec3T<float>& imagePoint) const {
    GLdouble objx, objy, objz;

    gluUnProject(imagePoint.x(),imagePoint.y(),imagePoint.z(),
                 (const double*)(_modelviewMatrix * model),
                 (const double*)(_projectionMatrix)
                 ,_viewport,
                 &objx,&objy,&objz);

    return Vec3T<float>(objx,objy,objz);
}


// ******************** camera modification ********************

void Camera::translate(const Vec3T<float> & t) { _trans += t;}

void Camera::translateWindow(float dx, float dy) { 
    _wincentery -= 2.0*dy*tan(0.5*_fovy*M_PI/180.0)*_znear*_winscale;
    _wincenterx -= 2.0*dx*tan(0.5*_fovy*M_PI/180.0)*_aspect*_znear*_winscale;
    computeProjection();
}

void Camera::scaleWindow( float sfactor) { 
    _winscale = _winscale*pow(2.0f, sfactor);
    computeProjection();
}


// ******************** input/output  ********************




  void Camera::read(istream& is) 
  {
  is>>_projectionMatrix;
  is>>_modelviewMatrix;
  is>>_trans;
  is>>_ballModel;
  is>>_znear;
  is>>_zfar;
  is>>_fovy;
  is>>_aspect;
  is>>_wincenterx;
  is>>_wincentery;
  is>>_winscale;
  is>>_viewport[0];
  is>>_viewport[1];
  is>>_viewport[2];
  is>>_viewport[3];
  }

  void Camera::write(ostream& os) const {
  os<<_projectionMatrix<<endl;
  os<<_modelviewMatrix<<endl;
  os<<_trans<<endl;
  os<<_ballModel<<endl;
  os<<_znear<<endl;
  os<<_zfar<<endl;
  os<<_fovy<<endl;
  os<<_aspect<<endl;
  os<<_wincenterx<<endl;
  os<<_wincentery<<endl;
  os<<_winscale<<endl;
  os<<_viewport[0]<<" ";
  os<<_viewport[1]<<" ";
  os<<_viewport[2]<<" ";
  os<<_viewport[3]<<endl;
  }
/*
  ostream& operator<<(ostream& os, const Camera& camera ) {
  os << "TRANSFORMATIONS:"    << endl;  
  os << "Modelview matrix: "  << camera._modelviewMatrix << endl;
  os << "Projection matrix: " << camera._projectionMatrix << endl;

  os << "CAMERA PARAMETERS:"  << endl;
  os << "ballModel : "       << camera._ballModel << endl;
  os << "trans: "            << camera._trans  << endl;
  os << "aspect: "           << camera._aspect << endl;
  os << "znear: "            << camera._znear  <<  endl; 
  os << "zfar: "             << camera._zfar   << endl;
  os << "fovy: "             << camera._fovy   << endl;
  return os;
  }

  istream& operator>>(istream& is, Camera& ) { 
  die();
  //TODO  
  return  is;
  }
*/

// ********************** current matrix project *********

/*
  void Camera::getMatricesC()const {
  glGetDoublev(GL_MODELVIEW_MATRIX,(double*)_view);
  glGetDoublev(GL_PROJECTION_MATRIX,(double*)_proj);
  glGetIntegerv(GL_VIEWPORT,_vp);
  }

  Vec3T<float> Camera::projectC(const Vec3T<float>&  worldPoint) const {
  GLdouble winx, winy, winz;
  gluProject(worldPoint.x(),worldPoint.y(),worldPoint.z(),
  _view,_proj,_vp,&winx,&winy,&winz);
  return Vec3T<float>(winx,winy,winz);
  }

  Vec3T<float> Camera::unprojectC(const Vec3T<float>& screenPoint) const {
  GLdouble wx, wy, wz;
  gluUnProject(screenPoint.x(),screenPoint.y(),screenPoint.z(),
  _view, _proj, _vp, &wx, &wy, &wz);
  return Vec3T<float>(wx,wy,wz);
  }
*/

