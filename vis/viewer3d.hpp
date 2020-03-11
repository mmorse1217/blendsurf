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

#ifndef _VIEWER3D_HPP_
#define _VIEWER3D_HPP_

#include "viewer.hpp"
#include "camera.hpp"

class GeoObject;

typedef Vec3T<float> vec3f;

class Viewer3D : public Viewer 
{
public:  
    Viewer3D(char* name = 0, int w = 512, int h = 512);
    virtual ~Viewer3D();
    // viewed object
    void setObject(GeoObject*);  //virtual void addObject(GeoObject* geoObject);  //int numObj() const { return _objvec.size(); }
    GeoObject* getObject() const { return _obj; }  //const vector<GeoObject*>& getObjVec() const { return _objvec; }
    vec3f minPoint() const;
    vec3f maxPoint() const;
    vec3f centerPoint() const;
    // camera
    Camera* getCamera() { return _camera; }
  void setCamera(Camera* cam); 

protected:
    virtual void display();
    virtual void idle();
    virtual void reshape(int w, int h);
public:
    // user interface components
    Camera* _camera;
    GeoObject* _obj;
protected:
    bool _needPosition;
    void positionCamera();
};

#endif
