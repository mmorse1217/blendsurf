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

#ifndef _UIACTION_HPP_
#define _UIACTION_HPP_

//#include "compat.h"
#include "arcball.hpp"
#include "camera.hpp"


class UiAction 
{
public:
    typedef Vec3T<float> vec3f;
public:
    virtual void update(int , int ) { ; }
    virtual ~UiAction() { ; }
};

class CameraRotateAction : public UiAction {
public:
    CameraRotateAction(int x, int y, ArcBall* arcball, Camera* camera) :
            _x(x), _y(y), _arcball(arcball), _camera(camera) { 

        vec3f worldMousePos = _camera->unproject(vec3f(x,y,0));
        _arcball->PickingRay(_camera->viewWorldPosition(), 
                             (worldMousePos - _camera->viewWorldPosition()).dir());
        _arcball->BeginDrag();
        _arcball->Update();
    }

    virtual ~CameraRotateAction() {
        vec3f cp = _arcball->Position();
        _camera->model() = _camera->model() * HMatrix::Translation(cp)
            * _arcball->Rotation() * HMatrix::Translation(-cp);

        _camera->computeModelview();
        _camera->loadMatrices();
    
        _arcball->SetRotation(HMatrix::Identity());
    }
  
    virtual void update(int x, int y) {

        vec3f worldMousePos = _camera->unproject(vec3f(x,y,0));
        _arcball->PickingRay(_camera->viewWorldPosition(), 
                             (worldMousePos - _camera->viewWorldPosition()).dir());
        _arcball->Update();
    }

  
private:
    int _x, _y;
    ArcBall* _arcball;
    Camera* _camera;
};


class CameraTransXYAction : public UiAction {
public:
  
    CameraTransXYAction(int x, int y, ArcBall* arcball, Camera* camera) :
            _x(x), _y(y), _arcball(arcball), _camera(camera) { 
    }

    virtual ~CameraTransXYAction() {
    }

    virtual void update(int x, int y) {
        float z = _camera->project(_arcball->Position()).z();
        vec3f p0 = _camera->unproject(vec3f(_x, _y, z));
        vec3f p1 = _camera->unproject(vec3f(x, y, z));
        _camera->model() = _camera->model() * HMatrix::Translation(p1-p0);
    
        _camera->computeModelview();
        _camera->loadMatrices();

        _x = x;
        _y = y;
    }

private:
    int _x, _y;
    ArcBall* _arcball;
    Camera* _camera;
};


class CameraTransZAction : public UiAction {
public:
  
    CameraTransZAction(int x, int y, ArcBall* arcball, Camera* camera) :
            _x(x), _y(y), _arcball(arcball), _camera(camera) { 
    }

    virtual ~CameraTransZAction() {
    }

    virtual void update(int x, int y) {
        vec3f p0 = _camera->unproject(vec3f(x, _y, _camera->project(_arcball->Position()).z()));
        vec3f p1 = _camera->unproject(vec3f(x,  y, _camera->project(_arcball->Position()).z()));
        if (y<_y)
            _camera->model() = _camera->model() * 
                HMatrix::Translation((p0-p1).l2() * (_camera->viewWorldPosition()-_arcball->Position()).dir());
        else
            _camera->model() = _camera->model() * 
                HMatrix::Translation(-(p0-p1).l2() * (_camera->viewWorldPosition()-_arcball->Position()).dir());
        _x = x;
        _y = y;
        _camera->computeModelview();
        _camera->loadMatrices();
    }
  
private:
    int _x, _y;
    ArcBall* _arcball;
    Camera* _camera;
};

#endif
