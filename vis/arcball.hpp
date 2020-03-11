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

#ifndef _ARCBALL_HPP_
#define _ARCBALL_HPP_

#include "quat.hpp"
#include "hmatrix.hpp"

#define MAX_SET_SIZE 3 

//: Arcball rotation for UI
class ArcBall 
{ 
public:  
    enum AxisSet {NoAxes, WorldAxes, BodyAxes, NSets};
    enum AxesNames { X = 0,Y = 1 ,Z = 2};
  
    void Update();
    void Init();

    ArcBall() { Init(); }


    void SetRadius( double r ) { radius = r; } 
    double Radius() const { return radius; }
    void SetPosition( const Vec3T<float>& c) { center = c; 
    combined = HMatrix::Translation(center)*mNow;}
    const Vec3T<float>& Position() const { return center;  }

    void SetRotation( const HMatrix& rot ) { mNow = rot; 
    combined = HMatrix::Translation(center)*mNow;} 
    const HMatrix&  Rotation() { return mNow; }

    void PickingRay( const Vec3T<float>& rs, const Vec3T<float>& rd)
        { rayStartNow = rs; rayDirNow = rd; }
  
    void UseSet(AxisSet aset) { 
        if (!dragging) axisSet = aset; 
    }
  
    // appearance
    void SetInactiveColor(const Vec3T<float>& color) { inactiveColor = color;  }

    void SetActiveColor  (const Vec3T<float>& color) { activeColor = color;  }

    void SetDragColor    (const Vec3T<float>& color) { dragColor = color;  }

    void SetOutlineColor (const Vec3T<float>& color) { outlineColor = color;  }

    void SetResultColor  (const Vec3T<float>& color) { resultColor = color;  }


    // Begin drawing arc for all drags combined
    void ShowResult() { showResult = GL_TRUE; }
    // Stop drawing arc for all drags combined. 
    void HideResult() { showResult = GL_FALSE; }

    void BeginDrag();
    void EndDrag();
    void Draw(const Vec3T<float>& camerapos );
    void DrawCompleteBall(const Vec3T<float>& camerapos );

private:
    Vec3T<float> RayOnSphere(const Vec3T<float>& ray_start, 
                             const Vec3T<float>& ray_dir);

    void DrawOuterRing(const Vec3T<float>& camerapos);
    void DrawConstraints();
    void DrawDragArc();
    void DrawResultArc();

    // data:
    double radius;
  
    // position  the world coord system
    Vec3T<float> center;
  
    //  combined transformation matrix
    HMatrix combined;

    Quat  qNow, qDown, qDrag;
    // current picking ray and the picking ray when drag was started
    Vec3T<float> rayStartNow, rayDirNow,rayStartDown, rayDirDown;
    Vec3T<float> vFrom, vTo, vrFrom, vrTo;
    HMatrix mNow, mDown;
    bool showResult, dragging;
    Vec3T<float> sets[NSets][MAX_SET_SIZE];
    int setSizes[NSets];
    AxisSet axisSet;
    int axisIndex;

    Vec3T<float> outlineColor;
    Vec3T<float> inactiveColor;
    Vec3T<float> activeColor;
    Vec3T<float> dragColor;
    Vec3T<float> resultColor;

};

#endif /* _ARCBALL_H_ */







