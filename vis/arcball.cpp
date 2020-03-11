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

#include <GL/gl.h>
#include "psOpenGL.h"
#include <GL/glu.h>

#include "arcball.hpp"

#define LG_NSEGS 6
#define NSEGS (1<<LG_NSEGS)

#define RIMCOLOR()    glColor3d(1.0, 1.0, 1.0)
#define FARCOLOR()    glColor3d(0.75, 0.5, 0.125)
#define NEARCOLOR()   glColor3d(1.0, 1.0, 0.25)
#define DRAGCOLOR()   glColor3d(0.5, 1.0, 1.0)
#define RESCOLOR()    glColor3d(0.75, 0.125, 0.125)


// ******************** auxilary math. functions ****************************

// Convert window coordinates to sphere coordinates. 


Vec3T<float> ArcBall::RayOnSphere(const Vec3T<float>& ray_start, const Vec3T<float>& ray_dir) {
    Vec3T<float> ballpoint;
    float ray_proj = (-center+ray_start) * (ray_dir);
  
    Vec3T<float> ray_perp = (-center+ray_start) - ray_proj*ray_dir;
    if( ray_perp.l2() > radius ) {
        ballpoint = ray_perp.dir();
    } else {
        float start_center_distsq = (center-ray_start) * (center-ray_start);
        float Discr = max(0.0f, float(ray_proj*ray_proj - 
                                      start_center_distsq + radius*radius) );
        float t1 = -ray_proj + sqrt(Discr);
        float t2 = -ray_proj - sqrt(Discr);
        float t = (fabs(t1) > fabs(t2))? t2:t1;
        ballpoint = ray_start + t*ray_dir - center;
    }
    return (ballpoint);
}


/* Force sphere point to be perpendicular to axis. */
static Vec3T<float> ConstrainToAxis(const Vec3T<float>& loose, const Vec3T<float>& axis)
{
    Vec3T<float> onPlane;
    float norm;
    onPlane = loose - axis*(axis*loose);
    norm = onPlane.l2();
    if (norm > 0.0) {
        return onPlane.dir();
    } /* else drop through */
    if (axis.z() == 1.0f) {
        onPlane = Vec3T<float>(1.0f, 0.0f, 0.0f);
    } else {
        onPlane = (Vec3T<float>(-axis.y(), axis.x(), 0.0f)).dir();
    }
    return (onPlane);
}

/* Find the index of nearest arc of axis set. */
static int NearestConstraintAxis(Vec3T<float> loose, Vec3T<float> *axes, int nAxes)
{
    Vec3T<float> onPlane;
    float max, dp;
    int i, nearest;
    max = -1.0f; nearest = 0;
    for (i=0; i<nAxes; i++) {
        onPlane = ConstrainToAxis(loose, axes[i]);
        dp = onPlane * loose;
	 
        if (dp > max) {
            max = dp; nearest = i;
        }
    }
    return (nearest);
}

// ******************** drawing helpers ********************

// Draw an arc defined by its ends.
static void DrawAnyArc(Vec3T<float> vFrom, Vec3T<float> vTo)
{
    int i;
    Vec3T<float> pts[NSEGS+1];
    float dp;
    pts[0] = vFrom;
    pts[1] = pts[NSEGS] = vTo;
    for (i=0; i<LG_NSEGS; i++) 
        pts[1] = (pts[0] + pts[1]).dir();

    dp = 2.0f* pts[0] * pts[1];

    for (i=2; i<NSEGS; i++) {
        pts[i] = pts[i-1]*dp - pts[i-2];
    }

    glBegin(GL_LINE_STRIP);
    for (i=0; i<=NSEGS; i++)
        glVertex3fv(pts[i]);
    glEnd();
}

// Draw a big circle in the plane perp. to the vector
static void DrawFullArc(Vec3T<float> n)
{
    Vec3T<float> p(0,0,0), m(0,0,0);
    p.z() = 0.0f;
    if ((n.x() == 0.0f) && (n.y() == 0.0f)) {
        p.x() = 0.0f; p.y() = 1.0f;
    } else {
        p.x() = n.y(); 
        p.y() = -n.x();
    }
    p = p.dir(); 
  
    m = cross(p,n);
  
    glBegin(GL_LINES);
    glVertex3f(0.0,0.0,0.0);
    glVertex3fv(float(2.0)*n );
    glEnd();
  
    DrawAnyArc(p, m);
    DrawAnyArc(m, -p);
    DrawAnyArc(-p, -m);
    DrawAnyArc(-m,p);
}



// ******************** ArcBall memeber functions ********************

void ArcBall::Init()
{
    static Vec3T<float> unitVecs[3] = { Vec3T<float>(1.0f, 0.0f, 0.0f), 
                                        Vec3T<float>(0.0f, 1.0f, 0.0f), 
                                        Vec3T<float>(0.0f, 0.0f, 1.0f) };
    radius = 1.0f;
    center = Vec3T<float>(0.0,0.0,0.0);
  
  
    rayStartDown = rayStartNow = Vec3T<float>(0.0,0.0,0.0);
    rayDirDown = rayDirNow = Vec3T<float>(0.0,0.0,1.0);

  
    qDown = qNow = Quat(0.0,0.0,0.0,1.0);
    mNow.setIdentity();
    mDown.setIdentity();
   
    combined.setIdentity();

    showResult = dragging = GL_FALSE;
    axisSet = NoAxes;
    sets[WorldAxes][X] = unitVecs[X];     
    sets[WorldAxes][Y] = unitVecs[Y];     
    sets[WorldAxes][Z] = unitVecs[Z];     
    setSizes[WorldAxes] = 3;
    sets[BodyAxes][X] = mDown.col(X); 
    sets[BodyAxes][Y] = mDown.col(Y); 
    sets[BodyAxes][Z] = mDown.col(Z); 
    setSizes[BodyAxes] = 3;

    outlineColor  = Vec3T<float>(1.0f, 0.0f, 1.0f);
    inactiveColor = Vec3T<float>(0.75f, 0.5f, 0.125f);
    activeColor   = Vec3T<float>(1.0f, 1.0f, 0.25f);
    dragColor     = Vec3T<float>(0.5f, 1.0f, 1.0f);
    resultColor   = Vec3T<float>(0.75f, 0.125f, 0.125f);

    Update();
}


void ArcBall::Update()
{
    int setSize = setSizes[axisSet];
    Vec3T<float>* axset = sets[axisSet];

    vFrom = RayOnSphere(rayStartDown,rayDirDown).dir();
    vTo   = RayOnSphere(rayStartNow,rayDirNow).dir();

    if (dragging) {
        if (axisSet != NoAxes) {
            vFrom = ConstrainToAxis(vFrom, axset[axisIndex]);
            vTo   = ConstrainToAxis(  vTo, axset[axisIndex]);
        }
        qDrag = Quat(vFrom, vTo);
        qNow = qDrag* qDown;
    } else {
        if (axisSet != NoAxes) {
            axisIndex = NearestConstraintAxis(vTo, axset, setSize);
        }
    }
    qDown.ballPoints(vrFrom, vrTo);
    SetRotation( HMatrix( qNow ) );


    if( X != axisIndex ) sets[BodyAxes][X] = mNow.col(X); 
    if( Y != axisIndex ) sets[BodyAxes][Y] = mNow.col(Y); 
    if( Z != axisIndex ) sets[BodyAxes][Z] = mNow.col(Z); 

}


void ArcBall::BeginDrag() { 
    dragging = GL_TRUE; 
    rayStartDown = rayStartNow; 
    rayDirDown = rayDirNow;
}

void ArcBall::EndDrag() {  
    dragging = GL_FALSE;
    qDown = qNow; mDown = mNow;
    sets[BodyAxes][X] = mDown.col(X); 
    sets[BodyAxes][Y] = mDown.col(Y); 
    sets[BodyAxes][Z] = mDown.col(Z); 
}

void ArcBall::DrawCompleteBall(const Vec3T<float>& camerapos) {
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glTranslatef(center.x(),center.y(), center.z());  

    DrawOuterRing(camerapos);
    glScalef(radius, radius, radius);
  
    DrawResultArc();
    DrawConstraints(); 
    DrawDragArc();
  
    glPopMatrix();

}

// Draw the controller with all its arcs.
void ArcBall::Draw(const Vec3T<float>& camerapos) {
  
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glTranslatef(center.x(),center.y(), center.z());

    if (axisSet != NoAxes)  DrawOuterRing(camerapos);

    glScalef(radius, radius, radius);
  
    DrawResultArc();
    DrawConstraints();
    DrawDragArc();

    glPopMatrix();

}


// Draw all constraint arcs.
void ArcBall::DrawConstraints() {

    if (axisSet==NoAxes) return;

    Vec3T<float>* axset = sets[axisSet];

    int axisI, setSize = setSizes[axisSet];

    for (axisI=0; axisI < setSize; axisI++) {
        if (axisIndex != axisI) {
            //        if (dragging) continue;
            glColor3fv(inactiveColor);
        } else  glColor3fv(activeColor);
        DrawFullArc(axset[axisI]);
    }
}

void ArcBall::DrawOuterRing(const Vec3T<float>& camerapos) {

    glColor3fv(outlineColor);
    Vec3T<float> cameradir = (camerapos - center).dir();
    float d = (camerapos - center).l2();


    glPushMatrix();
    Vec3T<float> displ = cameradir*float((radius*radius/d));
    glTranslatef(displ.x(), displ.y(), displ.z() );
    float newradius = radius*sqrt(1.0f - radius*radius/(d*d) );

    glScalef( newradius, newradius, newradius );
    DrawFullArc(cameradir);
    glPopMatrix();

}

// Draw "rubber band" arc during dragging.
void ArcBall::DrawDragArc() {
    if (dragging) { 
        glColor3fv(dragColor);
        DrawAnyArc(vFrom, vTo);
    }
}

// Draw arc for result of all drags.
void ArcBall::DrawResultArc()
{
    glColor3fv(resultColor);
    if (showResult) DrawAnyArc(vrFrom, vrTo);
}



