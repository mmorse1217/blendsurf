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

#ifndef _HMATRIX_HPP_
#define _HMATRIX_HPP_

#include "quat.hpp"


class HMatrix;
//inline ostream& operator<< ( ostream& os, const HMatrix& M );
//inline istream& operator>> ( istream& os,  HMatrix& M );

class HMatrix 
{
private:
    enum AxesNames { X = 0,Y = 1 ,Z = 2, W = 3};
public:
    HMatrix&  setIdentity() {     
        memset(_m,0, 16*sizeof(double));
        _m[0] = _m[5] = _m[10] = _m[15] = 1.0f;
        return *this;
    }
  
    static HMatrix Identity() {
        return HMatrix().setIdentity();
    }
  
    HMatrix() { 
        setIdentity();
    }
  
    HMatrix(const HMatrix& mat) { 
        memcpy(_m,mat._m,16*sizeof(double));
    }
  
    HMatrix( float m00,float m01,float m02,float m03,
             float m10,float m11,float m12,float m13,
             float m20,float m21,float m22,float m23,
             float m30,float m31,float m32,float m33) {
        _m[0]  = m00, _m[1]  = m10,_m[2]   = m20,_m[3]  = m30;
        _m[4]  = m01, _m[5]  = m11,_m[6]   = m21,_m[7]  = m31;
        _m[8]  = m02, _m[9]  = m12,_m[10]  = m22,_m[11] = m32;
        _m[12] = m03, _m[13] = m13,_m[14]  = m23,_m[15] = m33;
    }
  
    operator double*() { 
        return _m;
    }
  
    operator const double*() const { 
        return _m;
    }
  
    explicit HMatrix(double* m ) { 
        memcpy(_m, m,16*sizeof(double));
    }
  
    HMatrix& setTranslation( const Vec3T<float>& trans ) { 
        setIdentity(); 
        (*this)(0,3) = trans.x();
        (*this)(1,3) = trans.y();
        (*this)(2,3) = trans.z();
        return *this;
    }
  
    static HMatrix Translation( const Vec3T<float>& trans ) {
        return HMatrix().setTranslation(trans);
    }
  
    double operator()(int i, int j) const { 
        return _m[4*j + i ];
    }
  
    double& operator()(int i, int j)  { 
        return _m[4*j + i ];
    }
  
    HMatrix& operator=(const HMatrix& mat) { 
        memcpy(_m, mat._m, 16*sizeof(double));
        return *this;
    }
    Vec3T<float> col(int i) { 
        return Vec3T<float>(_m[4*i], _m[4*i+1],  _m[4*i+2]);
    }

    Vec3T<float> row(int i) { 
        return Vec3T<float>(_m[i], _m[i+4],  _m[i+8]);
    }
  
    // create a rotation matrix for 

    static HMatrix Rotation(const Vec3T<float>& v, float a) {
        //Vec3T<float> u = v.dir();
        Vec3T<float> u = v/v.l2();
        float u1 = u.x();
        float u2 = u.y();
        float u3 = u.z();
	 
        HMatrix U = HMatrix(u1*u1, u1*u2, u1*u3, 0,
                            u2*u1, u2*u2, u2*u3, 0,
                            u3*u1, u3*u2, u3*u3, 0,
                            0,     0,     0,     1);
        HMatrix S = HMatrix( 0, -u3, u2, 0,
                             u3,  0, -u1, 0,
                             -u2,  u1, 0,  0,
                             0,    0, 0,  1);

        HMatrix tmp = U + (Identity() - U) * cos(a)  + S * sin(a);
        tmp(3,3) = 1.0;
        return tmp;
    }
  
    // create a rotation matrix for 
    /*****************************************************************************
     * DESCRIPTION:                                                               M
     * Routine to compute the INVERSE of a given matrix M which is not modified.  M
     *   The matrix is assumed to be 4 by 4 (transformation matrix).		     M
     *   Return TRUE if inverted matrix (InvM) do exists.			     M
     *                                                                            *
     * PARAMETERS:                                                                M
     *   M:          Original matrix to invert.                                   M
     *   InvM:       Inverted matrix will be placed here.                         M
     *                                                                            *
     * RETURN VALUE:                                                              M
     *   int:        TRUE if inverse exists, FALSE otherwise.                     M
     *                                                                            *
     * KEYWORDS:                                                                  M
     *   MatInverseMatrix, transformations, matrix inverse                        M
     *****************************************************************************/
  
    bool setInv(const HMatrix& M) {
        HMatrix A;
        int i, j, k;
        double V;
    
        A = M;
        setIdentity();
	 
    
        for (i = 0; i < 4; i++) {
            V = A(i,i);				      /* Find the new pivot. */
            k = i;
            for (j = i + 1; j < 4; j++) 
                if (abs(A(j,i)) > abs(V)) {
                    /* Find maximum on col i, row i+1..n */
                    V = A(j,i);
                    k = j;
                }
            j = k;
		

            double tmp;
            if (i != j)
                for (k = 0; k < 4; k++) {
                    tmp = A(i,k); A(i,k) = A(j,k); A(j,k) = tmp;
                    tmp = (*this)(i,k); (*this)(i,k) = (*this)(j,k); (*this)(j,k) = tmp;
                }
		
		
            for (j = i + 1; j < 4; j++) {	 /* Eliminate col i from row i+1..n. */
                if(A(j,i) != 0) {
                    V = A(j,i) / A(i,i);
			 
                    for (k = 0; k < 4; k++) {
                        A(j,k)    -= V * A(i,k);
                        (*this)(j,k) -= V * (*this)(i,k);
                    }
                }
		  
            }
        }
	 
        for (i = 3; i >= 0; i--) {			       /* Back Substitution. */
            if (A(i,i) == 0)
                return false;					   /* Error. */
		
            for (j = 0; j < i; j++) {	 /* Eliminate col i from row 1..i-1. */
                V = A(j,i) / A(i,i);
		  
                for (k = 0; k < 4; k++) {
                    /* A[j][k] -= V * A[i][k]; */
                    (*this)(j,k) -= V * (*this)(i,k);
                }
            }
        }
    
        for (i = 0; i < 4; i++)		    /* Normalize the inverse Matrix. */
            for (j = 0; j < 4; j++)
                (*this)(i,j) /= A(i,i);
    
        return true;
    }  
  
    // Construct rotation matrix from (possibly non-unit) quaternion.
    // Works correctly for right-handed coordinate system
    // and right-handed rotations.
  
  
    explicit HMatrix( const Quat q )  {
        double Nq = q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w;
        double s = (Nq > 0.0) ? (2.0 / Nq) : 0.0;
        double xs = q.x*s,	      ys = q.y*s,	  zs = q.z*s;
        double wx = q.w*xs,	      wy = q.w*ys,	  wz = q.w*zs;
        double xx = q.x*xs,	      xy = q.x*ys,	  xz = q.x*zs;
        double yy = q.y*ys,	      yz = q.y*zs,	  zz = q.z*zs;
    
        (*this)(X,X) = 1.0 - (yy + zz); 
        (*this)(Y,X) = xy + wz; 
        (*this)(Z,X) = xz - wy;
        (*this)(X,Y) = xy - wz; 
        (*this)(Y,Y) = 1.0 - (xx + zz); 
        (*this)(Z,Y) = yz + wx;
        (*this)(X,Z) = xz + wy; 
        (*this)(Y,Z) = yz - wx; 
        (*this)(Z,Z) = 1.0 - (xx + yy);
        (*this)(X,W) = (*this)(Y,W) = (*this)(Z,W) = (*this)(W,X) =  (*this)(W,Y) = (*this)(W,Z) = 0.0;
        (*this)(W,W) = 1.0;
    }
  
    Vec3T<float> operator* (const Vec3T<float>& v ) const { 
        const HMatrix& M = *this;
        return Vec3T<float>( 
            v.x() * M(0,0) + v.y() * M(0,1) + v.z() * M(0,2) + M(0,3),
            v.x() * M(1,0) + v.y() * M(1,1) + v.z() * M(1,2) + M(1,3),
            v.x() * M(2,0) + v.y() * M(2,1) + v.z() * M(2,2) + M(2,3));
    }
  
    HMatrix operator* (float s) const {
        const HMatrix& M = *this;
        return HMatrix(s*M(0,0), s*M(0,1), s*M(0,2), s*M(0,3),
                       s*M(1,0), s*M(1,1), s*M(1,2), s*M(1,3),
                       s*M(2,0), s*M(2,1), s*M(2,2), s*M(2,3),
                       s*M(3,0), s*M(3,1), s*M(3,2), s*M(3,3));
    }
  
  
  
    HMatrix operator* (const HMatrix& M2 ) const {     
        const HMatrix& M1 = *this;
        return HMatrix(
            M1(0,0)*M2(0,0) + M1(0,1)*M2(1,0) + M1(0,2)*M2(2,0) + M1(0,3)*M2(3,0),
            M1(0,0)*M2(0,1) + M1(0,1)*M2(1,1) + M1(0,2)*M2(2,1) + M1(0,3)*M2(3,1),
            M1(0,0)*M2(0,2) + M1(0,1)*M2(1,2) + M1(0,2)*M2(2,2) + M1(0,3)*M2(3,2),
            M1(0,0)*M2(0,3) + M1(0,1)*M2(1,3) + M1(0,2)*M2(2,3) + M1(0,3)*M2(3,3),
						 
            M1(1,0)*M2(0,0) + M1(1,1)*M2(1,0) + M1(1,2)*M2(2,0) + M1(1,3)*M2(3,0),
            M1(1,0)*M2(0,1) + M1(1,1)*M2(1,1) + M1(1,2)*M2(2,1) + M1(1,3)*M2(3,1),
            M1(1,0)*M2(0,2) + M1(1,1)*M2(1,2) + M1(1,2)*M2(2,2) + M1(1,3)*M2(3,2),
            M1(1,0)*M2(0,3) + M1(1,1)*M2(1,3) + M1(1,2)*M2(2,3) + M1(1,3)*M2(3,3),
						 
            M1(2,0)*M2(0,0) + M1(2,1)*M2(1,0) + M1(2,2)*M2(2,0) + M1(2,3)*M2(3,0),
            M1(2,0)*M2(0,1) + M1(2,1)*M2(1,1) + M1(2,2)*M2(2,1) + M1(2,3)*M2(3,1),
            M1(2,0)*M2(0,2) + M1(2,1)*M2(1,2) + M1(2,2)*M2(2,2) + M1(2,3)*M2(3,2),
            M1(2,0)*M2(0,3) + M1(2,1)*M2(1,3) + M1(2,2)*M2(2,3) + M1(2,3)*M2(3,3),
						 
            M1(3,0)*M2(0,0) + M1(3,1)*M2(1,0) + M1(3,2)*M2(2,0) + M1(3,3)*M2(3,0),
            M1(3,0)*M2(0,1) + M1(3,1)*M2(1,1) + M1(3,2)*M2(2,1) + M1(3,3)*M2(3,1),
            M1(3,0)*M2(0,2) + M1(3,1)*M2(1,2) + M1(3,2)*M2(2,2) + M1(3,3)*M2(3,2),
            M1(3,0)*M2(0,3) + M1(3,1)*M2(1,3) + M1(3,2)*M2(2,3) + M1(3,3)*M2(3,3)
            );
    }
  
    HMatrix operator+ (const HMatrix& M2) const {
        const HMatrix& M1 = *this;
        return HMatrix(
            M1(0,0)+M2(0,0), M1(0,1)+M2(0,1), M1(0,2)+M2(0,2), M1(0,3)+M2(0,3),
            M1(1,0)+M2(1,0), M1(1,1)+M2(1,1), M1(1,2)+M2(1,2), M1(1,3)+M2(1,3),
            M1(2,0)+M2(2,0), M1(2,1)+M2(2,1), M1(2,2)+M2(2,2), M1(2,3)+M2(2,3),
            M1(3,0)+M2(3,0), M1(3,1)+M2(3,1), M1(3,2)+M2(3,2), M1(3,3)+M2(3,3));
    }
  
    HMatrix operator- (const HMatrix& M2) const {
        const HMatrix& M1 = *this;
        return HMatrix(
            M1(0,0)-M2(0,0), M1(0,1)-M2(0,1), M1(0,2)-M2(0,2), M1(0,3)-M2(0,3),
            M1(1,0)-M2(1,0), M1(1,1)-M2(1,1), M1(1,2)-M2(1,2), M1(1,3)-M2(1,3),
            M1(2,0)-M2(2,0), M1(2,1)-M2(2,1), M1(2,2)-M2(2,2), M1(2,3)-M2(2,3),
            M1(3,0)-M2(3,0), M1(3,1)-M2(3,1), M1(3,2)-M2(3,2), M1(3,3)-M2(3,3));
    }
  
    friend inline ostream& operator<< ( ostream& os, const HMatrix& M );
    friend inline istream& operator>> ( istream& is, HMatrix& M );
  
private:
    double _m[16];
};

  inline istream& operator>> ( istream& is, HMatrix& M ) {
  int i;
  char c;
  is>>c; // read '['
  for(i = 0; i < 4; ++i)
  is>>M._m[i];
  is>>c; // read ';'
  for(i = 4; i < 8; ++i)
  is>>M._m[i];
  is>>c; // read ';'
  for(i = 8; i < 12; ++i)
  is>>M._m[i];
  is>>c; // read ';'
  for(i = 12; i < 16; ++i)
  is>>M._m[i];
  is>>c; // read ']'
  return is;
  }
  inline ostream& operator<< ( ostream& os, const HMatrix& M ) { 
  os << "[ " << M._m[0]  << " " << M._m[1]  << " " << M._m[2]  << " " << M._m[3]  << "; ";
  os         << M._m[4]  << " " << M._m[5]  << " " << M._m[6]  << " " << M._m[7]  << "; ";
  os         << M._m[8]  << " " << M._m[9]  << " " << M._m[10] << " " << M._m[11] << "; ";
  os         << M._m[12] << " " << M._m[13] << " " << M._m[14] << " " << M._m[15] 
  << "] ";   
  return os;
  }

#endif



