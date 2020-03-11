/// \file
#ifndef _QUATERNION_HPP_
#define _QUATERNION_HPP_

#include "nummatrix.hpp"
#include "cvect.hpp"

// FREE SOURCE CODE
// http://www.magic-software.com/License/free.pdf


 
/// Implments quaritenions
class Quaternion {
public:
    double w, x, y, z;
    // construction and destruction
    Quaternion (double* f);
    Quaternion (double fW = 1.0, double fX = 0.0, double fY = 0.0,  double fZ = 0.0);
    Quaternion (const Quaternion& rkQ);
    Quaternion (const CPoint &p);
    void normalize();
  
    // conversion between quaternions, matrices, and angle-axes
    void fromRotationMatrix (const NumMatrix& kRot);
    void toRotationMatrix (NumMatrix& kRot) const;
    void fromAngleAxis (const double& rfAngle, const CPoint& rkAxis);
    void toAngleAxis (double& rfAngle, CPoint& rkAxis) const;
    void fromAxes (const CPoint* akAxis);
    void toAxes (CPoint* akAxis) const;
  
    /// arithmetic operations
    Quaternion& operator= (const Quaternion& rkQ);
    Quaternion operator+ (const Quaternion& rkQ) const;
    Quaternion operator- (const Quaternion& rkQ) const;
    Quaternion &operator+= (const Quaternion& rkQ);
    Quaternion &operator-= (const Quaternion& rkQ);
    Quaternion operator* (const Quaternion& rkQ) const;
    Quaternion operator* (double fScalar) const;
    friend Quaternion operator* (double fScalar, const Quaternion& rkQ);
    friend Quaternion operator* (const CPoint &rkVector, const Quaternion& rkQ);  
    Quaternion operator- () const;
  
    // functions of a quaternion
  
    double dot (const Quaternion& rkQ) const;  ///< dot product
    double norm () const;  ///< squared-length
    Quaternion inverse () const;  
    Quaternion unitInverse () const;  ///< apply to unit-length quaternion
    Quaternion exp () const;
    Quaternion log () const;

    // rotation of a vector by a quaternion
    CPoint operator* (const CPoint& rkVector) const;
  
    // spherical linear interpolation
    static Quaternion slerp (double fT, const Quaternion& rkP,
                             const Quaternion& rkQ);

    static Quaternion slerpExtraSpins (double fT,
                                       const Quaternion& rkP, const Quaternion& rkQ,
                                       int iExtraSpins);

    // setup for spherical quadratic interpolation
    static void intermediate (const Quaternion& rkQ0,
                              const Quaternion& rkQ1, const Quaternion& rkQ2,
                              Quaternion& rka, Quaternion& rkB);

    // spherical quadratic interpolation
    static Quaternion squad (double fT, const Quaternion& rkP,
                             const Quaternion& rkA, const Quaternion& rkB,
                             const Quaternion& rkQ);

    // special values
    static const Quaternion ZERO;
    static const Quaternion IDENTITY;
};

#endif



