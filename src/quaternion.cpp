/// \file

// FREE SOURCE CODE
// http://www.magic-software.com/License/free.pdf

#include "quaternion.hpp"

static double gs_fEpsilon = 1e-03;
const Quaternion Quaternion::ZERO(0.0,0.0,0.0,0.0);
const Quaternion Quaternion::IDENTITY(1.0,0.0,0.0,0.0);

//----------------------------------------------------------------------------
Quaternion::Quaternion(double* f)
{
    w = f[0];  x = f[1];  y = f[2];  z = f[3];
}

//----------------------------------------------------------------------------
Quaternion::Quaternion (double fW, double fX, double fY, double fZ)
{
    w = fW;
    x = fX;
    y = fY;
    z = fZ;
}

//----------------------------------------------------------------------------
Quaternion::Quaternion (const CPoint &p)
{
    w = 0;
    x = p[0];
    y = p[1];
    z = p[2];
}

//----------------------------------------------------------------------------
Quaternion::Quaternion (const Quaternion& rkQ)
{
    w = rkQ.w;
    x = rkQ.x;
    y = rkQ.y;
    z = rkQ.z;
}

//---------------------------------------------------------------------------
void Quaternion::normalize()
{
    double nrm= w*w+x*x+y*y+z*z;
    if( nrm > 0.0){
        nrm = 1.0/nrm;
        w *= nrm;
        x *= nrm;
        y *= nrm;
        z *= nrm;
    }
}

//----------------------------------------------------------------------------
void Quaternion::fromRotationMatrix (const NumMatrix& kRot)
{
    // Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
    // article "Quaternion Calculus and Fast Animation".

    //    double fTrace = kRot[0][0]+kRot[1][1]+kRot[2][2];
    double fTrace = kRot(0,0)+kRot(1,1)+kRot(2,2);
    double fRoot;

    if ( fTrace > 0.0 )
	{
	    // |w| > 1/2, may as well choose w > 1/2
	    fRoot = sqrt(fTrace + 1.0);  // 2w
	    w = 0.5f*fRoot;
	    fRoot = 0.5f/fRoot;  // 1/(4w)
	    x = (kRot(2,1)-kRot(1,2))*fRoot;
	    y = (kRot(0,2)-kRot(2,0))*fRoot;
	    z = (kRot(1,0)-kRot(0,1))*fRoot;
	}
    else
	{
	    // |w| <= 1/2
	    static int s_iNext[3] = { 1, 2, 0 };
	    int i = 0;
	    if ( kRot(1,1) > kRot(0,0) )
            i = 1;
	    if ( kRot(2,2) > kRot(i,i) )
            i = 2;
	    int j = s_iNext[i];
	    int k = s_iNext[j];

	    fRoot = sqrt(kRot(i,i)-kRot(j,j)-kRot(k,k) + 1.0);
	    double* apkQuat[3] = { &x, &y, &z };
	    *apkQuat[i] = 0.5f*fRoot;
	    fRoot = 0.5f/fRoot;
	    w = (kRot(k,j)-kRot(j,k))*fRoot;
	    *apkQuat[j] = (kRot(j,i)+kRot(i,j))*fRoot;
	    *apkQuat[k] = (kRot(k,i)+kRot(i,k))*fRoot;
	}
}
//----------------------------------------------------------------------------
void Quaternion::toRotationMatrix (NumMatrix& kRot) const
{
    kRot.resize(3,3);
    double fTx  = 2.0*x;
    double fTy  = 2.0*y;
    double fTz  = 2.0*z;
    double fTwx = fTx*w;
    double fTwy = fTy*w;
    double fTwz = fTz*w;
    double fTxx = fTx*x;
    double fTxy = fTy*x;
    double fTxz = fTz*x;
    double fTyy = fTy*y;
    double fTyz = fTz*y;
    double fTzz = fTz*z;
  
    kRot(0,0) = 1.0-(fTyy+fTzz);
    kRot(0,1) = fTxy-fTwz;
    kRot(0,2) = fTxz+fTwy;
    kRot(1,0) = fTxy+fTwz;
    kRot(1,1) = 1.0-(fTxx+fTzz);
    kRot(1,2) = fTyz-fTwx;
    kRot(2,0) = fTxz-fTwy;
    kRot(2,1) = fTyz+fTwx;
    kRot(2,2) = 1.0-(fTxx+fTyy);
}
//----------------------------------------------------------------------------
void Quaternion::fromAngleAxis (const double& rfAngle, const CPoint& rkAxis)
{
    // assert:  axis[] is unit length
    //
    // The quaternion representing the rotation is
    //   q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)

    double fHalfAngle = 0.5f*rfAngle;
    double fSin = sin(fHalfAngle);
    w = cos(fHalfAngle);
    x = fSin*rkAxis[0];
    y = fSin*rkAxis[1];
    z = fSin*rkAxis[2];
}
//----------------------------------------------------------------------------
void Quaternion::toAngleAxis (double& rfAngle, CPoint& rkAxis) const
{
    // The quaternion representing the rotation is
    //   q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)

    double fSqrLength = x*x+y*y+z*z;
    if ( fSqrLength > 0.0 )
	{
	    rfAngle = 2.0*acos(w);
	    double fInvLength = 1.0/sqrt(fSqrLength);
	    rkAxis[0] = x*fInvLength;
	    rkAxis[1] = y*fInvLength;
	    rkAxis[2] = z*fInvLength;
	}
    else
	{
	    // angle is 0 (mod 2*pi), so any axis will do
	    rfAngle = 0.0;
	    rkAxis[0] = 1.0;
	    rkAxis[1] = 0.0;
	    rkAxis[2] = 0.0;
	}
}
//----------------------------------------------------------------------------
void Quaternion::fromAxes (const CPoint* akAxis)
{
    NumMatrix kRot;
    kRot.resize(3,3);

    for (int iCol = 0; iCol < 3; iCol++)
	{
	    kRot(0,iCol) = akAxis[iCol][0];
	    kRot(1,iCol) = akAxis[iCol][1];
	    kRot(2,iCol) = akAxis[iCol][2];
	}

    fromRotationMatrix(kRot);
}
//----------------------------------------------------------------------------
void Quaternion::toAxes (CPoint* akAxis) const
{
    NumMatrix kRot;
    kRot.resize(3,3);	 

    toRotationMatrix(kRot);

    for (int iCol = 0; iCol < 3; iCol++)
	{
	    akAxis[iCol][0] = kRot(0,iCol);
	    akAxis[iCol][1] = kRot(1,iCol);
	    akAxis[iCol][2] = kRot(2,iCol);
	}
}
//----------------------------------------------------------------------------
Quaternion& Quaternion::operator= (const Quaternion& rkQ)
{
    w = rkQ.w;
    x = rkQ.x;
    y = rkQ.y;
    z = rkQ.z;
    return *this;
}
//----------------------------------------------------------------------------
Quaternion Quaternion::operator+ (const Quaternion& rkQ) const
{
    return Quaternion(w+rkQ.w,x+rkQ.x,y+rkQ.y,z+rkQ.z);
}
//----------------------------------------------------------------------------
Quaternion Quaternion::operator- (const Quaternion& rkQ) const
{
    return Quaternion(w-rkQ.w,x-rkQ.x,y-rkQ.y,z-rkQ.z);
}

//----------------------------------------------------------------------------
Quaternion &Quaternion::operator+=(const Quaternion& rkQ) 
{
    w+=rkQ.w; x+=rkQ.x; y+=rkQ.y; z+=rkQ.z;
    return *this;
}

//----------------------------------------------------------------------------
Quaternion &Quaternion::operator-=(const Quaternion& rkQ) 
{
    w-=rkQ.w; x-=rkQ.x; y-=rkQ.y; z-=rkQ.z;
    return *this;
}


//----------------------------------------------------------------------------
Quaternion Quaternion::operator* (const Quaternion& rkQ) const
{
    // NOTE:  Multiplication is not generally commutative, so in most
    // cases p*q != q*p.
  
    return Quaternion
        (
            w*rkQ.w-x*rkQ.x-y*rkQ.y-z*rkQ.z,
            w*rkQ.x+x*rkQ.w+y*rkQ.z-z*rkQ.y,
            w*rkQ.y+y*rkQ.w+z*rkQ.x-x*rkQ.z,
            w*rkQ.z+z*rkQ.w+x*rkQ.y-y*rkQ.x
            );
}
//----------------------------------------------------------------------------
Quaternion Quaternion::operator* (double fScalar) const
{
    return Quaternion(fScalar*w,fScalar*x,fScalar*y,fScalar*z);
}
//----------------------------------------------------------------------------
Quaternion operator* (double fScalar, const Quaternion& rkQ)
{
    return Quaternion(fScalar*rkQ.w,fScalar*rkQ.x,fScalar*rkQ.y,
                      fScalar*rkQ.z);
}

//----------------------------------------------------------------------------
Quaternion operator* (const CPoint &rkVector, const Quaternion& rkQ)
{
    return  Quaternion(rkVector) * rkQ;
}
//----------------------------------------------------------------------------
Quaternion Quaternion::operator- () const
{
    return Quaternion(-w,-x,-y,-z);
}
//----------------------------------------------------------------------------
double Quaternion::dot (const Quaternion& rkQ) const
{
    return w*rkQ.w+x*rkQ.x+y*rkQ.y+z*rkQ.z;
}
//----------------------------------------------------------------------------
double Quaternion::norm () const
{
    return w*w+x*x+y*y+z*z;
}
//----------------------------------------------------------------------------
Quaternion Quaternion::inverse () const
{
    double fNorm = w*w+x*x+y*y+z*z;
    if ( fNorm > 0.0 )
	{
	    double fInvNorm = 1.0/fNorm;
	    return Quaternion(w*fInvNorm,-x*fInvNorm,-y*fInvNorm,-z*fInvNorm);
	}
    else
	{
	    // return an invalid result to flag the error
	    return ZERO;
	}
}
//----------------------------------------------------------------------------
Quaternion Quaternion::unitInverse () const
{
    // assert:  'this' is unit length
    return Quaternion(w,-x,-y,-z);
}
//----------------------------------------------------------------------------
Quaternion Quaternion::exp () const
{
    // If q = A*(x*i+y*j+z*k) where (x,y,z) is unit length, then
    // exp(q) = cos(A)+sin(A)*(x*i+y*j+z*k).  If sin(A) is near zero,
    // use exp(q) = cos(A)+A*(x*i+y*j+z*k) since A/sin(A) has limit 1.

    double fAngle = sqrt(x*x+y*y+z*z);
    double fSin = sin(fAngle);

    Quaternion kResult;
    kResult.w = cos(fAngle);

    if ( fabs(fSin) >= gs_fEpsilon )
	{
	    double fCoeff = fSin/fAngle;
	    kResult.x = fCoeff*x;
	    kResult.y = fCoeff*y;
	    kResult.z = fCoeff*z;
	}
    else
	{
	    kResult.x = x;
	    kResult.y = y;
	    kResult.z = z;
	}

    return kResult;
}
//----------------------------------------------------------------------------
Quaternion Quaternion::log () const
{
    // If q = cos(A)+sin(A)*(x*i+y*j+z*k) where (x,y,z) is unit length, then
    // log(q) = A*(x*i+y*j+z*k).  If sin(A) is near zero, use log(q) =
    // sin(A)*(x*i+y*j+z*k) since sin(A)/A has limit 1.

    Quaternion kResult;
    kResult.w = 0.0;

    if ( fabs(w) < 1.0 )
	{
	    double fAngle = acos(w);
	    double fSin = sin(fAngle);
	    if ( fabs(fSin) >= gs_fEpsilon )
		{
		    double fCoeff = fAngle/fSin;
		    kResult.x = fCoeff*x;
		    kResult.y = fCoeff*y;
		    kResult.z = fCoeff*z;
		    return kResult;
		}
	}

    kResult.x = x;
    kResult.y = y;
    kResult.z = z;

    return kResult;
}
//----------------------------------------------------------------------------
CPoint Quaternion::operator* (const CPoint& rkVector) const
{
    // Given a vector u = (x0,y0,z0) and a unit length quaternion
    // q = <w,x,y,z>, the vector v = (x1,y1,z1) which represents the
    // rotation of u by q is v = q*u*q^{-1} where * indicates quaternion
    // multiplication and where u is treated as the quaternion <0,x0,y0,z0>.
    // Note that q^{-1} = <w,-x,-y,-z>, so no real work is required to
    // invert q.  Now
    //
    //   q*u*q^{-1} = q*<0,x0,y0,z0>*q^{-1}
    //     = q*(x0*i+y0*j+z0*k)*q^{-1}
    //     = x0*(q*i*q^{-1})+y0*(q*j*q^{-1})+z0*(q*k*q^{-1})
    //
    // As 3-vectors, q*i*q^{-1}, q*j*q^{-1}, and 2*k*q^{-1} are the columns
    // of the rotation matrix computed in Quaternion::toRotationMatrix.
    // The vector v is obtained as the product of that rotation matrix with
    // vector u.  As such, the quaternion representation of a rotation
    // matrix requires less space than the matrix and more time to compute
    // the rotated vector.  Typical space-time tradeoff...
    NumMatrix kRot(3,3);
    toRotationMatrix(kRot);
    double tmp[3];
    tmp[0]  = kRot(0,0)*rkVector[0] + kRot(0,1)*rkVector[1] + kRot(0,2)*rkVector[2];
    tmp[1]  = kRot(1,0)*rkVector[0] + kRot(1,1)*rkVector[1] + kRot(1,2)*rkVector[2];
    tmp[2]  = kRot(2,0)*rkVector[0] + kRot(2,1)*rkVector[1] + kRot(2,2)*rkVector[2];	 
    CPoint res(DIM3,tmp);
    return res;
}
//----------------------------------------------------------------------------
Quaternion Quaternion::slerp (double fT, const Quaternion& rkP,
                              const Quaternion& rkQ)
{
    double fCos = rkP.dot(rkQ);
    double fAngle = acos(fCos);

    if ( fabs(fAngle) < gs_fEpsilon )
        return rkP;

    double fSin = sin(fAngle);
    double fInvSin = 1.0/fSin;
    double fCoeff0 = sin((1.0-fT)*fAngle)*fInvSin;
    double fCoeff1 = sin(fT*fAngle)*fInvSin;
    return fCoeff0*rkP + fCoeff1*rkQ;
}
//----------------------------------------------------------------------------
Quaternion Quaternion::slerpExtraSpins (double fT,
                                        const Quaternion& rkP, const Quaternion& rkQ, int iExtraSpins)
{
    double fCos = rkP.dot(rkQ);
    double fAngle = acos(fCos);

    if ( fabs(fAngle) < gs_fEpsilon )
        return rkP;

    double fSin = sin(fAngle);
    double fPhase = M_PI*iExtraSpins*fT;
    double fInvSin = 1.0/fSin;
    double fCoeff0 = sin((1.0-fT)*fAngle - fPhase)*fInvSin;
    double fCoeff1 = sin(fT*fAngle + fPhase)*fInvSin;
    return fCoeff0*rkP + fCoeff1*rkQ;
}
//----------------------------------------------------------------------------
void Quaternion::intermediate (const Quaternion& rkQ0,
                               const Quaternion& rkQ1, const Quaternion& rkQ2, Quaternion& rkA,
                               Quaternion& rkB)
{
    // assert:  q0, q1, q2 are unit quaternions

    Quaternion kQ0inv = rkQ0.unitInverse();
    Quaternion kQ1inv = rkQ1.unitInverse();
    Quaternion rkP0 = kQ0inv*rkQ1;
    Quaternion rkP1 = kQ1inv*rkQ2;
    Quaternion kArg = 0.25*(rkP0.log()-rkP1.log());
    Quaternion kMinusArg = -kArg;

    rkA = rkQ1*kArg.exp();
    rkB = rkQ1*kMinusArg.exp();
}
//----------------------------------------------------------------------------
Quaternion Quaternion::squad (double fT, const Quaternion& rkP,
                              const Quaternion& rkA,	 const Quaternion& rkB, const Quaternion& rkQ)
{
    double fSlerpT = 2.0*fT*(1.0-fT);
    Quaternion kSlerpP = slerp(fT,rkP,rkQ);
    Quaternion kSlerpQ = slerp(fT,rkA,rkB);
    return slerp(fSlerpT,kSlerpP,kSlerpQ);
}
//----------------------------------------------------------------------------
