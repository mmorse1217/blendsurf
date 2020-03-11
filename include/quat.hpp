// $Id$
// $Source$

#ifndef _QUAT_HPP_
#define _QUAT_HPP_

#include "common.hpp"
#include "vec3t.hpp"


class HMatrix;
class Quat;

class Quat 
{ 
public:
    friend class HMatrix;
    Quat() { }
    Quat(float xc,float yc,float zc, float wc): 
            x(xc), y(yc), z(zc), w(wc) {}
    Quat( const Quat& q ) : x(q.x), y(q.y), z(q.z), w(q.w) {}
    // Construct a unit quaternion from two points on unit sphere
    Quat(const Vec3T<float>& from, const Vec3T<float>& to) {
        x = from.y()*to.z() - from.z()*to.y();
        y = from.z()*to.x() - from.x()*to.z();
        z = from.x()*to.y() - from.y()*to.x();
        w = from.x()*to.x() + from.y()*to.y() + from.z()*to.z();
    }
    // Convert a unit quaternion to two points on unit sphere 
    void ballPoints(Vec3T<float>& arcFrom, Vec3T<float>& arcTo) {
        double s;
        s = sqrt(x*x + y*y);
        arcFrom =  (s == 0.0)? Vec3T<float>(0.0f, 1.0f, 0.0f): Vec3T<float>(-y/s,x/s, 0.0f);
        if (w < 0.0f) arcFrom = -arcFrom; 
        arcTo = Vec3T<float>(w*arcFrom.x() - z*arcFrom.y(),
                             w*arcFrom.y() + z*arcFrom.x(),
                             x*arcFrom.y() - y*arcFrom.x());
    }
    Quat& operator=(const Quat& q) { 
        x = q.x; y = q.y; z = q.z;  w = q.w; return *this; 
    } 
    Quat operator-() const { 
        return Quat(-x,-y,-z,-w); 
    }
    Quat operator+(const Quat& q) const {
        return Quat(x + q.x,y + q.y,z + q.z,w + q.w);  
    } 
    Quat operator-(const Quat& q) { 
        return Quat(x - q.x, y - q.y, z - q.z, w - q.w);  
    } 
    Quat operator*( const float s ) const {
        return Quat(s*x,s*y,s*z,s*w);
    }
  
    //To combine rotations, use the product qSecond*qFirst,
    //which gives the effect of rotating by qFirst then qSecond
    Quat operator*(const Quat& q)  {
        return Quat(
		    w*q.x + x*q.w + y*q.z - z*q.y,
		    w*q.y + y*q.w + z*q.x - x*q.z,
		    w*q.z + z*q.w + x*q.y - y*q.x,
		    w*q.w - x*q.x - y*q.y - z*q.z
		    );
    }
    Quat conj() { return Quat(-x,-y,-z,w); } 
  
private:
    float x,y,z,w;
};

inline  Quat operator*(const float s, const Quat& q)  { return q*s;  } 

#endif  /* __QUAT_H__ */


