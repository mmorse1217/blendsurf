#ifndef _MAT3T_HPP_
#define _MAT3T_HPP_

#include "vec3t.hpp"
#include "vecmatop.hpp"

template<class F> 
class Mat3T {
private:
    F _v[9];
public:
    enum{ X = 0, Y = 1, Z = 2 };
    //--------constructors and destructors
    Mat3T() { ; }
    Mat3T(F f) { 
        for(int i=0; i<3*3; i++)
            _v[i] = f;
    }
    Mat3T(const F* m) { 
        for(int i=0; i<3*3; i++)
            _v[i] = m[i];
    }
    Mat3T(const Mat3T& m)  {
        for(int i=0; i<3*3; i++)
            _v[i] = m._v[i];
    }
    Mat3T(const Vec3T<F>& v0, const Vec3T<F>& v1, const Vec3T<F>& v2) { //column by column
        (*this)(0,0) = v0(0);	 (*this)(0,1) = v1(0);	 (*this)(0,2) = v2(0);
        (*this)(1,0) = v0(1);	 (*this)(1,1) = v1(1);	 (*this)(1,2) = v2(1);
        (*this)(2,0) = v0(2);	 (*this)(2,1) = v1(2);	 (*this)(2,2) = v2(2);
    }
    //-------pointer
    operator F*(void)       { return &_v[0]; }
    operator const F*(void) { return &_v[0]; }
    F* array()              { return &_v[0]; }
    const F* array() const  { return &_v[0]; }
    F& operator () (int i, int j)             { return _v[i + j*3]; }
    const F& operator () (int i, int j) const { return _v[i + j*3]; }
    //------assign
    Mat3T& operator=  (const Mat3T& m) {
        for(int i=0; i<3*3; i++)
            _v[i] = m._v[i];
        return *this;
    }
    Mat3T& operator+= (const Mat3T& m) {
        for(int i=0; i<3*3; i++)
            _v[i] += m._v[i];
        return *this;
    }
    Mat3T& operator-= (const Mat3T& m)  {
        for(int i=0; i<3*3; i++)
            _v[i] -= m._v[i];
        return *this;
    }
    Mat3T& operator*= (F d)  {
        for(int i=0; i<3*3; i++)
            _v[i] *= d;
        return *this;
    }
    Mat3T& operator/= (F d)  {
        for(int i=0; i<3*3; i++)
            _v[i] /= d;
        return *this;
    }
    //------------norm stuff
  
    //------------inverse
    Mat3T inverse();                                            // inverse
};

//-----------NUMERICAL OPS
template<class F> inline Vec3T<F> operator * (const Mat3T<F>& m, const Vec3T<F>& v) { 
    Vec3T<F> tmp(0.0);
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            tmp(i) += m(i,j)*v(j);
    return tmp;
}

//---------------INOUT
template<class F> istream& operator >> (istream& s, Mat3T<F>& m) { 
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            s>>m(i,j);
    return s;
}
template<class F> ostream& operator <<(ostream& s, const Mat3T<F>& m) { 
    for(int i=0; i<3; i++) {
        for(int j=0; j<3; j++)
            s<<m(i,j)<<" ";
        s<<endl;
    }
    return s;
}

//---------------INVERSE
template<class F> Mat3T<F> Mat3T<F>::inverse()
{
    if(       sizeof(F)==4) {//float
        assert(0);
    } else if(sizeof(F)==8) {//double
        Mat3T<F> res;
        NumMatrix in(3,3,false, this->array());
        NumMatrix ou(3,3,false, res.array());
        inv(in, ou); //out matrix
        return res;
    } else {
        assert(0);
    }
}

//---------------------------------------------------------
/// MOST COMMONLY USED
typedef Mat3T<double> Matrix3;


#endif
