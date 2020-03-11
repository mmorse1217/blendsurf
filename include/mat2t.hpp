#ifndef _MAT2T_HPP_
#define _MAT2T_HPP_

#include "vec2t.hpp"
#include "vecmatop.hpp"

template<class F> 
class Mat2T {
private:
    F _v[4];
public:
    enum{ X=0, Y=1 };
    //--------constructors and destructors
    Mat2T() {
        for(int i=0; i<2*2; i++)		_v[i] = 0;
    }
    Mat2T(F f) { 
        for(int i=0; i<2*2; i++)		_v[i] = f;
    }
    Mat2T(const F* m) { 
        for(int i=0; i<2*2; i++)		_v[i] = m[i];
    }
    Mat2T(const Mat2T& m)  {
        for(int i=0; i<2*2; i++)		_v[i] = m._v[i];
    }
    Mat2T(const Vec2T<F>& v0, const Vec2T<F>& v1) { //column by column
        (*this)(0,0) = v0(0);	 (*this)(0,1) = v1(0);
        (*this)(1,0) = v0(1);	 (*this)(1,1) = v1(1);
    }
    //-------pointer
    operator F*(void)       { return &_v[0]; }
    operator const F*(void) { return &_v[0]; }
    F* array()              { return &_v[0]; }
    const F* array() const  { return &_v[0]; }
    F& operator () (int i, int j)             { return _v[i + j*2]; }
    const F& operator () (int i, int j) const { return _v[i + j*2]; }
    //------assign
    Mat2T& operator=  (const Mat2T& m) {
        for(int i=0; i<2*2; i++)
            _v[i] = m._v[i];
        return *this;
    }
    Mat2T& operator+= (const Mat2T& m) {
        for(int i=0; i<2*2; i++)
            _v[i] += m._v[i];
        return *this;
    }
    Mat2T& operator-= (const Mat2T& m)  {
        for(int i=0; i<2*2; i++)
            _v[i] -= m._v[i];
        return *this;
    }
    Mat2T& operator*= (F d)  {
        for(int i=0; i<2*2; i++)
            _v[i] *= d;
        return *this;
    }
    Mat2T& operator/= (F d)  {
        for(int i=0; i<2*2; i++)
            _v[i] /= d;
        return *this;
    }
    //------------norm stuff
  
    //------------inverse
    Mat2T inverse();                                            // inverse
};

//-----------NUMERICAL OPS
template<class F> inline Mat2T<F> operator- (const Mat2T<F>& a) {
    Mat2T<F> r;  for(int i=0; i<2; i++)	 for(int j=0; j<2; j++)		r(i,j) = -a(i,j);
    return r;
}
template<class F> inline Mat2T<F> operator+ (const Mat2T<F>& a, const Mat2T<F>& b) {
    Mat2T<F> r;  for(int i=0; i<2; i++)	 for(int j=0; j<2; j++)		r(i,j) = a(i,j)+b(i,j);
    return r;
}
template<class F> inline Mat2T<F> operator- (const Mat2T<F>& a, const Mat2T<F>& b) {
    Mat2T<F> r;  for(int i=0; i<2; i++)	 for(int j=0; j<2; j++)		r(i,j) = a(i,j)-b(i,j);
    return r;
}
template<class F> inline Mat2T<F> operator* (F scl, const Mat2T<F>& a) {
    Mat2T<F> r;  for(int i=0; i<2; i++)	 for(int j=0; j<2; j++)		r(i,j) = scl*a(i,j);
    return r;
}
template<class F> inline Mat2T<F> operator* (const Mat2T<F>& a, F scl) {
    Mat2T<F> r;  for(int i=0; i<2; i++)	 for(int j=0; j<2; j++)		r(i,j) = scl*a(i,j);
    return r;
}
template<class F> inline Mat2T<F> operator/ (const Mat2T<F>& a, F scl) {
    Mat2T<F> r;  for(int i=0; i<2; i++)	 for(int j=0; j<2; j++)		r(i,j) = a(i,j) / scl;
    return r;
}
template<class F> inline Vec2T<F> operator* (const Mat2T<F>& m, const Vec2T<F>& v) { 
    Vec2T<F> tmp;
    for(int i=0; i<2; i++)
        for(int j=0; j<2; j++)
            tmp(i) += m(i,j)*v(j);
    return tmp;
}
template<class F> inline Mat2T<F> operator* (const Mat2T<F>& a, const Mat2T<F>& b) {
    Mat2T<F> r;
    for(int i=0; i<2; i++)
        for(int j=0; j<2; j++)
            for(int k=0; k<2; k++)
                r(i,j) += a(i,k) * b(k,j);
    return r;
}

//---------------INOUT
template<class F> istream& operator >> (istream& s, Mat2T<F>& m) { 
    for(int i=0; i<2; i++)
        for(int j=0; j<2; j++)
            s>>m(i,j);
    return s;
}
template<class F> ostream& operator <<(ostream& s, const Mat2T<F>& m) { 
    for(int i=0; i<2; i++) {
        for(int j=0; j<2; j++)
            s<<m(i,j)<<" ";
        s<<endl;
    }
    return s;
}

//---------------INVERSE
template<class F> Mat2T<F> Mat2T<F>::inverse()
{
    if(       sizeof(F)==4) {//float
        assert(0);
    } else if(sizeof(F)==8) {//double
        Mat2T<F> res;
        NumMatrix in(2,2,false, this->array());
        NumMatrix ou(2,2,false, res.array());
        inv(in, ou); //out matrix
        return res;
    } else {
        assert(0);
    }
}


//---------------------------------------------------------
/// MOST COMMONLY USED
typedef Mat2T<double> Matrix2;

#endif
