#ifndef  _CVECT_HPP_
#define  _CVECT_HPP_

#include "common.hpp"


///Common VECtor Template
template <class F>
class CVecT {
private:
    DimType _dim; //dimension;
    F _v[MAXDIM];
public:
    //------------CONSTRUCTOR AND DESTRUCTOR 
    //CVecT()                      { _dim = DIM0; } //invalid state
    CVecT(DimType d)             { _dim = d;      for(int i=0; i<_dim; i++) _v[i]=F(0); }
    CVecT(DimType d, F f)        { _dim = d;      for(int i=0; i<_dim; i++) _v[i]=f; }
    CVecT(DimType d, const F* f) { _dim = d;      for(int i=0; i<_dim; i++) _v[i]=f[i]; }
    //CVecT(DimType d, F a,F b)    { assert(d==DIM2); _dim = d;    _v[0]=a; _v[1]=b; }
    //CVecT(           F a,F b)    {                  _dim = DIM2; _v[0]=a; _v[1]=b; }
    //CVecT(DimType d, F a,F b,F c){ assert(d==DIM3); _dim = d;    _v[0]=a; _v[1]=b; _v[2]=c; }
    //CVecT(           F a,F b,F c){                  _dim = DIM3; _v[0]=a; _v[1]=b; _v[2]=c; }
    CVecT(const CVecT& c)        { _dim = c._dim; for(int i=0; i<_dim; i++) _v[i]=c._v[i]; }
    ~CVecT() {}
  
    //------------POINTER and ACCESS
    DimType dim() const       { return _dim; }
    operator F*()             { return &_v[0]; }
    operator const F*() const { return &_v[0]; }
    F* array()                { return &_v[0]; }  //access array
    F& operator()(int i)             { assert(i<_dim); return _v[i]; }
    const F& operator()(int i) const { assert(i<_dim); return _v[i]; }
    F& operator[](int i)             { assert(i<_dim); return _v[i]; }
    const F& operator[](int i) const { assert(i<_dim); return _v[i]; }
    F& x()             { return _v[0];}
    F& y()             { return _v[1];}
    F& z()             { return _v[2];}
    const F& x() const { return _v[0];}
    const F& y() const { return _v[1];}
    const F& z() const { return _v[2];}
  
    //------------ASSIGN
    CVecT& operator= ( const CVecT& c ) { _dim = c._dim;        for(int i=0; i<_dim; i++) _v[i] =c._v[i]; return *this; }
    CVecT& operator+=( const CVecT& c ) { assert(_dim==c._dim); for(int i=0; i<_dim; i++) _v[i]+=c._v[i]; return *this; }
    CVecT& operator-=( const CVecT& c ) { assert(_dim==c._dim); for(int i=0; i<_dim; i++) _v[i]-=c._v[i]; return *this; }
    CVecT& operator*=( const F& s )     {                       for(int i=0; i<_dim; i++) _v[i]*=s;       return *this; }
    CVecT& operator/=( const F& s )     {                       for(int i=0; i<_dim; i++) _v[i]/=s;       return *this; }
  
    //-----------LENGTH...
    F l1( void )     const  { F sum=F(0); for(int i=0; i<_dim; i++) sum=sum+abs(_v[i]); return sum; }
    F linfty( void ) const  { F cur=F(0); for(int i=0; i<_dim; i++) cur=max(cur,abs(_v[i])); return cur; }
    F l2( void )     const  { F sum=F(0); for(int i=0; i<_dim; i++) sum=sum+_v[i]*_v[i]; return sqrt(sum); }
    F length( void ) const  { return l2(); }
    CVecT dir( void )    const  { F a=l2(); return (*this)/a; }
};

//-----------BOOLEAN OPS
template <class F> inline bool operator==(const CVecT<F>& a, const CVecT<F>& b) {
    assert(a.dim()==b.dim());
    bool res = true;  for(int i=0; i<a.dim(); i++)   res = res && (a(i)==b(i));  return res;
}
template <class F> inline bool operator!=(const CVecT<F>& a, const CVecT<F>& b) {
    return !(a==b);
}
template <class F> inline bool operator> (const CVecT<F>& a, const CVecT<F>& b) {
    assert(a.dim()==b.dim());  
    bool res = true;  for(int i=0; i<a.dim(); i++)   res = res && (a(i)> b(i));  return res; 
}
template <class F> inline bool operator< (const CVecT<F>& a, const CVecT<F>& b) {
    assert(a.dim()==b.dim());  
    bool res = true;  for(int i=0; i<a.dim(); i++)   res = res && (a(i)< b(i));  return res; 
}
template <class F> inline bool operator>=(const CVecT<F>& a, const CVecT<F>& b) {
    assert(a.dim()==b.dim());
    bool res = true;  for(int i=0; i<a.dim(); i++)	res = res && (a(i)>=b(i));  return res; 
}
template <class F> inline bool operator<=(const CVecT<F>& a, const CVecT<F>& b) {
    assert(a.dim()==b.dim());
    bool res = true;  for(int i=0; i<a.dim(); i++)   res = res && (a(i)<=b(i));  return res; 
}

//-----------NUMERICAL OPS
template <class F> inline CVecT<F> operator- (const CVecT<F>& a) {
    CVecT<F> r(a.dim());  for(int i=0; i<a.dim(); i++) r[i] = -a[i]; return r;
}
template <class F> inline CVecT<F> operator+ (const CVecT<F>& a, const CVecT<F>& b) {
    assert(a.dim()==b.dim()); 
    CVecT<F> r(a.dim());  for(int i=0; i<a.dim(); i++) r[i] = a[i]+b[i]; return r; 
}
template <class F> inline CVecT<F> operator- (const CVecT<F>& a, const CVecT<F>& b) {
    assert(a.dim()==b.dim()); 
    CVecT<F> r(a.dim());  for(int i=0; i<a.dim(); i++) r[i] = a[i]-b[i]; return r;
}
template <class F> inline CVecT<F> operator* (F scl, const CVecT<F>& a) {
    CVecT<F> r(a.dim());  for(int i=0; i<a.dim(); i++) r[i] = scl*a[i];  return r;
}
template <class F> inline CVecT<F> operator* (const CVecT<F>& a, F scl) {
    CVecT<F> r(a.dim());  for(int i=0; i<a.dim(); i++) r[i] = scl*a[i];  return r;
}
template <class F> inline CVecT<F> operator/ (const CVecT<F>& a, F scl) {
    CVecT<F> r(a.dim());  for(int i=0; i<a.dim(); i++) r[i] = a[i]/scl;  return r;
}

template <class F> inline F operator* (const CVecT<F>& a, const CVecT<F>& b) {
    assert(a.dim()==b.dim());
    F sum=F(0); for(int i=0; i<a.dim(); i++) sum=sum+a(i)*b(i); return sum;
}
template <class F> inline F dot       (const CVecT<F>& a, const CVecT<F>& b) {
    return a*b;
}
template <class F> inline CVecT<F> operator^ (const CVecT<F>& a, const CVecT<F>& b) {
    assert(a.dim()==DIM3 && b.dim()==DIM3);
    return CVecT<F>(DIM3, a(1)*b(2)-a(2)*b(1), a(2)*b(0)-a(0)*b(2), a(0)*b(1)-a(1)*b(0)); 
}
template <class F> inline CVecT<F> cross     (const CVecT<F>& a, const CVecT<F>& b) { 
    return a^b; 
}

//-------------ew OPS
template <class F> inline CVecT<F> min(const CVecT<F>& a, const CVecT<F>& b) {
    assert(a.dim()==b.dim()); 
    CVecT<F> r(a.dim());  for(int i=0; i<a.dim(); i++) r[i] = min(a[i], b[i]); return r;
}
template <class F> inline CVecT<F> max(const CVecT<F>& a, const CVecT<F>& b) {
    assert(a.dim()==b.dim());
    CVecT<F> r(a.dim());  for(int i=0; i<a.dim(); i++) r[i] = max(a[i], b[i]); return r;
}
template <class F> inline CVecT<F> abs(const CVecT<F>& a) {
    CVecT<F> r(a.dim());  for(int i=0; i<a.dim(); i++) r[i] = abs(a[i]); return r;
}
template <class F> inline CVecT<F> ewmul(const CVecT<F>&a, const CVecT<F>& b) {
    assert(a.dim()==b.dim());
    CVecT<F> r(a.dim());  for(int i=0; i<a.dim(); i++) r[i] = a[i]*b[i]; return r;
}
template <class F> inline CVecT<F> ewdiv(const CVecT<F>&a, const CVecT<F>& b) { 
    assert(a.dim()==b.dim()); 
    CVecT<F> r(a.dim());  for(int i=0; i<a.dim(); i++) r[i] = a[i]/b[i]; return r;
}

//---------------SHIFT
template <class F> inline CVecT<F> operator<<(const CVecT<F>& a, int num) {
    CVecT<F> r(a.dim());  for(int i=0; i<a.dim(); i++) r[i] = a[i]<<num; return r;
}
template <class F> inline CVecT<F> operator>>(const CVecT<F>& a, int num) {
    CVecT<F> r(a.dim());  for(int i=0; i<a.dim(); i++) r[i] = a[i]>>num; return r;
}
//---------------INOUT
template <class F> istream& operator>>(istream& is, CVecT<F>& a) {
    for(int i=0; i<a.dim(); i++) is>>a[i]; return is;
}
template <class F> ostream& operator<<(ostream& os, const CVecT<F>& a) { 
    for(int i=0; i<a.dim(); i++) os<<a[i]<<" "; return os;
}

//---------------------------------------------------------
/// MOST COMMONLY USED
typedef CVecT<double> CPoint;
typedef CVecT<int>    CIndex;

//--------------------------------------------------------------------------------------

/// Convert a CIndex to CPoint
inline CPoint ci2cp(CIndex x) {
    CPoint r(x.dim()); for(int i=0; i<x.dim(); i++) r[i]=double(x[i]); return r;  
}

inline  void CPointToDouble(const CPoint &p, double *arr){
    for(int i=0; i<p.dim(); i++) arr[i]=p[i];  //arr[0] = p[0];  //arr[1] = p[1];  //if(p.dim()==DIM3) arr[2]=p[2];
}
inline  void CPointToDoubleAdd(const CPoint &p, double *arr){
    for(int i=0; i<p.dim(); i++) arr[i]+=p[i];  //arr[0] = p.x(); arr[1] = p.y();  if(p.dim()==DIM3) arr[2]+=p.z();
}
inline  CPoint doubleToCPoint(double *arr, DimType d){
    return CPoint(d, arr);
}

/// \todo GB: SUN compiler confuses this abs with the abs on ebi_uti.hpp. 
///The old version doesn't work in 2D (because F(0) in abs of ebi_uti.hpp defaults 3D)

#endif
