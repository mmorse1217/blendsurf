/*! \file */
#ifndef _NUMMATRIX_HPP_
#define _NUMMATRIX_HPP_

#include "numvector.hpp"

template <class F>
class Matrix
{
public:
    int _m, _n;
    bool _owndata;
    F* _data;
public:
    Matrix(int m=0, int n=0): _m(m), _n(n), _owndata(true) {
        if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); } else _data=NULL;
    }
    Matrix(int m, int n, bool owndata, F* data): _m(m), _n(n), _owndata(owndata) {
        if(_owndata) {
            if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); } else _data=NULL;
            if(_m>0 && _n>0) memcpy( _data, data, _m*_n*sizeof(F) );
        } else {
            _data = data;
        }
    }
    Matrix(const Matrix& C): _m(C._m), _n(C._n), _owndata(C._owndata) {
        if(_owndata) {
            if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); } else _data=NULL;
            if(_m>0 && _n>0) memcpy( _data, C._data, _m*_n*sizeof(F) );
        } else {
            _data = C._data;
        }
    }
    ~Matrix() { 
        if(_owndata) { 
            if(_m>0 && _n>0) { delete[] _data; _data = NULL; } 
        }
    }
    Matrix& operator=(const Matrix& C) {
        if(_owndata) { 
            if(_m>0 && _n>0) { delete[] _data; _data = NULL; } 
        }
        _m = C._m; _n=C._n; _owndata=C._owndata;
        if(_owndata) {
            if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); } else _data=NULL;
            if(_m>0 && _n>0) memcpy( _data, C._data, _m*_n*sizeof(F) );
        } else {
            _data = C._data;
        }
        return *this;
    }
    void resize(int m, int n)  {
        //	 assert( _owndata==true );
        if(_m!=m || _n!=n) {
            if(_m>0 && _n>0) { delete[] _data; _data = NULL; } 
            _m = m; _n = n;
            if(_m>0 && _n>0) { _data = new F[_m*_n]; assert( _data!=NULL ); } else _data=NULL;
        }
    }
    const F& operator()(int i, int j) const  { 
        assert( i>=0 && i<_m && j>=0 && j<_n );
        return _data[i+j*_m];
    }
    F& operator()(int i, int j)  { 
        assert( i>=0 && i<_m && j>=0 && j<_n );
        return _data[i+j*_m];
    }
  
    void getColumn(int j, Vector<F>& vec)  {
        assert( j>=0 && j<n() );
        vec.resize(m());
        for(int i=0; i<m(); i++)
            vec(i) = (*this)(i,j);
    }
    void getRow(int i, Vector<F>& vec)  {
        assert( i>=0 && i<m() );
        vec.resize(n());
        for(int j=0; j<n(); j++)
            vec(j) = (*this)(i,j);
    }
    void setColumn(int j, Vector<F>& vec)  {
        assert( j>=0 && j<n() );
        assert( vec.length() == m() );
        for(int i=0; i<m(); i++)
            (*this)(i,j) = vec(i);
    }
    void setRow(int i, Vector<F>& vec)  {
        assert( i>=0 && i<m());
        assert( vec.length() == n());
        for(int j=0; j<n(); j++)
            (*this)(i,j) = vec(j);
    }
  
    F* data() const { return _data; }
    F* clmdata(int j) { return &(_data[j*_m]); }
    int m() const { return _m; }
    int n() const { return _n; }
};

template <class F> inline ostream& operator<<( ostream& os, const Matrix<F>& mat)
{
    os<<mat.m()<<" "<<mat.n()<<endl;
    for(int i=0; i<mat.m(); i++) {
        for(int j=0; j<mat.n(); j++)
            os<<" "<<mat(i,j);
        os<<endl;
    }
    return os;
}
template <class F> void setvalue(Matrix<F>& M, F val)
{
    for(int i=0; i<M.m(); i++)
        for(int j=0; j<M.n(); j++)
            M(i,j) = val;
}
template <class F> inline void multvalue(Matrix<F>& M, F val)
{
    for(int i=0; i<M.m(); i++)
        for(int j=0; j<M.n(); j++)
            M(i,j) *= val;
}
template <class F> inline void clear(Matrix<F>& M)
{
    memset(M.data(), 0, M.m()*M.n()*sizeof(F));
}

typedef Matrix<bool>   BolMatrix;
typedef Matrix<int>    IntMatrix;
typedef Matrix<double> NumMatrix; //typedef Matrix<double> SclMatrix;
//typedef Matrix<double> DblNumMat; //typedef Matrix<double> SclMatrix;

/* ********************************************************************** */
inline void NumMatrixToDouble( const NumMatrix &In, double *ar){
    if(ar==NULL) return;
    memcpy( ar, In.data(), In.m()*In.n()*sizeof(double) );
}

/* ********************************************************************** */
inline void doubleToNumMatrix( const double *ar, NumMatrix &In){
    if(ar==NULL) return;
    memcpy( In._data, ar, In.m()*In.n()*sizeof(double) );
}

/// Matmat
inline void NumMatrixMult(const NumMatrix &M1, const NumMatrix &M2, NumMatrix &M)
{
    assert( M1.n()==M2.m() );
    assert( M.m() == M1.m() && M.n() == M2.n());
    for( int i=0; i<M1.m();i++){
        for( int j=0; j<M2.n(); j++){
            M(i,j) = 0;
            for( int k=0;k<M1.n();k++) M(i,j) += M1(i,k)*M2(k,j);
        }
    }
}

/* ********************************************************************** */
inline void NumMatrixAdd( NumMatrix &M1, NumMatrix &M2, NumMatrix &M)
{
    assert( M1.n()==M2.n()&& M1.n()==M.n() && M1.m()==M2.m() && M1.m() == M.m() );	
    for( int i=0; i<M1.m(); i++){
        for( int j=0; j<M2.n(); j++) M(i,j) = M1(i,j) + M2(i,j);
    }
}
/* ********************************************************************** */
/// In place transpose for square matrices
inline void NumMatrixTranspose( NumMatrix &M)
{
    assert( M.n() == M.m()); // need implementation for non-square matrices 
        double swap;
        for( int i=0;i<M.n();i++){
            for( int j=0; j<i; j++){
                swap = M(i,j);
                M(i,j) = M(j,i);
                M(j,i) = swap;
            }
        }
    }
inline void NumMatrixTranspose( NumMatrix M, NumMatrix& M_T)
{
    //NumMatrix M_T(M.n(), M.m());
    for( int i=0;i<M.n();i++)
        for( int j=0; j<M.m(); j++)
            M_T(i,j) = M(j,i);
    //return M_T;
}

inline void NumMatrixScale(double scale, NumMatrix& M){
    for(int i=0; i<M.n(); i++){
        for( int j=0; j<M.m(); j++){
            M(i,j) *= scale;
        }
    }
}



#endif




