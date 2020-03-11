/*! \file */
#ifndef _NUMVECTOR_HPP_
#define _NUMVECTOR_HPP_


#include "common.hpp"
template <class F>
class Vector
{
public:
    int  _length;
    bool _owndata;
    F* _data;
public:
    Vector(int length=0): _length(length), _owndata(true)  {
        if(_length>0) { _data = new F[_length]; assert(_data!=NULL); } else _data=NULL;
    }
    Vector(int length, bool owndata, F* data): _length(length), _owndata(owndata) {
        if(_owndata) {
            if(_length>0) { _data = new F[_length]; assert(_data!=NULL); } else _data=NULL;
            if(_length>0) memcpy( _data, data, _length*sizeof(F) );
        } else {
            _data = data;
        }
    }
    Vector(const Vector& C): _length(C._length), _owndata(C._owndata)  {
        if(_owndata) {
            if(_length>0) { _data = new F[_length]; assert(_data!=NULL); } else _data=NULL;
            if(_length>0) memcpy( _data, C._data, _length*sizeof(F) );
        } else {
            _data = C._data;
        }
    }
    ~Vector() {
        if(_owndata) {
            if(_length>0) { delete[] _data; _data = NULL; }
        }
    }
    Vector& operator=(const Vector& C)  {
        if(_owndata) { 
            if(_length>0) { delete[] _data; _data = NULL; }
        }
        _length = C._length; _owndata=C._owndata;
        if(_owndata) {
            if(_length>0) { _data = new F[_length]; assert(_data!=NULL); } else _data=NULL;
            if(_length>0) memcpy( _data, C._data, _length*sizeof(F) );
        } else {
            _data =C._data;
        }
        return *this;
    }
    void resize(int length)  {
        //	 assert(_owndata==true);
        if(length !=_length) {
            if(_length>0) { delete[] _data; _data = NULL; }
            _length = length;
            if(_length>0) { _data = new F[_length]; assert(_data!=NULL); } else _data=NULL;
        }
    }
    const F& operator()(int i) const  {
        assert(i>=0 && i<_length);
        return _data[i]; 
    }
    F& operator()(int i)  {
        assert(i>=0 && i<_length);
        return _data[i]; 
    }
  
    F* data() const { return _data; }
    int length () const { return _length; }
};

template <class F> inline ostream& operator<<( ostream& os, const Vector<F>& vec)
{
    os<<vec.length()<<endl;
    for(int i=0; i<vec.length(); i++)
        os<<" "<<vec(i);
    os<<endl;
    return os;
}
template <class F> inline void setvalue(Vector<F>& V, F val)
{
    for(int i=0; i<V.length(); i++)
        V(i) = val;
}
template <class F> inline void clear(Vector<F>& V)
{
    memset(V.data(), 0, V.length()*sizeof(F));
}
template <class F> 
inline double dot(Vector<F> V1, Vector<F> V2){
    assert(V1.length() == V2.length());
    double dot_product = 0.;
    for(int i =0; i < V1.length(); i++){
        dot_product += V1(i)*V2(i);
    }
    return dot_product;
}


typedef Vector<bool>   BolVector;
typedef Vector<int>    IntVector;
typedef Vector<double> NumVector; //typedef Vector<double> SclVector;

#endif


