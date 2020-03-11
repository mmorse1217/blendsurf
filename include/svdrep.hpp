/*! \file */
#ifndef _SVDREP_HPP_
#define _SVDREP_HPP_

#include "nummatrix.hpp"

class SVDRep
{
public:
    SVDRep()   {}
    SVDRep(const SVDRep& C): _matU(C._matU), _matS(C._matS), _matVT(C._matVT)  {}
    ~SVDRep()  {}

    SVDRep& operator=(const SVDRep& c)  { _matU = c._matU; _matS = c._matS; _matVT = c._matVT; return *this; }
    //access
    NumMatrix& U() { return _matU; }
    NumVector& S() { return _matS; }
    NumMatrix& VT(){ return _matVT; }
    //ops
    int construct(double epsilon, const NumMatrix& M);

    //  int dgemv(double alpha, const NumVector& X, double beta, NumVector& Y, double tol=0.0); // y <- a Mx + b y
    //  int dgemv(double alpha, double* X, double beta, double* Y, double tol=0.0);
  
    int m() const { return _matU.m(); }
    int k() const { return _matS.length(); }
    int n() const { return _matVT.n(); }
  
protected:
    NumMatrix _matU;
    NumVector _matS;
    NumMatrix _matVT;
    static int _wssize;
    static double _wsbuf[];
    //static int _wisize;
    //static int _wibuf[];
};	

inline ostream& operator<<( ostream& os, SVDRep& svdrep)
{
    os<<svdrep.U().m()<<" "<<svdrep.S().length()<<" "<<svdrep.VT().n()<<endl;
    os<<svdrep.U()<<svdrep.S()<<svdrep.VT()<<endl;
    return os;
}

#endif
