/*! \file */
#ifndef _VECMATOP_HPP_
#define _VECMATOP_HPP_

#include "nummatrix.hpp"
#include "vec3t.hpp"

void matvecmult3(const NumMatrix& A, const vector<Point3>& v, vector<Point3>& ret );

//--------------------------------------------------
// matrix operation, call lapack



// c = alpha*a*b + beta*c
int dgemm(double alpha, const NumMatrix& A, const NumMatrix& B, double beta, NumMatrix& C);
int dgemm(int m, int n, int k, double alpha, double* A, double* B, double beta, double* C);
// R <= tran(M)
int tran(const NumMatrix& M, NumMatrix& R);
// R <= pinv(M, epsilon)
int pinv(const NumMatrix& M, double epsilon, NumMatrix& R);
// R <= inv(M);
int inv(const NumMatrix& M, NumMatrix& R);
//--------------------------------------------------
// interpolation, etc.
//evaluation flags
enum {  EVFLAG_VL = 1,  EVFLAG_FD = 2,  EVFLAG_SD = 4 };
//domain flag
enum {  DMFLAG_PERIOD = 0,  DMFLAG_CLOSED = 1 };

// cubic spline interpolation
int spev1d( int evflag, int dmflag, int dof, double* M, int n,   double e,   int i,   double u,   double* res);
int spev2d( int evflag, int dmflag, int dof, double* M, int* mn, double* ef, int* ij, double* uv, double* res);
int spcoef( int evflag, double u, double* us);

#endif
