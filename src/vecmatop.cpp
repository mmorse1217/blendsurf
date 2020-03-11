/*! \file */
#include "blas_lapack.h"
#include "svdrep.hpp"
#include "vecmatop.hpp"

// ---------------------------------------------------------------------- 
void matvecmult3(const NumMatrix& A, const vector<Point3>& v, vector<Point3>& ret )
{
  assert(A.n() == v.size()); assert(A.m() == ret.size());

  for(int i=0; i<A.m(); i++){
    ret[i] = Point3(0.0);
    for(int j=0; j<A.n(); j++)
      ret[i] += A(i,j) *v[j];
  }

}



//Y <- a M X + b Y
// ---------------------------------------------------------------------- 
int dgemm(double alpha, const NumMatrix& A, const NumMatrix& B, double beta, NumMatrix& C)
{
  
    assert( A.m() == C.m() );  assert( A.n() == B.m() );  assert( B.n() == C.n() );
    dgemm(C.m(), C.n(), A.n(), alpha, A.data(), B.data(), beta, C.data());
    return 0;
}
// ---------------------------------------------------------------------- 
int dgemm(int m, int n, int k, double alpha, double* A, double* B, double beta, double* C)
{
  
    char transa = 'N';
    char transb = 'N';
    assert(m!=0 && n!=0 && k!=0);
    DGEMM(&transa, &transb, &m, &n, &k,
          &alpha, A, &m, B, &k, &beta, C, &m);
    return 0;
}
// ---------------------------------------------------------------------- 
int tran(const NumMatrix& M, NumMatrix& R)
{
  
    assert(R.m()==M.n() && R.n()==M.m());  //R.resize(M.n(), M.m());
    for(int i=0; i<M.m(); i++)
        for(int j=0; j<M.n(); j++)
            R(j,i) = M(i,j);
    return 0;
}
// ----------------------------------------------------------------------
int pinv(const NumMatrix& M, double epsilon, NumMatrix& R)
{
  
    assert(M.m() == R.n());  assert(M.n() == R.m());
    SVDRep svd;
    svd.construct(epsilon, M);
    //invert Svd
    double cutoff = epsilon * svd.S()(0);
 
    for(int i=0; i<svd.S().length(); i++) {
        if( svd.S()(i) >= cutoff) {
            svd.S()(i) = 1.0/(svd.S()(i));
        } else {
            assert(0);      //svd.S()(i) = 0.0;
        }
    }
    NumMatrix UT(svd.U().n(),  svd.U().m());
    NumMatrix V( svd.VT().n(), svd.VT().m());
    tran(svd.U(), UT);
    tran(svd.VT(), V);
    for(int i=0; i<V.m(); i++)
        for(int j=0; j<V.n(); j++) {
            V(i,j) = V(i,j) * svd.S()(j);
        }
    //  PetscLogFlops(V.m()*V.n());
    char transa = 'N';
    char transb = 'N';
    double alpha = 1.0;
    double beta = 0.0;
    int m = V.m();
    int n = UT.n();
    int k = V.n();
    DGEMM(&transa, &transb, &m, &n, &k, &alpha,
          V.data(), &m, UT.data(), &k, 
          &beta, R.data(), &m);  
    //  PetscLogFlops( 2*m*n*k );
    return 0;
}

// ---------------------------------------------------------------------- 
int inv(const NumMatrix& M, NumMatrix& R) //Gaussian Elimination
{
    //OR pinv(M, 0.0, R);
    assert(M.m()==M.n() && R.m()==R.n() && M.m()==R.m());
    memcpy(R.data(), M.data(), M.m()*M.n()*sizeof(double));
    int info;
    int m = M.m();
    int* ipiv = new int[m];
    DGETRF(&m, &m, R.data(), &m, ipiv, &info); assert(info==0);
    int lwork = m;
    double* work = new double[lwork];
    DGETRI(&m, R.data(), &m, ipiv, work, &lwork, &info);  assert(info==0);
    delete [] ipiv;
    delete [] work;
    return 0;
}
// -------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------- 
int spev1d(int evflag, int dmflag, int dof, double* data, int m, double e, int i, double u, double* res)
{
  
    int is[4];
    if(dmflag==DMFLAG_PERIOD) {
        for(int k=0; k<4; k++) is[k]=(i+k-1 + m) % m;
    } else {
        assert(i>=1 && i<=m-3);
        for(int k=0; k<4; k++) is[k]=(i+k-1);
    }
    NumMatrix M(dof,m,false,data); //assert(M.n()==n && M.m()==res.m()); //double dof = M.m();  //int cnt = 0;
    //---------------------------
    if(evflag & EVFLAG_VL) {
        double scl = 1.0;
        double us[4]; spcoef(EVFLAG_VL, u, us);
	 
        for(int d=0; d<dof; d++)
            res[d] = 0;
        for(int d=0; d<dof; d++) 
            for(int a=0; a<4; a++) {
                res[d] += us[a] * M(d,is[a]);
            }
        for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
        res+=dof;
    }
    //---------------------------
    if(evflag & EVFLAG_FD) {
        double scl = double(m) / e;
        double us[4]; spcoef(EVFLAG_FD, u, us);
        for(int d=0; d<dof; d++)		res[d] = 0;
        for(int d=0; d<dof; d++)
            for(int a=0; a<4; a++) {
                res[d] += us[a] * M(d,is[a]);
            }
        for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
        res+=dof;
	
    }
    //---------------------------
    if(evflag & EVFLAG_SD) {
        double scl = double(m*m)/(e*e);
        double us[4]; spcoef(EVFLAG_SD, u, us);
        for(int d=0; d<dof; d++)		res[d] = 0;
        for(int d=0; d<dof; d++)
            for(int a=0; a<4; a++) {
                res[d] += us[a] * M(d,is[a]);
            }
        for(int d=0; d<dof; d++)		res[d] *= scl; //scaling
        res+=dof;
    }
    return 0;
}
// ---------------------------------------------------------------------- 
int spev2d(int evflag, int dmflag, int dof, double* data, int* mn, double* ef, int* ij, double* uv, double* res)
{
  
  
    int m = mn[0];  int n = mn[1];
    double e = ef[0];  double f = ef[1];
    int i = ij[0];  int j = ij[1];
    double u = uv[0];  double v = uv[1];
  
    int is[4]; int js[4];
    if(dmflag==DMFLAG_PERIOD) {
        for(int k=0; k<4; k++) is[k]=(i+k-1 + m) % m;
        for(int k=0; k<4; k++) js[k]=(j+k-1 + n) % n;
    } else {
        assert(i>=1 && i<=m-3);
        for(int k=0; k<4; k++)	is[k]=(i+k-1);
        assert(j>=1 && j<=n-3);
        for(int k=0; k<4; k++) js[k]=(j+k-1);
    }
    NumMatrix M(dof,m*n,false,data);
    double scl;
    double us[4], vs[4];
    //---------------------------
    if(evflag & EVFLAG_VL) {
        scl = 1.0;
        spcoef(EVFLAG_VL, u, us);
        spcoef(EVFLAG_VL, v, vs);
        for(int d=0; d<dof; d++)		res[d] = 0;
        for(int a=0; a<4; a++)
            for(int b=0; b<4; b++) {
                double coef = us[a]*vs[b]; 
                for(int d=0; d<dof; d++)
                    res[d] += coef * M(d, is[a]+js[b]*m);
            }
        for(int d=0; d<dof; d++)		res[d] *= scl;
        res+=dof;
    }
    //---------------------------
    if(evflag & EVFLAG_FD) {
        scl = double(m)/e;
        spcoef(EVFLAG_FD, u, us);
        spcoef(EVFLAG_VL, v, vs);
        for(int d=0; d<dof; d++)		res[d] = 0;
        for(int a=0; a<4; a++)
            for(int b=0; b<4; b++) {
                double coef = us[a]*vs[b];
                for(int d=0; d<dof; d++)
                    res[d] += coef * M(d, is[a]+js[b]*m);
            }
        for(int d=0; d<dof; d++)		res[d] *= scl;
        res+=dof;
        //...
        scl = double(n)/f;
        spcoef(EVFLAG_VL, u, us);
        spcoef(EVFLAG_FD, v, vs);
        for(int d=0; d<dof; d++)		res[d] = 0;
        for(int a=0; a<4; a++)
            for(int b=0; b<4; b++) {
                double coef = us[a]*vs[b]; 
                for(int d=0; d<dof; d++)
                    res[d] += coef * M(d, is[a]+js[b]*m);
            }
        for(int d=0; d<dof; d++)		res[d] *= scl;
        res+=dof;
    }
    //---------------------------
    if(evflag & EVFLAG_SD) {
        scl = double(m*m)/(e*e);
        spcoef(EVFLAG_SD, u, us);
        spcoef(EVFLAG_VL, v, vs);
        for(int d=0; d<dof; d++)		res[d] = 0;
        for(int a=0; a<4; a++)
            for(int b=0; b<4; b++) {
                double coef = us[a]*vs[b]; 
                for(int d=0; d<dof; d++)
                    res[d] += coef * M(d, is[a]+js[b]*m);
            }
        for(int d=0; d<dof; d++)		res[d] *= scl;
        res+=dof;
        //...
        scl = double(m*n)/(e*f);
        spcoef(EVFLAG_FD, u, us);
        spcoef(EVFLAG_FD, v, vs);
        for(int d=0; d<dof; d++)		res[d] = 0;
        for(int a=0; a<4; a++)
            for(int b=0; b<4; b++) {
                double coef = us[a]*vs[b]; 
                for(int d=0; d<dof; d++)
                    res[d] += coef * M(d, is[a]+js[b]*m);
            }
        for(int d=0; d<dof; d++)		res[d] *= scl;
        res+=dof;
        //...
        scl = double(n*n)/(f*f);
        spcoef(EVFLAG_VL, u, us);
        spcoef(EVFLAG_SD, v, vs);
        for(int d=0; d<dof; d++)		res[d] = 0;
        for(int a=0; a<4; a++)
            for(int b=0; b<4; b++) {
                double coef = us[a]*vs[b]; 
                for(int d=0; d<dof; d++)
                    res[d] += coef * M(d, is[a]+js[b]*m);
            }
        for(int d=0; d<dof; d++)		res[d] *= scl;
        res+=dof;
    }
    return 0;
}
// ---------------------------------------------------------------------- 
int spcoef(int evflag, double u, double* us)
{
  
    double u1 = u;
    double u2 = u*u;
    double u3 = u*u*u;
    if(       evflag==EVFLAG_VL) {
        us[0] = (  1 - 3*u1 + 3*u2 -   u3)/6.0;
        us[1] = (  4        - 6*u2 + 3*u3)/6.0;
        us[2] = (  1 + 3*u1 + 3*u2 - 3*u3)/6.0;
        us[3] = (                  +   u3)/6.0;
    } else if(evflag==EVFLAG_FD) {
        us[0] = (- 3 + 6*u1 - 3*u2)/6.0;
        us[1] = (    -12*u1 + 9*u2)/6.0;
        us[2] = (  3 + 6*u1 - 9*u2)/6.0;
        us[3] = (             3*u2)/6.0;
    } else if(evflag==EVFLAG_SD) {
        us[0] = (  6 - 6*u1 ) / 6.0;
        us[1] = (-12 +18*u1 ) / 6.0;
        us[2] = (  6 -18*u1 ) / 6.0;
        us[3] = (      6*u1 ) / 6.0;	 //assert(0); //TODO;
    }  
    return 0;
}
