/*! \file */
#include "blas_lapack.h"
#include "svdrep.hpp"

int    SVDRep::_wssize = 4194304;
double SVDRep::_wsbuf[4194304];

/* ********************************************************************** */
int SVDRep::construct(double epsilon, const NumMatrix& K)
{
  
  
    int m = K.m();
    int n = K.n();
    int k = std::min(m, n);
  
    NumMatrix tU(m, k);
    NumVector tS(k);
    NumMatrix tVT(k, n);
    int wssize = _wssize;

    //SVD
    int INFO;
    char JOBU  = 'S';
    char JOBVT = 'S';
    assert(_wssize >= max(3*min(m,n)+max(m,n), 5*min(m,n)));
    double* wsbuf = _wsbuf;
    DGESVD(&JOBU, &JOBVT, &m, &n, K.data(), &m, tS.data(), tU.data(), &m, tVT.data(), &k, wsbuf, &wssize, &INFO);
    //assert(INFO==0);

    //cutoff
    double cutoff = epsilon*tS(0);
    int cnt=0;
    while(cnt< k)
        if(std::abs(tS(cnt)) >= cutoff)
            cnt++;
        else
            break;
  
    _matU.resize(m, cnt);
    _matS.resize(cnt);	
    _matVT.resize(cnt,n);
  
    for(int i=0; i<m; i++)
        for(int j=0; j<cnt; j++)
            _matU(i,j) = tU(i,j);
    for(int i=0; i<cnt; i++)
        _matS(i) = tS(i);
    for(int i=0; i<cnt; i++)
        for(int j=0; j<n; j++)
            _matVT(i,j) = tVT(i,j);
  
    return 0;
}
