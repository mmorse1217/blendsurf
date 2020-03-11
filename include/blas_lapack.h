#ifndef _BLAS_LAPACK_H_
#define _BLAS_LAPACK_H_

#include "common.hpp"

//blas and lapack workspace for the code


#define DGESVD dgesvd_
#define DGESDD dgesdd_
#define DGETRF dgetrf_
#define DGETRI dgetri_
#define DGEMM dgemm_

//EXTERN_C_BEGIN
extern "C"
{
    extern void DGESVD(char *JOBU, char *JOBVT, int *M, int *N, double *A, int *LDA, 
                       double *S, double *U, int *LDU, double *VT, int *LDVT, double *WORK, int *LWORK, int *INFO);
    extern void DGETRF(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
    extern void DGETRI(int *N, double *A, int *LDA, int *IPIV, double *WORK, int *LWORK, int *INFO);
    extern void DGEMM(char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* ALPHA, double* A,
                      int* LDA, double* B, int* LDB, double* BETA, double* C, int* LDC);
  
}
//EXTERN_C_END

#endif
