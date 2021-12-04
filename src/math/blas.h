// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/**
 This contains C front-ends to some functions of BLAS
 see http://www.netlib.org/blas
  
 Functions are renamed depending on what type is REAL:
 
 xcopy() calls scopy() if REAL_IS_FLOAT and dcopy() if REAL_IS_DOUBLE.
*/

#ifndef BLAS_H
#define BLAS_H

//#include "cblas.h"
#define CBLAS_INDEX size_t  /* this may vary between platforms */

#include "real.h"

/// macro will expand to the FORTRAN function name
#if REAL_IS_DOUBLE
#   define BLAS(x) d##x##_
#   define iBLAS(x) id##x##_
#else
#   define BLAS(x) s##x##_
#   define iBLAS(x) is##x##_
#endif


extern "C"
{
// BLAS - Level 1
float  sdot_(int*, const float*, int*, const float*, int*);
double ddot_(int*, const double*, int*, const double*, int*);
double dsdot_(int*, const float*, int*, const float*, int*);
float sdsdot_(int*, const float*, const float*, int*, const float*, int*);

real BLAS(nrm2)(int*, const real*, int*);
real BLAS(asum)(int*, const real*, int*);
real BLAS(sum)(int*, const real*, int*);
CBLAS_INDEX iBLAS(amax)(int*, const real*, int*);
CBLAS_INDEX iBLAS(max)(int*, const real*, int*);
CBLAS_INDEX iBLAS(amin)(int*, const real*, int*);
CBLAS_INDEX iBLAS(min)(int*, const real*, int*);
void BLAS(swap)(int*, real*, int*, real*, int*);
void BLAS(copy)(int*, const real*, int*, real*, int*);
void BLAS(axpy)(int*, real*, const real*, int*, real*, int*);
void BLAS(rotg)(real*, real*, real*, real*);
void BLAS(rotmg)(const real*, const real*, const real*, real*, real*);
void BLAS(rot)(int*, real*, int*, real*, int*, real*, real*);
void BLAS(rotm)(int*, real*, int*, real*, int*, real*);
void BLAS(scal)(int*, real*, real*, int*);

// BLAS - Level 2
void BLAS(gemv)(char*, int*, int*, real*, const real*, int*, const real*, int*, real*, real*, int*);
void BLAS(trmv)(char*, char*, char*, int*, const real*, int*, real*, int*);
void BLAS(trmsv)(char*, char*, char*, int*, const real*, int*, real*, int*);
void BLAS(gbmv)(char*, int*, int*, int*,  int*, real*, const real*, int*, const real*, int*, real*, real*, int*);
void BLAS(tbmv)(char*, char*, char*, int*, int*, const real*, int*, real*, int*);
void BLAS(tbsv)(char*, char*, char*, int*, int*, const real*, int*, real*, int*);

void BLAS(tpsv)(char*, char*, char*, int*, const real*, real*, int*);
void BLAS(tpmv)(char*, char*, char*, int*, const real*, real*, int*);
void BLAS(ger)(int*, int*, real* alpha, const real*, int*, const real*, int*, real*, int*);
void BLAS(symv)(char*, int*, real*, const real*, int*, const real*, int*, real*, real*, int*);
void BLAS(sbmv)(char*, int*, int*, real*, const real*, int*, const real*, int*, real*, real*, int*);
void BLAS(spmv)(char*, int*, real*, const real*, const real*, int*, real*, real*, int*);
void BLAS(syr)(char*, int*, real*, const real*, int*, real*, int*);
void BLAS(syr2)(char*, int*, real*, const real*, int*, const real*, int*, real*, int*);
void BLAS(spr)(char*, int*, real*, const real*, int*, real*);
void BLAS(spr2)(char*, int*, real*, const real*, int*, const real*, int*, real*);

// BLAS - Level 3
void BLAS(gemm)(char*, char*, int*, int*, int*, real*, const real*, int*, const real*, int*, real*, real*, int*);
void BLAS(symm)(char*, char*, int*, int*, real*, const real*, int*, const real*, int*, real*, real*, int*);
void BLAS(syrk)(char*, char*, int*, int*, real*, const real*, int*, real*, real*, int*);
void BLAS(xsyr2k)(char*, char*, int*, int*, real*, const real*, int*, const real*, int*, real*, real*, int*);
void BLAS(xtrmm)(char*, char*, char*, int*, int*, real*, const real*, int*, real*, int*);
void BLAS(trsm)(char*, char*, char*, char*, int*, int*, real*, const real*, int*, real*, int*);
}


namespace blas
{
#pragma mark - Level 1

/**
 This follows William Kahan's summation to calculate the dot product
 with minimal numerical error.
 https://en.wikipedia.org/wiki/Kahan_summation_algorithm
 */
template < typename FLOAT >
inline FLOAT dot_kahan(int N, const real* X, const real* Y)
{
    FLOAT sum = 0.0;
    FLOAT c = 0.0;
    
    for ( int i = 0; i < N; ++i )
    {
        FLOAT y = X[i] * Y[i] - c;
        FLOAT t = sum + y;
        c = ( t - sum ) - y;
        sum = t;
    }
    return sum;
}


inline double dot(int N, const real* X, const real* Y)
{
    int ONE = 1;
    double res;
#if REAL_IS_DOUBLE
    res = ddot_(&N, X, &ONE, Y, &ONE);
#else
    res = dsdot_(&N, X, &ONE, Y, &ONE);
#endif
#if 0
    real k = dot_kahan<double>(N, X, Y);
    printf("   sum-kahan = %+8.3e (%+8.3e)\n", res-k, (res-k)/k);
#endif
    return res;
}

/**
 We always use double precision to accumulate the dot product of two vectors:
 */
inline double xdot(int N, const real* X, int incX, const real* Y, int incY)
{
#if REAL_IS_DOUBLE
    return ddot_(&N, X, &incX, Y, &incY);
#else
    return dsdot_(&N, X, &incX, Y, &incY);
#endif
}

inline real ddot(int N, const double* X, int incX, const double* Y, int incY)
{
    return ddot_(&N, X, &incX, Y, &incY);
}

inline double dsdot(int N, const float* X, int incX, const float* Y, int incY)
{
    return dsdot_(&N, X, &incX, Y, &incY);
}

inline float sdsdot(int N, float SB, const float* X, int incX, const float* Y, int incY)
{
    return sdsdot_(&N, &SB, X, &incX, Y, &incY);
}

// use 'blas::nrm2' defined above if applicable
inline real xnrm2(int N, const real*X, int incX)
{
    return BLAS(nrm2)(&N, X, &incX);
}

inline real xasum(int N, const real*X, int incX)
{
    return BLAS(asum)(&N, X, &incX);
}

inline real xsum(int N, const real*X, int incX)
{
    return BLAS(sum)(&N, X, &incX);
}

inline CBLAS_INDEX ixamax(int N, const real*X, int incX)
{
    return iBLAS(amax)(&N, X, &incX);
}

inline CBLAS_INDEX ixmax(int N, const real*X, int incX)
{
    return iBLAS(max)(&N, X, &incX);
}

inline CBLAS_INDEX ixamin(int N, const real*X, int incX)
{
    return iBLAS(amin)(&N, X, &incX);
}

inline CBLAS_INDEX ixmin(int N, const real*X, int incX)
{
    return iBLAS(min)(&N, X, &incX);
}

inline void xswap(int N, real*X, int incX, real*Y, int incY)
{
    BLAS(swap)(&N, X, &incX, Y, &incY);
}

inline void xcopy(int N, const real*X, int incX, real*Y, int incY)
{
    BLAS(copy)(&N, X, &incX, Y, &incY);
}

inline void copy(int N, const real* X, real* Y)
{
    //copy_real(N, X, Y);
    blas::xcopy(N, X, 1, Y, 1);
}

inline void xaxpy(int N, real alpha, const real*X, int incX, real*Y, int incY)
{
    BLAS(axpy)(&N, &alpha, X, &incX, Y, &incY);
}

inline void xrotg(real*a, real*b, real*c, real*s)
{
    BLAS(rotg)(a, b, c, s);
}

inline void xrotmg(const real*d1, const real*d2, const real*b1, real b2, real*P)
{
    BLAS(rotmg)(d1, d2, b1, &b2, P);
}

inline void xrot( int N, real*X, int incX, real*Y, int incY, real c, real s)
{
    BLAS(rot)(&N, X, &incX, Y, &incY, &c, &s);
}

inline void xrotm( int N, real*X, int incX, real*Y, int incY, real*P)
{
    BLAS(rotm)(&N, X, &incX, Y, &incY, P);
}

inline void xscal(int N, real alpha, real*X, int incX)
{
    BLAS(scal)(&N, &alpha, X, &incX);
}


#pragma mark - Level 2


inline void xgemv(char TransA, int M, int N, real alpha, const real*A, int lda, const real*X, int incX, real beta, real*Y, int incY)
{
    BLAS(gemv)(&TransA, &M, &N, &alpha, A, &lda, X, &incX, &beta, Y, &incY);
}

inline void xtrmv(char Uplo, char TransA, char Diag, int N, const real*A, int lda, real*X, int incX)
{
    BLAS(trmv)(&Uplo, &TransA, &Diag, &N, A, &lda, X, &incX);
}

inline void xtrsv(char Uplo, char TransA, char Diag, int N, const real*A, int lda, real*X, int incX);

inline void xgbmv(char TransA, int M, int N, int Kl,  int Ku, real alpha, const real*A, int lda, const real*X, int incX, real beta, real*Y, int incY);

inline void xtbmv(char Uplo, char Trans, char Diag, int N, int K, const real*A, int lda, real*X, int incX)
{
    BLAS(tbmv)(&Uplo, &Trans, &Diag, &N, &K, A, &lda, X, &incX);
}

inline void xtbsv(char Uplo, char Trans, char Diag, int N, int K, const real* A, int lda, real*X, int incX)
{
    BLAS(tbsv)(&Uplo, &Trans, &Diag, &N, &K, A, &lda, X, &incX);
}

inline void xtpsv(char Uplo, char TransA, char Diag, int N, const real*A, real*X, int incX);

inline void xtpmv(char Uplo, char TransA, char Diag, int N, const real*A, real*X, int incX);

inline void xger(int M, int N, real alpha, const real*X, int incX, const real*Y, int incY, real*A, int lda)
{
    BLAS(ger)(&M, &N, &alpha, X, &incX, Y, &incY, A, &lda);
}

inline void xsymv(char Uplo, int N, real alpha,  const real*A, int lda, const real*X, int incX, real beta, real*Y, int incY)
{
    BLAS(symv)(&Uplo, &N, &alpha, A, &lda, X, &incX, &beta, Y, &incY);
}

inline void xsbmv(char Uplo, int N, int K, real alpha, const real*A, int lda, const real*X, int incX, real beta, real*Y, int incY)
{
    BLAS(sbmv)(&Uplo, &N, &K, &alpha, A, &lda, X, &incX, &beta, Y, &incY);
}

inline void xspmv(char Uplo, int N, real alpha, const real*A, const real*X, int incX, real beta, real*Y, int incY)
{
    BLAS(spmv)(&Uplo, &N, &alpha, A, X, &incX, &beta, Y, &incY);
}

inline void xsyr(char Uplo, int N, real alpha, const real*X, int incX, real*A, int lda)
{
    BLAS(syr)(&Uplo, &N, &alpha, X, &incX, A, &lda);
}

inline void xsyr2(char Uplo, int N, real alpha, const real*X, int incX, const real*Y, int incY, real* A, int lda)
{
    BLAS(syr2)(&Uplo, &N, &alpha, X, &incX, Y, &incY, A, &lda);
}

inline void xspr(char Uplo, int N, real alpha, const real*X, int incX, real*A)
{
    BLAS(spr)(&Uplo, &N, &alpha, X, &incX, A);
}

inline void xspr2(char Uplo, int N, real alpha, const real*X, int incX, const real*Y, int incY, real*A)
{
    BLAS(spr2)(&Uplo, &N, &alpha, X, &incX, Y, &incY, A);
}


#pragma mark - Level 3


inline void xgemm(char TransA, char TransB, int M, int N, int K, real alpha, const real*A, int lda, const real*B, int ldb, real beta, real*C, int ldc)
{
    BLAS(gemm)(&TransA, &TransB, &M, &N, &K, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

inline void xsymm(char Side, char Uplo, int M, int N, real alpha, const real*A, int lda, const real*B, int ldb, real beta, real*C, int ldc)
{
    BLAS(symm)(&Side, &Uplo, &M, &N, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

inline void xsyrk(char Uplo, char Trans, int N, int K, real alpha, const real*A, int lda, real beta, real*C, int ldc)
{
    BLAS(syrk)(&Uplo, &Trans, &N, &K, &alpha, A, &lda, &beta, C, &ldc);
}

inline void xsyr2k(char Uplo, char Trans, int N, int K, real alpha, const real*A, int lda, const real*B, int ldb, real beta, real*C, int ldc);

inline void xtrmm(char Uplo, char TransA, char Diag, int M, int N, real alpha, const real*A, int lda, real*B, int ldb);

inline void xtrsm(char side, char uplo, char transA, char diag, int M, int N, real alpha, const real*A, int lda, real*B, int ldb)
{
    BLAS(trsm)(&side, &uplo, &transA, &diag, &M, &N, &alpha, A, &lda, B, &ldb);
}

}
#undef CBLAS_INDEX

#endif
