// Cytosim was created by Francois Nedelec.
// Copyright 2020 Sainsbury Laboratory, Cambridge University

#ifndef XTBSV_H
#define XTBSV_H

#if defined(__AVX__)
#  include "simd.h"
#  include "simd_float.h"
#elif defined(__SSE3__)
#  include "simd_float.h"
#endif

/**
 Enables an optimization that speeds-up some routines on Broadwell architecture
 by expanding multiplications to reduce the dependency chain of the calculation
 ---> compare the two implementations by running 'test_xtbsv'
 */
#ifndef DEVELOP_XTBSV
#   define DEVELOP_XTBSV 1
#endif


/**
 This is a C-translation of the BLAS reference implementation of DTBSV
 FJN 28.04.2020
 
 DTBSV  solves one of the systems of equations
 
     A*x = b,   or   A**T*x = b,
 
 where b and x are n element vectors and A is an n by n unit, or non-unit,
 upper or lower triangular band matrix, with ( k + 1 ) diagonals.
 */
template < char diag >
void blas_xtbsvUN(const int N, const int KD, const real* A, const int lda, real* X, const int incX)
{
    assert_true( lda > KD );
    assert_true( incX != 0 );
    int kx = 0;
    if ( incX > 0 )
        kx = (N-1) * incX;
    int jx = kx;
    for (int j = N-1; j >= 0; --j)
    {
        kx -= incX;
        if (X[jx] != 0.)
        {
            int ix = kx;
            const real * pA = A + KD + j * lda;
            if ( diag == 'N' ) X[jx] /= pA[0];
            else if ( diag == 'I' ) X[jx] *= pA[0];
            real tmp = X[jx];
            const int inf = std::max(0, j-KD);
            for (int i = j - 1; i >= inf; --i)
            {
                X[ix] -= tmp * pA[i-j];
                ix -= incX;
            }
        }
        jx -= incX;
    }
}

template < char diag >
void blas_xtbsvUN(const int N, const int KD, const real* A, const int lda, real* X)
{
    assert_true( lda > KD );
    for (int j = N-1; j >= 0; --j)
    {
        if (X[j] != 0.)
        {
            const real * pA = A + KD + j * lda;
            if ( diag == 'N' ) X[j] /= pA[0];
            else if ( diag == 'I' ) X[j] *= pA[0];
            real tmp = X[j];
            const int inf = std::max(0, j-KD);
            for (int i = j - 1; i >= inf; --i)
                X[i] -= tmp * pA[i-j];
        }
    }
}


template < char diag >
void blas_xtbsvLN(const int N, const int KD, const real* A, const int lda, real* X, const int incX)
{
    assert_true( lda > KD );
    assert_true( incX != 0 );
    int kx = 0;
    if ( incX <= 0 )
        kx = - (N-1) * incX;
    int jx = kx;
    for (int j = 0; j < N; ++j)
    {
        kx += incX;
        if (X[jx] != 0.)
        {
            int ix = kx;
            real tmp = X[jx];
            const real * pA = A + j * lda;
            if ( diag == 'N' )
            { tmp /= pA[0]; X[jx] = tmp; }
            else if ( diag == 'I' )
            { tmp *= pA[0]; X[jx] = tmp; }
            else if ( diag == 'C' )
            { X[jx] = tmp * pA[0]; }
            const int sup = std::min(N-1, j+KD);
            for (int i = j + 1; i <= sup; ++i)
            {
                X[ix] -= tmp * pA[i-j];
                ix += incX;
            }
        }
        jx += incX;
    }
}

template < char diag >
void blas_xtbsvLN(const int N, const int KD, const real* A, const int lda, real* X)
{
    assert_true( lda > KD );
    for (int j = 0; j < N; ++j)
    {
        if ( X[j] != 0 )
        {
            const real * pA = A + j * lda;
            real tmp = X[j];
            if ( diag == 'N' )
            { tmp /= pA[0]; X[j] = tmp; }
            else if ( diag == 'I' )
            { tmp *= pA[0]; X[j] = tmp; }
            else if ( diag == 'C' )
            { X[j] = tmp * pA[0]; }
            const int sup = std::min(N-1, j+KD);
            for (int i = j + 1; i <= sup; ++i)
                X[i] -= tmp * pA[i-j];
        }
    }
}


template < char diag >
void blas_xtbsvUT(const int N, const int KD, const real* A, const int lda, real* X, const int incX)
{
    assert_true( lda > KD );
    assert_true( incX != 0 );
    int kx = 0;
    if ( incX <= 0 )
        kx = - (N-1) * incX;
    int jx = kx;
    for (int j = 0; j < N; ++j)
    {
        real tmp = X[jx];
        int ix = kx;
        const real * pA = A + KD + j * lda;
        for (int i = std::max(0, j-KD); i < j; ++i)
        {
            tmp -= pA[i-j] * X[ix];
            ix += incX;
        }
        if ( diag == 'N' ) tmp /= pA[0];
        else if ( diag == 'I' ) tmp *= pA[0];
        X[jx] = tmp;
        jx += incX;
        if (j >= KD)
            kx += incX;
    }
}

template < char diag >
void blas_xtbsvUT(const int N, const int KD, const real* A, const int lda, real* X)
{
    assert_true( lda > KD );
    for (int j = 0; j < N; ++j)
    {
        real tmp = X[j];
        const real * pA = A + KD + j * lda;
        for (int i = std::max(0, j-KD); i < j; ++i)
            tmp -= pA[i-j] * X[i];
        if ( diag == 'N' ) tmp /= pA[0];
        else if ( diag == 'I' ) tmp *= pA[0];
        X[j] = tmp;
    }
}


template < char diag >
void blas_xtbsvLT(const int N, const int KD, const real* A, const int lda, real* X, const int incX)
{
    assert_true( lda > KD );
    assert_true( incX != 0 );
    int kx = 0;
    if ( incX > 0 )
        kx = (N-1) * incX;
    int jx = kx;
    for (int j = N-1; j >= 0; --j)
    {
        real tmp = X[jx];
        int ix = kx;
        const real * pA = A + j * lda;
        if ( diag == 'C' ) tmp *= pA[0];
        const int sup = std::min(N-1, j+KD);
        for (int i = sup; i > j; --i)
        {
            tmp -= pA[i-j] * X[ix];
            ix -= incX;
        }
        if ( diag == 'N' ) tmp /= pA[0];
        else if ( diag == 'I' ) tmp *= pA[0];
        X[jx] = tmp;
        jx -= incX;
        if ( j < N-KD )
            kx -= incX;
    }
}

template < char diag >
void blas_xtbsvLT(const int N, const int KD, const real* A, const int lda, real* X)
{
    assert_true( lda > KD );
    for (int j = N-1; j >= 0; --j)
    {
        real tmp = X[j];
        const real * pA = A + j * lda;
        if ( diag == 'C' ) tmp *= pA[0];
        const int sup = std::min(N-1, j+KD);
        for (int i = sup; i > j; --i)
            tmp -= pA[i-j] * X[i];
        if ( diag == 'N' ) tmp /= pA[0];
        else if ( diag == 'I' ) tmp *= pA[0];
        X[j] = tmp;
    }
}

    
void blas_xtbsv(char Uplo, char Trans, char Diag, const int N, const int KD, const real* A, const int lda, real* X, const int incX)
{
    if ( Uplo == 'U' )
    {
        if ( Trans == 'N' ) {
            if ( Diag == 'N' )
                blas_xtbsvUN<'N'>(N, KD, A, lda, X, incX);
            else
                blas_xtbsvUN<'U'>(N, KD, A, lda, X, incX);
        } else {
            if ( Diag == 'N' )
                blas_xtbsvUT<'N'>(N, KD, A, lda, X, incX);
            else
                blas_xtbsvUT<'U'>(N, KD, A, lda, X, incX);
        }
    }
    else if ( Uplo == 'L' )
    {
        if ( Trans == 'N' ) {
            if ( Diag == 'N' )
                blas_xtbsvLN<'N'>(N, KD, A, lda, X, incX);
            else
                blas_xtbsvLN<'U'>(N, KD, A, lda, X, incX);
        } else {
            if ( Diag == 'N' )
                blas_xtbsvLT<'N'>(N, KD, A, lda, X, incX);
            else
                blas_xtbsvLT<'U'>(N, KD, A, lda, X, incX);
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - Alsatian factorization xPBTF2

/*
  SUBROUTINE DSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
 
    IF (LSAME(UPLO,'L')) THEN
    IF (INCX.EQ.1) THEN

       DO 60 J = 1,N
            IF (X(J).NE.ZERO) THEN
                TEMP = ALPHA*X(J)
                DO 50 I = J,N
                    A(I,J) = A(I,J) + X(I)*TEMP
50                 CONTINUE
            END IF
60         CONTINUE
    END IF
    END IF
 */
void blas_xsyrL(int N, real ALPHA, const real* X, real* A, int LDA)
{
    for ( int J = 0; J < N; ++J )
    {
        real tmp = ALPHA * X[J];
        for ( int I = J; I < N; ++I )
            A[I] += X[I] * tmp;
        A += LDA;
    }
}


/**
 This is equivalent to calling the standard lapack::pbtf2()
 and then *** inverting *** the diagonal terms
 
 SUBROUTINE DPBTF2( UPLO='L', N, KD, AB, LDAB, INFO )
*/
void alsatian_xpbtf2L(const int N, const int KD, real* AB, const int LDAB, int* INFO)
{
    int KLD = std::max(1, LDAB-1);
    //Compute the Cholesky factorization A = L*L**T.
    for ( int J = 0; J < N; ++J )
    {
        //Compute L(J,J) and test for non-positive-definiteness.
        real dia = AB[0];
        if ( dia <= 0 )
        {
            *INFO = J+1;
            return;
        }
        dia = real(1) / dia;
        AB[0] = std::sqrt(dia);
        int KN = std::min(KD, N-1-J);  // N-1-J < KD iff J >= N-KD
        // update the trailing submatrix within the band:
        real* A = AB + LDAB;
        for ( int K = 0; K < KN; ++K )
        {
            // update column J+K of AB
            real tmp = AB[K+1] * dia;
            for ( int I = K; I < KN; ++I )
                A[I] -= AB[I+1] * tmp;
            A += KLD;
        }
        // scale off diagonal terms in column:
        for ( int K = 1; K <= KN; ++K )
            AB[K] *= dia;
        AB += LDAB;
    }
    *INFO = 0;
}

//------------------------------------------------------------------------------
#pragma mark - Alsatian DTBSV with KD as argument

#if 1
//The original LAPACK implementation loads and stores X[] at every operation
void alsatian_xtbsvLNN_ORI(const int N, const int KD, const real* A, const int lda, real* X)
{
    for ( int j = 0; j < N; ++j )
    {
        const real tmp = X[j];
        X[j] = tmp * A[0];
        const int sup = std::min(N-1, j+KD);  // sup=N-1 if N-1<j+KD : j >= N-KD
        for ( int i = j + 1; i <= sup; ++i )
            X[i] -= tmp * A[i];
        A += lda;
    }
}

void alsatian_xtbsvLTN_ORI(const int N, const int KD, const real* A, const int lda, real* X)
{
    A += (N-1)*(lda-1);
    for ( int j = N-1; j >= 0; --j )
    {
        real tmp = X[j] * A[0];
        const int sup = std::min(N-1, j+KD);  // sup=N-1 if N-1<j+KD : j >= N-KD
        for ( int i = sup; i > j; --i )
            tmp -= A[i] * X[i];
        X[j] = tmp;
        A -= lda;
    }
}
#endif

void alsatian_xtbsvLNN(const int N, const int KD, const real* A, const int lda, real* X)
{
    const int nkd = std::max(0, N-KD);
    // general case:
    for ( int j = 0; j < nkd; ++j )
    {
        const real tmp = X[j];
        X[j] = tmp * A[0];
        for ( int i = 1; i <= KD; ++i )
            X[i+j] -= tmp * A[i];
        A += lda;
    }
    // process truncated cases:
    for ( int j = nkd; j < N; ++j )
    {
        const real tmp = X[j];
        X[j] = tmp * A[0];
        for ( int i = j + 1; i < N; ++i )
            X[i] -= tmp * A[i-j];
        A += lda;
    }
}

void alsatian_xtbsvLTN(const int N, const int KD, const real* A, const int lda, real* X)
{
    A += (N-1)*lda;
    const int nkd = std::max(0, N-KD);
    // process KD truncated cases:
    for ( int j = N-1; j >= nkd; --j )
    {
        real tmp = A[0] * X[j];
        for ( int i = N-1; i > j; --i )
            tmp -= A[i-j] * X[i];
        X[j] = tmp;
        A -= lda;
    }
    // general case, downward:
    for ( int j = nkd-1; j >= 0; --j )
    {
        real tmp = A[0] * X[j];
        for ( int i = KD; i > 0; --i )
            tmp -= A[i] * X[i+j];
        X[j] = tmp;
        A -= lda;
    }
}


//------------------------------------------------------------------------------
#pragma mark - Alsatian, KD as template argument, using C-array as buffer

template < int KD >
void alsatian_xtbsvLNNK(const int N, const real* A, const int lda, real* X)
{
    assert_true( lda > KD );
    const int nkd = std::max(0, N-KD);
    real buf[KD] = { 0 };
    if ( 0 < N )
    {
        // 'general case' with j = 0, initializing buf[]
        real t = X[0];
        X[0] = t * A[0];
        for ( int i = 1; i <= std::min(N,KD); ++i )
            buf[i-1] = X[i] - t * A[i];
        A += lda;
        X += 1;
    }
    // general case:
    for ( int j = 1; j < nkd; ++j )
    {
        /*
        const real t = X[j];
        X[j] = t * A[0];
        for ( int i = 1; i <= KD; ++i )
            X[i+j] -= t * A[i];
         */
        real t = buf[0];
        X[0] = t * A[0];
        for ( int i = 1; i < KD; ++i )
            buf[i-1] = buf[i] - t * A[i]; // buf[i] contains X[i+j]
        buf[KD-1] = X[KD] - t * A[KD];
        A += lda;
        X += 1;
    }
    // process KD truncated cases:
    for ( int j = nkd; j < N; ++j )
    {
        /*
        const real t = X[j] * A[0];
        X[j] = t;
        for ( int i = j + 1; i < N; ++i )
            X[i] -= t * A[i-j];
         */
        real t = buf[0];
        X[0] = t * A[0];
        for ( int i = 1; i < std::min(KD, N-j); ++i )
            buf[i-1] = buf[i] - t * A[i]; // buf[i] contains X[i+j]
        A += lda;
        X += 1;
    }
}


template < int KD >
void alsatian_xtbsvLTNK(const int N, const real* A, const int lda, real* X)
{
    assert_true( lda > KD );
    A += ( N - 1 ) * lda;
    X += N - 1;
    const int nkd = std::max(0, N-KD);
    real buf[KD];
    // process KD truncated cases:
    for ( int j = N-1; j >= nkd; --j )
    {
        real tmp = A[0] * X[0];
/*
        for ( int i = j + 1; i < N; ++i )
            tmp -= A[i-j] * X[i];
*/
        for ( int i = N-1; i > j; --i )
            tmp -= A[i-j] * buf[i-nkd]; // buf[] is X[i];
        buf[j-nkd] = tmp;
        X[0] = buf[j-nkd];
        A -= lda;
        X -= 1;
    }
    // at this stage, buf[1] = X[N-KD] and buf[KD] = X[N-1]:
    // general case, downward!
    for ( int j = nkd-1; j >= 0; --j )
    {
        /*
         real tmp = X[j];
         for ( int i = 1; i <= KD; ++i )
            tmp -= A[i] * X[i+j];
         */
        real tmp = A[0] * X[0] - A[KD] * buf[KD-1];
        for ( int i = KD-1; i > 0; --i )
        {
            tmp -= A[i] * buf[i-1];
            buf[i] = buf[i-1];
        }
        buf[0] = tmp;
        X[0] = buf[0];
        A -= lda;
        X -= 1;
    }
}


//------------------------------------------------------------------------------
#pragma mark - Specialized versions for KD==6


/// Specialized version for KD==6
void alsatian_xtbsvLNN6K(const int N, const real* A, const int lda, real* X)
{
    const real * end = X + N - 6;
    //constexpr int KD = 6;
    assert_true( lda > 6 );
    assert_true( N >= 6 );
    real tt = X[0];
    X[0] = tt * A[0];
    real x0 = X[1] - tt * A[1];
    real x1 = X[2] - tt * A[2];
    real x2 = X[3] - tt * A[3];
    real x3 = X[4] - tt * A[4];
    real x4 = X[5] - tt * A[5];
    real x5 = X[6] - tt * A[6];
    A += lda;
    ++X;
    for ( ; X < end; ++X )
    {
        tt = x0;
        X[0] = x0 * A[0];
        x0 = x1 - tt * A[1];
        x1 = x2 - tt * A[2];
        x2 = x3 - tt * A[3];
        x3 = x4 - tt * A[4];
        x4 = x5 - tt * A[5];
        x5 = X[6] - tt * A[6];
        A += lda;
    }
    for ( ; X < end+6; ++X )
    {
        tt = x0;
        X[0] = x0 * A[0];
        x0 = x1 - tt * A[1];
        x1 = x2 - tt * A[2];
        x2 = x3 - tt * A[3];
        x3 = x4 - tt * A[4];
        x4 = x5 - tt * A[5];
        x5 = 0;
        A += lda;
    }
}


/// Specialized version for KD==6
void alsatian_xtbsvLTN6K(const int N, const real* A, const int lda, real* X)
{
    const real * end = X;
    //constexpr int KD = 6;
    assert_true( lda > 6 );
    assert_true( N >= 6 );
    A += ( N - 1 ) * lda;
    X += N - 1;
    real x0 = 0, x1 = 0, x2 = 0, x3 = 0, x4 = 0, x5 = 0;
    for ( ; X >= end; --X )
    {
        real tt = x5 * A[6];
        tt += x4 * A[5];
        x5 = x4;
        x4 = x3;
        tt += x3 * A[4];
        tt += x2 * A[3];
        x3 = x2;
        x2 = x1;
        tt += x1 * A[2];
        tt += x0 * A[1];
        x1 = x0;
        x0 = A[0] * X[0] - tt;
        X[0] = x0;
        A -= lda;
    }
}


//------------------------------------------------------------------------------
#pragma mark - Optimized SSE versions for KD==6

#if USE_SIMD
/**
 Optmized version for KD == 6
 Beware: this works assuming that N >= KD, and it will in particular write
 to X[i] for i = 0 ... 6 for any value of N.
*/
void alsatian_xtbsvLNN6K_SSE(const int N, const double* A, const int lda, double* X)
{
    constexpr int KD = 6;
    assert_true( lda > KD );
    assert_true( N >= KD );
    const double * end = X + N - KD;
    vec2 x0, x2, x4, tt;
    // first iteration with j = 0
    tt = loaddup2(X);
    x0 = fnmadd2(tt, loadu2(A+1), loadu2(X+1));
    x2 = fnmadd2(tt, loadu2(A+3), loadu2(X+3));
    x4 = fnmadd2(tt, loadu2(A+5), loadu2(X+5)); // may load garbage if N <= KD
    store1(X, mul1(tt, load1(A)));
    A += lda;
    X += 1;
    // general case:
    for ( ; X < end; ++X ) //for ( int j = 1; j < nkd; ++j )
    {
        tt = unpacklo2(x0, x0);
        store1(X, mul1(load1(A), x0));
        x0 = catshift(x0, x2);
        x2 = catshift(x2, x4);
        x4 = catshift(x4, load1(X+KD));
        x0 = fnmadd2(tt, loadu2(A+1), x0);
        x2 = fnmadd2(tt, loadu2(A+3), x2);
        x4 = fnmadd2(tt, loadu2(A+5), x4);
        A += lda;
    }
    /*
     The ending sequence avoids loading elements of X beyond X+N, and elements
     of A in the last column that should normally be equal to zero, but otherwise
     it performs the same calculation than the normal iteration six times.
     */
    // process truncated case: j = N - 6
    tt = unpacklo2(x0, x0);
    vec2 yy = mul1(load1(A), x0);
    x0 = catshift(x0, x2);
    x2 = catshift(x2, x4);
    x4 = unpackhi2(x4, setzero2());
    x0 = fnmadd2(tt, loadu2(A+1), x0);
    x2 = fnmadd2(tt, loadu2(A+3), x2);
    x4 = fnmadd1(tt, load1(A+5), x4);
    A += lda;
    // process truncated case: j = N - 5
    tt = unpacklo2(x0, x0);
    storeu2(X, unpacklo2(yy, mul1(load1(A), x0)));
    x0 = catshift(x0, x2);
    x2 = catshift(x2, x4);
    x0 = fnmadd2(tt, loadu2(A+1), x0);
    x2 = fnmadd2(tt, loadu2(A+3), x2);
    A += lda;
    // process truncated case: j = N - 4
    tt = unpacklo2(x0, x0);
    yy = mul1(load1(A), x0);
    x0 = catshift(x0, x2);
    x2 = catshift(x2, setzero2());
    x0 = fnmadd2(tt, loadu2(A+1), x0);
    x2 = fnmadd1(tt, load1(A+3), x2);
    A += lda;
    // process truncated case: j = N - 3
    tt = unpacklo2(x0, x0);
    storeu2(X+2, unpacklo2(yy, mul1(load1(A), x0)));
    x0 = catshift(x0, x2);
    x0 = fnmadd2(tt, loadu2(A+1), x0);
    A += lda;
    // process last two scalars:
    if ( N <= KD ) {
        yy = mul1(load1(A), x0);
        store1(X+4, yy);
        return;
    }
    yy = mul1(load1(A), x0);
    x2 = unpackhi2(x0, x0);
    x2 = fnmadd1(x0, load1(A+1), x2);
    x2 = mul1(load1(A+lda), x2);
    storeu2(X+4, unpacklo2(yy, x2));
}


/**
 Optmized version for KD == 6
 Beware: this works assuming that N >= KD, and it will in particular write
 to X[i] for i = 0 ... 6 for any value of N.
*/
void alsatian_xtbsvLTN6K_SSE(const int N, const double* A, const int lda, double* X)
{
    const double * end = X;
    //constexpr int KD = 6;
    assert_true( lda > 6 );
    assert_true( N >= 6 );
    A += ( N - 1 ) * lda;
    X += N - 1;
    vec2 x0 = setzero2();
    vec2 x1 = setzero2();
    vec2 x2 = setzero2();
    vec2 x3 = setzero2();
    vec2 x4 = setzero2();
    vec2 x5 = setzero2();
#if 1
    /*
     The starting sequence avoids loading elements of X beyond X+N, and elements
     of A in the last column that should normally be equal to zero, but otherwise
     it performs the same calculation than the normal iteration six times.
     */
    // truncated round at j = N-1
    {
        x0 = mul1(load1(A), load1(X));
        A -= lda;
    }
    // truncated round at j = N-2
    {
        vec2 tt = mul1(load1(A), load1(X-1));
        x1 = x0;
        x0 = fnmadd1(x0, load1(A+1), tt);
        store2(X-1, unpacklo2(x0, x1));
        A -= lda;
    }
    // truncated round at j = N-3
    {
        vec2 tt = mul1(load1(A), load1(X-2));
        vec2 ab = loadu2(A+1);
        x2 = x1;
        tt = fnmadd1(x1, unpackhi2(ab, ab), tt); x1 = x0;
        x0 = fnmadd1(x0, unpacklo2(ab, ab), tt);
        A -= lda;
    }
    // truncated round at j = N-4
    {
        vec2 tt = mul1(load1(A), load1(X-3));
        x3 = x2;
        tt = fnmadd1(x2, load1(A+3), tt); x2 = x1;
        vec2 ab = loadu2(A+1);
        tt = fnmadd1(x1, unpackhi2(ab, ab), tt); x1 = x0;
        x0 = fnmadd1(x0, unpacklo2(ab, ab), tt);
        store2(X-3, unpacklo2(x0, x1));
        A -= lda;
    }
    // truncated round at j = N-5
    {
        vec2 tt = mul1(load1(A), load1(X-4));
        x4 = x3;
        vec2 ab = loadu2(A+3);
        tt = fnmadd1(x3, unpackhi2(ab, ab), tt); x3 = x2;
        tt = fnmadd1(x2, unpacklo2(ab, ab), tt); x2 = x1;
        ab = loadu2(A+1);
        tt = fnmadd1(x1, unpackhi2(ab, ab), tt); x1 = x0;
        x0 = fnmadd1(x0, unpacklo2(ab, ab), tt);
        A -= lda;
    }
    // truncated round at j = N-6
    {
        vec2 tt = mul1(load1(A), load1(X-5));
        x5 = x4;
        tt = fnmadd1(x4, load1(A+5), tt); x4 = x3;
        vec2 ab = loadu2(A+3);
        tt = fnmadd1(x3, unpackhi2(ab, ab), tt); x3 = x2;
        tt = fnmadd1(x2, unpacklo2(ab, ab), tt); x2 = x1;
        ab = loadu2(A+1);
        tt = fnmadd1(x1, unpackhi2(ab, ab), tt); x1 = x0;
        x0 = fnmadd1(x0, unpacklo2(ab, ab), tt);
        store2(X-5, unpacklo2(x0, x1));
        A -= lda;
    }
    // update pointer
    X -= 6;
#endif
    for ( ; X >= end; --X )
    {
        vec2 tt = mul1(load1(A), load1(X));
        vec2 ab = loadu2(A+5);
        tt = fnmadd1(x5, unpackhi2(ab, ab), tt); x5 = x4;
        tt = fnmadd1(x4, unpacklo2(ab, ab), tt); x4 = x3;
        ab = loadu2(A+3);
        tt = fnmadd1(x3, unpackhi2(ab, ab), tt); x3 = x2;
        tt = fnmadd1(x2, unpacklo2(ab, ab), tt); x2 = x1;
        ab = loadu2(A+1);
        tt = fnmadd1(x1, unpackhi2(ab, ab), tt); x1 = x0;
        x0 = fnmadd1(x0, unpacklo2(ab, ab), tt);
        store1(X, x0);
        A -= lda;
    }
}

void alsatian_xtbsvLTN6K_SSE_ALT(const int N, const double* A, const int lda, double* X)
{
    const double * end = X;
    //constexpr int KD = 6;
    assert_true( lda > 6 );
    assert_true( N >= 6 );
    A += ( N - 1 ) * lda;
    X += N - 1;
    vec2 tt;
    vec2 x0 = setzero2();
    vec2 x2 = setzero2();
    vec2 x4 = setzero2();
    // general case, downward!
    for ( ; X >= end; --X ) //for ( int j = nkd-1; j >= 0; --j )
    {
        tt = mul1(load1(A), load1Z(X));
        tt = fnmadd2(x4, loadu2(A+5), tt);
        tt = fnmadd2(x2, loadu2(A+3), tt);
        x4 = catshift(x2, x4);
        x2 = catshift(x0, x2);
        tt = fnmadd2(x0, loadu2(A+1), tt);
        tt = add1(tt, swap2(tt));
        x0 = unpacklo2(tt, x0);
        store1(X, tt);
        A -= lda;
    }
}
#endif

//------------------------------------------------------------------------------
#pragma mark - Isotropic Alsatian Solve xTBSV, using C-array buffer

template < int ORD >
void alsatian_iso_xtbsvLNN(const int N, const int KD, const real* A, const int lda, real* X)
{
    for ( int j = 0; j < N; ++j )
    {
        real buf[ORD];
        for ( int d = 0; d < ORD; ++d )
        {
            /// X[j] /= pA[0]; // buf = X[j];
            buf[d] = X[d];
            X[d] = buf[d] * A[0];
        }
        const int sup = std::min(N-1-j, KD);
        for ( int ij = 1; ij <= sup; ++ij )
        {
            for ( int d = 0; d < ORD; ++d )
                X[ORD*ij+d] -= buf[d] * A[ij]; // X[i] -= buf * pA[i-j];
        }
        A += lda;
        X += ORD;
    }
}

template < int ORD >
void alsatian_iso_xtbsvLTN(const int N, const int KD, const real* A, const int lda, real* X)
{
    A += N * lda;
    X += N * ORD;
    for ( int j = N-1; j >= 0; --j )
    {
        A -= lda;
        X -= ORD;
        real buf[ORD];
        for ( int d = 0; d < ORD; ++d )
            buf[d] = X[d] * A[0];
        const int sup = std::min(N-1-j, KD);
        for ( int ij = sup; ij > 0; --ij )
        {
            for ( int d = 0; d < ORD; ++d )
                buf[d] -= A[ij] * X[ORD*ij+d];
        }
        for ( int d = 0; d < ORD; ++d )
            X[d] = buf[d];
    }
}


//------------------------------------------------------------------------------
#pragma mark - 3D-Isotropic Symmetric Banded Solve xTBSV for KD==2

#if defined(__AVX__)
/// specialized version for KD==2 and ORD==3
void alsatian_iso3_xtbsvLNN2K_AVX(const int N, const double* pA, const int lda, double* pX)
{
    constexpr int ORD = 3;
    const double*const end = pA + (N-2) * lda;
    vec4 x1 = loadu4(pX); //may load garbage
    vec4 x2 = loadu4(pX+ORD); //may load garbage @ pX+7
    while ( pA < end ) // for ( int j = 0; j < N-2; ++j )
    {
#if DEVELOP_XTBSV
        vec4 x0 = x1;
        vec4 aa = broadcast1(pA);
        x1 = fnmadd4(mul4(aa, broadcast1(pA+1)), x0, x2);
        x2 = fnmadd4(mul4(aa, broadcast1(pA+2)), x0, loadu4(pX+2*ORD));
        storeu4(pX, mul4(aa, x0));
#else
        vec4 x0 = mul4(broadcast1(pA), x1);      // x1 = loadu4(pX);
        x1 = fnmadd4(broadcast1(pA+1), x0, x2);  // x2 = loadu4(pX+ORD);
        x2 = fnmadd4(broadcast1(pA+2), x0, loadu4(pX+2*ORD));
        storeu4(pX, x0); //only 3 values saved: last one will be overwritten
#endif
        pA += lda;
        pX += ORD;
    }
    if ( N >= 2 ) // j = N-2
    {
        vec4 x0 = mul4(broadcast1(pA), x1);
        x1 = fnmadd4(broadcast1(pA+1), x0, x2);
        storeu4(pX, x0);
        pA += lda;
        pX += ORD;
    }
    if ( N >= 1 ) // j = N-1
    {
        vec4 x0 = mul4(broadcast1(pA), x1);
        store3(pX, x0);
        //pA += lda;
        //pX += ORD;
    }
}
#elif USE_SIMD
/// specialized SSE & ARM's NEON version for KD==2 and ORD==3 (21.08.2022)
void alsatian_iso3_xtbsvLNN2K_SIMD(const int N, const double* A, const int lda, double* X)
{
    constexpr int ORD = 3;
    const double * end = X + ORD * ( N - 2 );
    vec2 x2 = load2(X);   //may load garbage if N == 0
    vec2 z2 = load1(X+2); //may load garbage if N == 0
    vec2 x4 = load2(X+3); //may load garbage if N < 1
    vec2 z4 = load1(X+5); //may load garbage if N < 1
    while ( X < end ) // ( int j = 2; j < N; ++j )
    {
        vec2 x0 = x2; // x1 = load2(X);
        vec2 z0 = z2;
        vec2 aa = load2(A+1);
        vec2 a1 = unpacklo2(aa, aa);
        vec2 a2 = unpackhi2(aa, aa);
        aa = loaddup2(A);
        x2 = fnmadd2(a1, x0, x4); // x4 = load2(X+2);
        z2 = fnmadd1(a1, z0, z4);
        x4 = fnmadd2(a2, x0, load2(X+6));
        z4 = fnmadd1(a2, z0, load1(X+8));
        store2(X, mul2(aa, x0));
        store1(X+2, mul1(aa, z0));
        A += lda;
        X += ORD;
    }
    if ( N > 1 )
    {
        vec2 aa = load2(A);
        vec2 a0 = unpacklo2(aa, aa);
        store2(X, mul2(a0, x2));
        store1(X+2, mul2(a0, z2));
        vec2 a1 = unpackhi2(aa, aa);
        x2 = fnmadd2(a1, x2, x4); // x4 = load2(X+2);
        z2 = fnmadd1(a1, z2, z4); // x4 = load2(X+2);
        vec2 aL = loaddup2(A+lda);
        store2(X+3, mul2(aL, x2));
        store1(X+5, mul2(aL, z2));
    }
    else if ( N > 0 )
    {
        vec2 aa = loaddup2(A);
        store2(X, mul2(aa, x2));
        store1(X+2, mul1(aa, z2));
    }
}
#endif

#if defined(__AVX__)
/// specialized version for KD==2 and ORD==3
void alsatian_iso3_xtbsvLTN2K_AVX(const int N, const double* pA, const int lda, double* pX)
{
    constexpr int ORD = 3;
    const double*const end = pA + lda;
    pX += ( N - 1 ) * ORD;
    pA += ( N - 1 ) * lda;
    vec4 x0, x1, x2;
    if ( N <= 2 )
    {
        x1 = setzero4();
        if ( N > 0 ) // j = N-1
        {
            x0 = load3(pX);   // would load garbage in 4th position
            x1 = mul4(broadcast1(pA), x0);
            store3(pX, x1); //storeu4(pX, blend31(x1, x0));
            pA -= lda;
            pX -= ORD;
        }
        x2 = x1;
        if ( N > 1 ) // j = N-2
        {
            x0 = loadu4(pX);
            x0 = fnmadd4(broadcast1(pA+1), x1, x0);
            x1 = mul4(broadcast1(pA), x0);
            storeu4(pX, blend31(x1, x0));
            pA -= lda;
            pX -= ORD;
        }
        return;
    }
    // j = N-1
    x0 = load3(pX);   // would load garbage in 4th position
    x1 = mul4(broadcast1(pA), x0);
    store3(pX, x1); //storeu4(pX, blend31(x1, x0));
    pA -= lda;
    pX -= ORD;
    x2 = x1;
    // j = N-2
    x0 = loadu4(pX);
    x0 = fnmadd4(broadcast1(pA+1), x1, x0);
    x1 = mul4(broadcast1(pA), x0);
    storeu4(pX, blend31(x1, x0));
    pA -= lda;
    pX -= ORD;
    vec4 tt = broadcastX(x1);
    vec4 af = loadu4(pX);
    while ( pA >= end ) // for ( int j = N-3; j > 0; --j )
    {
#if DEVELOP_XTBSV
        vec4 aa = broadcast1(pA);
        x0 = fnmadd4(mul4(aa, broadcast1(pA+2)), x2, mul4(aa, af));
        af = loadu4(pX-ORD);
        x2 = x1;
        x1 = fnmadd4(mul4(aa, broadcast1(pA+1)), x1, x0);
#else
        x0 = fnmadd4(broadcast1(pA+2), x2, af);
        af = loadu4(pX-ORD);
        x2 = x1;
        x0 = fnmadd4(broadcast1(pA+1), x1, x0);
        x1 = mul4(broadcast1(pA), x0);
#endif
        // restore 4th position that was saved in 'tt'
        storeu4(pX, blend31(x1, tt));
        tt = broadcastX(x1); // save 4th position for next round
        pA -= lda;
        pX -= ORD;
    }
    // j = 0
    x0 = fnmadd4(broadcast1(pA+2), x2, af);
    x0 = fnmadd4(broadcast1(pA+1), x1, x0);
    // restore 4th position that was saved in 'tt'
    x0 = blend31(mul4(broadcast1(pA), x0), tt);
    storeu4(pX, x0);
}
#elif USE_SIMD
/// specialized version for KD==2 and ORD==3, using SSE or ARM's NEON (21.08.2022)
void alsatian_iso3_xtbsvLTN2K_SIMD(const int N, const double* A, const int lda, double* X)
{
    constexpr int ORD = 3;
    const double * end = X;
    X += N * ORD;
    A += N * lda;
    vec2 x4, z4, x2, z2;
    if ( N > 0 ) // j = N-1
    {
        A -= lda;
        X -= ORD;
        vec2 aa = loaddup2(A);
        x4 = mul2(aa, load2(X));
        z4 = mul1(aa, load1(X+2));
        store2(X, x4);
        store1(X+2, z4);
    }
    if ( N > 1 ) // j = N-2
    {
        A -= lda;
        X -= ORD;
        vec2 aa = load2(A);
        vec2 a0 = unpacklo2(aa, aa);
        vec2 a1 = unpackhi2(aa, aa);
        vec2 x0 = mul2(a0, load2(X));
        vec2 z0 = mul1(a0, load1(X+2));
        x2 = fnmadd2(a1, x4, x0); // x4 = load2(X+2)
        z2 = fnmadd1(a1, z4, z0);
        store2(X, x2);
        store1(X+2, z2);
    }
    while ( X > end )
    {
        A -= lda;
        X -= ORD;
        vec2 aa = loaddup2(A);
        vec2 x0 = mul2(aa, load2(X));
        vec2 z0 = mul1(aa, load1(X+2));
        aa = load2(A+1);
        vec2 a1 = unpacklo2(aa, aa);
        vec2 a2 = unpackhi2(aa, aa);
        x0 = fnmadd2(a1, x2, x0); // x2 = load2(X+2)
        z0 = fnmadd1(a1, z2, z0);
        x0 = fnmadd2(a2, x4, x0); // x4 = load2(X+4)
        z0 = fnmadd1(a2, z4, z0);
        x4 = x2;
        z4 = z2;
        x2 = x0;
        z2 = z0;
        store2(X, x0);
        store1(X+2, z0);
    }
}
#endif

//------------------------------------------------------------------------------
#pragma mark - Single precision 3D-Isotropic Symmetric Banded Solve xTBSV for KD==2

#if defined(__SSE3__)
/// single precision 3D-specialized version for KD==2 and ORD==3
void alsatian_iso3_xtbsvLNN2K_SSE(const int N, const float* pA, const int lda, float* pX)
{
    constexpr int ORD = 3;
    const float*const end = pA + (N-2) * lda;
    vec4f x1 = loadu4f(pX);     //may load garbage
    vec4f x2 = loadu4f(pX+ORD); //may load garbage @ pX+7
    while ( pA < end ) // for ( int j = 0; j < N-2; ++j )
    {
#if DEVELOP_XTBSV
        vec4f x0 = x1;
        vec4f aa = broadcast1f(pA);
        x1 = fnmadd4f(mul4f(aa, broadcast1f(pA+1)), x0, x2);
        x2 = fnmadd4f(mul4f(aa, broadcast1f(pA+2)), x0, loadu4f(pX+2*ORD));
        storeu4f(pX, mul4f(aa, x0));
#else
        vec4f x0 = mul4f(broadcast1f(pA), x1);      // x1 = loadu4(pX);
        x1 = fnmadd4f(broadcast1f(pA+1), x0, x2);  // x2 = loadu4(pX+ORD);
        x2 = fnmadd4f(broadcast1f(pA+2), x0, loadu4f(pX+2*ORD));
        storeu4f(pX, x0);
#endif
        pA += lda;
        pX += ORD;
    }
    if ( N >= 2 ) // j = N-2
    {
        vec4f x0 = mul4f(broadcast1f(pA), x1);
        x1 = fnmadd4f(broadcast1f(pA+1), x0, x2);
        storeu4f(pX, x0);
        pA += lda;
        pX += ORD;
    }
    if ( N >= 1 ) // j = N-1
    {
        vec4f x0 = mul4f(broadcast1f(pA), x1);
        store3f(pX, x0);
        //pA += lda;
        //pX += ORD;
    }
}
#endif

#if defined(__SSE3__)
/// single precision 3D-specialized version for KD==2 and ORD==3
void alsatian_iso3_xtbsvLTN2K_SSE(const int N, const float* pA, const int lda, float* pX)
{
    constexpr int ORD = 3;
    const float*const end = pA + lda;
    pX += ( N - 1 ) * ORD;
    pA += ( N - 1 ) * lda;
    vec4f x0, x1, x2;
    if ( N <= 2 )
    {
        x1 = setzero4f();
        if ( N > 0 ) // j = N-1
        {
            x0 = load3f(pX);   // would load garbage in 4th position
            x1 = mul4f(broadcast1f(pA), x0);
            store3f(pX, x1); // storeu4f(pX, blend31f(x1, x0));
            pA -= lda;
            pX -= ORD;
        }
        if ( N > 1 ) // j = N-2
        {
            x0 = loadu4f(pX);
            x0 = fnmadd4f(broadcast1f(pA+1), x1, x0);
            x1 = mul4f(broadcast1f(pA), x0);
            storeu4f(pX, blend31f(x1, x0));
        }
        return;
    }
    // j = N-1
    x0 = load3f(pX);   // would load garbage in 4th position
    x1 = mul4f(broadcast1f(pA), x0);
    store3f(pX, x1); // storeu4f(pX, blend31f(x1, x0));
    pA -= lda;
    pX -= ORD;
    x2 = x1;
    // j = N-2
    x0 = loadu4f(pX);
    x0 = fnmadd4f(broadcast1f(pA+1), x1, x0);
    x1 = mul4f(broadcast1f(pA), x0);
    storeu4f(pX, blend31f(x1, x0));
    pA -= lda;
    pX -= ORD;
    vec4f tt = broadcastXf(x1);
    vec4f af = loadu4f(pX);
    while ( pA >= end ) // for ( int j = N-3; j > 0; --j )
    {
#if DEVELOP_XTBSV
        vec4f aa = broadcast1f(pA);
        x0 = fnmadd4f(mul4f(aa, broadcast1f(pA+2)), x2, mul4f(aa, af));
        af = loadu4f(pX-ORD);
        x2 = x1;
        x1 = fnmadd4f(mul4f(aa, broadcast1f(pA+1)), x1, x0);
#else
        x0 = fnmadd4f(broadcast1f(pA+2), x2, af);
        af = loadu4f(pX-ORD);
        x2 = x1;
        x0 = fnmadd4f(broadcast1f(pA+1), x1, x0);
        x1 = mul4f(broadcast1f(pA), x0);
#endif
        // restore 4th position that was saved in 'tt'
        storeu4f(pX, blend31f(x1, tt));
        tt = broadcastXf(x1); // save 4th position for next round
        pA -= lda;
        pX -= ORD;
    }
    // j = 0
    x0 = fnmadd4f(broadcast1f(pA+2), x2, af);
    x0 = fnmadd4f(broadcast1f(pA+1), x1, x0);
    // restore 4th position that was saved in 'tt'
    x0 = blend31f(mul4f(broadcast1f(pA), x0), tt);
    storeu4f(pX, x0);
}
#endif

//------------------------------------------------------------------------------
#pragma mark - 2D-Isotropic Symmetric Banded Solve xTBSV for KD==2

#if USE_SIMD
/// specialized version for KD==2 and ORD==2
void alsatian_iso2_xtbsvLNN2K_SIMD(const int N, const double* A, const int lda, double* X)
{
    constexpr int ORD = 2;
    const double * end = X + ORD * ( N - 2 );
    vec2 x2 = load2(X);   //may load garbage if N == 0
    vec2 x4 = load2(X+2); //may load garbage if N < 1
    while ( X < end ) // ( int j = 2; j < N; ++j )
    {
        vec2 x0 = x2; // x1 = load2(X);
        x2 = fnmadd2(loaddup2(A+1), x2, x4); // x4 = load2(X+2);
        x4 = fnmadd2(loaddup2(A+2), x0, load2(X+4));
        store2(X, mul2(loaddup2(A), x0));
        A += lda;
        X += ORD;
    }
    if ( N > 1 )
    {
        store2(X, mul2(loaddup2(A), x2));
        x2 = fnmadd2(loaddup2(A+1), x2, x4); // x4 = load2(X+2);
        store2(X+2, mul2(loaddup2(A+lda), x2));
    }
    else if ( N > 0 )
    {
        store2(X, mul2(loaddup2(A), x2));
    }
}


/// specialized version for KD==2 and ORD==2
void alsatian_iso2_xtbsvLTN2K_SIMD(const int N, const double* A, const int lda, double* X)
{
    constexpr int ORD = 2;
    const double * end = X;
    X += N * ORD;
    A += N * lda;
    vec2 x4, x2;
    if ( N > 0 ) // j = N-1
    {
        A -= lda;
        X -= ORD;
        x4 = mul2(loaddup2(A), load2(X));
        store2(X, x4);
    }
    if ( N > 1 ) // j = N-2
    {
        A -= lda;
        X -= ORD;
        vec2 x0 = mul2(loaddup2(A), load2(X));
        x2 = fnmadd2(loaddup2(A+1), x4, x0); // x4 = load2(X+2)
        store2(X, x2);
    }
    while ( X > end )
    {
        A -= lda;
        X -= ORD;
        vec2 x0 = mul2(loaddup2(A), load2(X));
        x0 = fnmadd2(loaddup2(A+1), x2, x0); // x2 = load2(X+2)
        x0 = fnmadd2(loaddup2(A+2), x4, x0); // x4 = load2(X+4)
        x4 = x2;
        x2 = x0;
        store2(X, x0);
    }
}

#endif

//------------------------------------------------------------------------------
#pragma mark - 1D-Isotropic Symmetric Banded Solve xTBSV for KD==2

/// specialized version for KD==2 and ORD==1
void alsatian_xtbsvLNN2K(const int N, const real* A, const int lda, real* X)
{
    const double * end = X + ( N - 2 );
    real x2 = X[0]; //may load garbage
    real x4 = X[1]; //may load garbage
    while ( X < end )
    {
        real x0 = x2;
        x2 = x4 - A[1] * x2;
        x4 = X[2] - A[2] * x0;
        X[0] = A[0] * x0;
        A += lda;
        X += 1;
    }
    if ( N > 1 ) // j = N-2
    {
        X[0] = A[0] * x2;
        x2 = x4 - A[1] * x2;
        A += lda;
        X += 1;
    }
    if ( N > 0 ) // j = N-1
    {
        X[0] = A[0] * x2;
        //A += lda;
        //X += 1;
    }
}


/// specialized version for KD==2 and ORD==1
void alsatian_xtbsvLTN2K(const int N, const real* A, const int lda, real* X)
{
    const double * end = X;
    X += N;
    A += N * lda;
    real x4, x2;
    if ( N > 0 ) // j = N-1
    {
        A -= lda;
        X -= 1;
        x4 = A[0] * X[0];
        X[0] = x4;
    }
    if ( N > 1 ) // j = N-2
    {
        A -= lda;
        X -= 1;
        real x0 = A[0] * X[0];
        x2 = x0 - A[1] * x4;
        X[0] = x2;
    }
    while ( X > end )
    {
        A -= lda;
        X -= 1;
        real x0 = A[0] * X[0];
        x0 = x0 - A[1] * x2;
        x0 = x0 - A[2] * x4;
        x4 = x2;
        x2 = x0;
        X[0] = x0;
    }
}

//------------------------------------------------------------------------------
#pragma mark - LAPACK-STYLE DPBTRS


void lapack_xpbtrs(char UPLO, int N, int KD, int NRHS, real const* AB, int LDAB, real* B, int LDB, int* INFO)
{
    *INFO = 0;
    if ( UPLO == 'U' )
    {
        for ( int i = 0; i < NRHS; ++i )
        {
            blas_xtbsvUT<'N'>(N, KD, AB, LDAB, B+i*LDB, 1); //blas::xtbsv('U', 'T', 'N', N, KD, AB, LDAB, B+i*LDB, 1);
            blas_xtbsvUN<'N'>(N, KD, AB, LDAB, B+i*LDB, 1); //blas::xtbsv('U', 'N', 'N', N, KD, AB, LDAB, B+i*LDB, 1);
        }
    }
    else if ( UPLO == 'L' )
    {
        for ( int i = 0; i < NRHS; ++i )
        {
            blas_xtbsvLN<'N'>(N, KD, AB, LDAB, B+i*LDB); //blas::xtbsv('L', 'N', 'N', N, KD, AB, LDAB, B+i*LDB, 1);
            blas_xtbsvLT<'N'>(N, KD, AB, LDAB, B+i*LDB); //blas::xtbsv('L', 'T', 'N', N, KD, AB, LDAB, B+i*LDB, 1);
        }
    }
    else
        *INFO = 1;
}


template < int ORD >
void iso_xpbtrs(char UPLO, int N, int KD, real const* AB, int LDAB, real* B, int LDB, int* INFO)
{
    *INFO = 0;
    if ( UPLO == 'U' )
    {
        for ( int d = 0; d < ORD; ++d )
        {
            blas_xtbsvUT<'N'>(N, KD, AB, LDAB, B+d, ORD); //blas::xtbsv('U', 'T', 'N', N, KD, AB, LDAB, B+d, ORD);
            blas_xtbsvUN<'N'>(N, KD, AB, LDAB, B+d, ORD); //blas::xtbsv('U', 'N', 'N', N, KD, AB, LDAB, B+d, ORD);
        }
    }
    else if ( UPLO == 'L' )
    {
        for ( int d = 0; d < ORD; ++d )
        {
            blas_xtbsvLN<'N'>(N, KD, AB, LDAB, B+d, ORD); //blas::xtbsv('L', 'N', 'N', N, KD, AB, LDAB, B+d, ORD);
            blas_xtbsvLT<'N'>(N, KD, AB, LDAB, B+d, ORD); //blas::xtbsv('L', 'T', 'N', N, KD, AB, LDAB, B+d, ORD);
        }
    }
    else
        *INFO = 1;
}


template < int ORD >
void alsatian_iso_xpbtrs(char UPLO, int N, int KD, real const* AB, int LDAB, real* B, int, int* INFO)
{
    *INFO = 0;
    if ( UPLO == 'U' )
    {
        ABORT_NOW("unfinished alsatian_xpbtrs('U', ...)");
        //alsatian_iso_xtbsvUTN<ORD>(N, KD, AB, LDAB, B);
        //alsatian_iso_xtbsvUNN<ORD>(N, KD, AB, LDAB, B);
    }
    else if ( UPLO == 'L' )
    {
        alsatian_iso_xtbsvLNN<ORD>(N, KD, AB, LDAB, B);
        alsatian_iso_xtbsvLTN<ORD>(N, KD, AB, LDAB, B);
    }
    else
        *INFO = 1;
}


template < int ORD >
void alsatian_iso_xpbtrsL(const int N, real const* AB, int LDAB, real* B)
{
    /* use routines for KD=2, and interleaved vectors of size `ORD*N` */
    if ( ORD == 3 )
    {
#if defined(__AVX__) && REAL_IS_DOUBLE
        alsatian_iso3_xtbsvLNN2K_AVX(N, AB, LDAB, B);
        alsatian_iso3_xtbsvLTN2K_AVX(N, AB, LDAB, B);
#elif USE_SIMD && REAL_IS_DOUBLE
        alsatian_iso3_xtbsvLNN2K_SIMD(N, AB, LDAB, B);
        alsatian_iso3_xtbsvLTN2K_SIMD(N, AB, LDAB, B);
#elif defined(__SSE3__) && !REAL_IS_DOUBLE
        alsatian_iso3_xtbsvLNN2K_SSE(N, AB, LDAB, B);
        alsatian_iso3_xtbsvLTN2K_SSE(N, AB, LDAB, B);
#else
        alsatian_iso_xtbsvLNN<3>(N, 2, AB, LDAB, B);
        alsatian_iso_xtbsvLTN<3>(N, 2, AB, LDAB, B);
#endif
    }
    else if ( ORD == 2 )
    {
#if USE_SIMD && REAL_IS_DOUBLE
        alsatian_iso2_xtbsvLNN2K_SIMD(N, AB, LDAB, B);
        alsatian_iso2_xtbsvLTN2K_SIMD(N, AB, LDAB, B);
#else
        alsatian_iso_xtbsvLNN<2>(N, 2, AB, LDAB, B);
        alsatian_iso_xtbsvLTN<2>(N, 2, AB, LDAB, B);
#endif
    }
    else if ( ORD == 1 )
    {
        alsatian_xtbsvLNN2K(N, AB, LDAB, B);
        alsatian_xtbsvLTN2K(N, AB, LDAB, B);
    }
    else
        ABORT_NOW("unexpected dimension!");
}

//#include "vecprint.h"

template < int KD >
void alsatian_xpbtrsLK(const int N, real const* AB, int LDAB, real* B)
{
#if 0
    // Comparing two implementations
    real* tmp = new_real(N);
    copy_real(N, B, tmp);
    alsatian_xtbsvLNNK<KD>(N, AB, LDAB, B);
    alsatian_xtbsvLNN(N, KD, AB, LDAB, tmp);
    VecPrint::print("t", N, tmp, 5);
    VecPrint::print("L", N, B, 5);
    copy_real(N, B, tmp);
    alsatian_xtbsvLTNK<KD>(N, AB, LDAB, B);
    alsatian_xtbsvLTN(N, KD, AB, LDAB, tmp);
    VecPrint::print("-", N, tmp, 5);
    VecPrint::print("L", N, B, 5);
    printf("\n");
    free_real(tmp);
#else
    if ( KD == 6 )
    {
#if REAL_IS_DOUBLE && USE_SIMD
        alsatian_xtbsvLNN6K_SSE(N, AB, LDAB, B);
        alsatian_xtbsvLTN6K_SSE(N, AB, LDAB, B);
#else
        alsatian_xtbsvLNN6K(N, AB, LDAB, B);
        alsatian_xtbsvLTN6K(N, AB, LDAB, B);
#endif
    }
    else
    {
        alsatian_xtbsvLNNK<KD>(N, AB, LDAB, B);
        alsatian_xtbsvLTNK<KD>(N, AB, LDAB, B);
    }
#endif
}

#endif
