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
            const real * pA = A + j * lda;
            if ( diag == 'N' ) X[jx] /= pA[0];
            else if ( diag == 'I' ) X[jx] *= pA[0];
            real tmp = X[jx];
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
            if ( diag == 'N' ) X[j] /= pA[0];
            else if ( diag == 'I' ) X[j] *= pA[0];
            real tmp = X[j];
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
#pragma mark - Alsatian factorization xPBTF2 that inverse the diagonal terms

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
 
 SUBROUTINE DPBTF2( UPLO, N, KD, AB, LDAB, INFO )
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
            *INFO = J;
            return;
        }
        dia = real(1) / std::sqrt(dia);
        AB[0] = dia;
        // Compute elements J+1:J+KN of column J
        int KN = std::min(KD, N-1-J);  // N-1-J < KD iff J >= N-KD
        // scale off diagonal terms in column:
        for ( int K = 1; K <= KN; ++K )
            AB[K] *= dia;
        // update the trailing submatrix within the band.
        real* A = AB + LDAB; // next column of AB
        for ( int K = 0; K < KN; ++K )
        {
            real tmp = AB[K+1];
            for ( int I = K; I < KN; ++I )
                A[I] -= AB[I+1] * tmp;
            A += KLD;
        }
        AB += LDAB;
    }
    *INFO = 0;
}

/**
 This calls the standard lapack::pbtf2()
 and then *** inverts *** the diagonal terms
 
 SUBROUTINE DPBTF2( UPLO, N, KD, AB, LDAB, INFO )
*/
template < int KD >
void alsatian_xpbtf2L(const int N, real* AB, const int LDAB, int* INFO)
{
    /*
    lapack::xpbtf2('L', N, KD, AB, LDAB, INFO);
    if ( 0 == *INFO )
    {
        for ( int u = 0; u < N; ++u )
            AB[LDAB*u] = 1.0 / AB[LDAB*u];
    }
     */
    assert_true( LDAB > KD );
    const int KLD = std::max(1, LDAB-1);
    const int nkd = std::max(0, N-KD);
    //Compute the Cholesky factorization A = L*L**T.
    for ( int J = 0; J < nkd; ++J )
    {
        //Compute L(J,J) and test for non-positive-definiteness.
        real dia = AB[0];
        if ( dia <= 0 )
        {
            *INFO = J;
            return;
        }
        /* In the Alsatian version of the factorization, we invert
         the diagonal term, such that only multiplications are needed
         for applying the factorization to a vector. */
        dia = std::sqrt(dia) * (real(1) / dia);
        AB[0] = dia;
        /* Compute elements J+1:J+KN of column J and update the
         trailing submatrix within the band.*/
        for ( int K = 1; K <= KD; ++K )
            AB[K] *= dia;
        real* nA = AB + LDAB; // next column of AB
        for ( int K = 0; K < KD; ++K )
        {
            real tmp = AB[K+1];
            for ( int I = K; I < KD; ++I )
                nA[I] -= AB[I+1] * tmp;
            nA += KLD;
        }
        AB += LDAB;
    }
    // process remaining columns with < KD terms
    for ( int J = nkd; J < N; ++J )
    {
        //Compute L(J,J) and test for non-positive-definiteness.
        real dia = AB[0];
        if ( dia <= 0 )
        {
            *INFO = J;
            return;
        }
        /* In the Alsatian version of the factorization, we invert
         the diagonal term, such that only multiplications are needed
         for applying the factorization to a vector. */
        dia = std::sqrt(dia) * (real(1) / dia);
        AB[0] = dia;
        /* Compute elements J+1:J+KN of column J and update the
         trailing submatrix within the band.*/
        int KN = std::min(KD, N-1-J); // always N-1-J
        for ( int K = 1; K <= KN; ++K )
            AB[K] *= dia;
        real* nA = AB + LDAB; // next column of AB
        for ( int K = 0; K < KN; ++K )
        {
            real tmp = AB[K+1];
            for ( int I = K; I < KN; ++I )
                nA[I] -= AB[I+1] * tmp;
            nA += KLD;
        }
        AB += LDAB;
    }
    *INFO = 0;
}

//------------------------------------------------------------------------------
#pragma mark - Alsatian DTBSV with various KD (inverted diagonal terms)

#if 0
//The original LAPACK implementation load and store X at every operation
void alsatian_xtbsvLNN(const int N, const int KD, const real* A, const int lda, real* X)
{
    for ( int j = 0; j < N; ++j )
    {
        const real tmp = X[j] * A[j];
        X[j] = tmp;
        const int sup = std::min(N-1, j+KD);  // sup=N-1 if N-1<j+KD : j >= N-KD
        for ( int i = j + 1; i <= sup; ++i )
            X[i] -= tmp * A[i];
        A += lda - 1;
    }
}

void alsatian_xtbsvLTN(const int N, const int KD, const real* A, const int lda, real* X)
{
    A += (N-1)*(lda-1);
    for ( int j = N-1; j >= 0; --j )
    {
        real tmp = X[j];
        const int sup = std::min(N-1, j+KD);  // sup=N-1 if N-1<j+KD : j >= N-KD
        for ( int i = j + 1; i <= sup; ++i )
        //for ( int i = sup; i > j; --i )
            tmp -= A[i] * X[i];
        X[j] = tmp * A[j];
        A -= lda - 1;
    }
}
#endif

void alsatian_xtbsvLNN(const int N, const int KD, const real* A, const int lda, real* X)
{
    const int nkd = std::max(0, N-KD);
    // general case:
    for ( int j = 0; j < nkd; ++j )
    {
        const real tmp = X[j] * A[0];
        X[j] = tmp;
        for ( int i = 1; i <= KD; ++i )
            X[i+j] -= tmp * A[i];
        A += lda;
    }
    // process truncated cases:
    for ( int j = nkd; j < N; ++j )
    {
        const real tmp = X[j] * A[0];
        X[j] = tmp;
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
        real tmp = X[j];
        for ( int i = j + 1; i < N; ++i )
        //for ( int i = N-1; i > j; --i )
            tmp -= A[i-j] * X[i];
        X[j] = A[0] * tmp;
        A -= lda;
    }
    // general case, downward:
    for ( int j = nkd-1; j >= 0; --j )
    {
        real tmp = X[j];
        for ( int i = 1; i <= KD; ++i )
        //for ( int i = KD; i > 0; --i )
            tmp -= A[i] * X[i+j];
        X[j] = A[0] * tmp;
        A -= lda;
    }
}


template < int KD >
void alsatian_xtbsvLNNK(const int N, const real* A, const int lda, real* X)
{
    assert_true( lda > KD );
    const int nkd = std::max(0, N-KD);
    real buf[KD] = { 0 };
    if ( 0 < N )
    {
        // 'general case' with j = 0, initializing buf[]
        real t = X[0] * A[0];
        X[0] = t;
        for ( int i = 1; i <= std::min(N,KD); ++i )
            buf[i-1] = X[i] - t * A[i];
        A += lda;
        X += 1;
    }
    // general case:
    for ( int j = 1; j < nkd; ++j )
    {
        /*
        const real t = X[j] * A[0];
        X[j] = t;
        for ( int i = 1; i <= KD; ++i )
            X[i+j] -= t * A[i];
         */
        real t = buf[0] * A[0];
        X[0] = t;
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
        real t = buf[0] * A[0];
        X[0] = t;
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
        real tmp = X[0];
/*
        for ( int i = j + 1; i < N; ++i )
            tmp -= A[i-j] * X[i];
*/
        for ( int i = j + 1; i < N; ++i )
            tmp -= A[i-j] * buf[i-nkd]; // buf[] is X[i];
        buf[j-nkd] = A[0] * tmp;
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
        real tmp = X[0] - A[KD] * buf[KD-1];
        for ( int i = KD-1; i > 0; --i )
            buf[i] = buf[i-1];
        for ( int i = 1; i < KD; ++i )
            tmp -= A[i] * buf[i];
        buf[0] = A[0] * tmp;
        X[0] = buf[0];
        A -= lda;
        X -= 1;
    }
}


//------------------------------------------------------------------------------
#pragma mark - Specialized versions for KD==6


/// Specialized version for KD==6
void alsatian_xtbsvLNN6(const int N, const real* A, const int lda, real* X)
{
    const real * end = X + N - 6;
    //constexpr int KD = 6;
    assert_true( lda > 6 );
    assert_true( N >= 6 );
    real tt = X[0];
    real aa = A[0];
    X[0] = tt * aa;
    real x0 = X[1] - tt * (aa * A[1]);
    real x1 = X[2] - tt * (aa * A[2]);
    real x2 = X[3] - tt * (aa * A[3]);
    real x3 = X[4] - tt * (aa * A[4]);
    real x4 = X[5] - tt * (aa * A[5]);
    real x5 = X[6] - tt * (aa * A[6]);
    A += lda;
    ++X;
    for ( ; X < end; ++X )
    {
        tt = x0;
        aa = A[0];
        X[0] = x0 * aa;
        x0 = x1 - tt * (aa * A[1]);
        x1 = x2 - tt * (aa * A[2]);
        x2 = x3 - tt * (aa * A[3]);
        x3 = x4 - tt * (aa * A[4]);
        x4 = x5 - tt * (aa * A[5]);
        x5 = X[6] - tt * (aa * A[6]);
        A += lda;
    }
    for ( ; X < end+6; ++X )
    {
        tt = x0;
        aa = A[0];
        X[0] = x0 * aa;
        x0 = x1 - tt * (aa * A[1]);
        x1 = x2 - tt * (aa * A[2]);
        x2 = x3 - tt * (aa * A[3]);
        x3 = x4 - tt * (aa * A[4]);
        x4 = x5 - tt * (aa * A[5]);
        x5 = 0;
        A += lda;
    }
}


/// Specialized version for KD==6
void alsatian_xtbsvLTN6(const int N, const real* A, const int lda, real* X)
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
        real aa = A[0];
        real tt = aa * ( X[0] - x5 * A[6] );
        tt -= x4 * (aa * A[5]);
        x5 = x4;
        x4 = x3;
        tt -= x3 * (aa * A[4]);
        tt -= x2 * (aa * A[3]);
        x3 = x2;
        x2 = x1;
        tt -= x1 * (aa * A[2]);
        tt -= x0 * (aa * A[1]);
        x1 = x0;
        x0 = tt;
        X[0] = tt;
        A -= lda;
    }
}


//------------------------------------------------------------------------------
#pragma mark - Optimized SSE versions for KD==6

#if defined(__SSE3__)
/**
 Optmized version for KD == 6
 Beware: this works assuming that N >= KD, and it will in particular write
 to X[i] for i = 0 ... 6 for any value of N.
*/
void alsatian_xtbsvLNN6SSE(const int N, const double* A, const int lda, double* X)
{
    constexpr int KD = 6;
    assert_true( lda > KD );
    assert_true( N >= KD );
    const double * end = X + N - KD;
    vec2 x0, x2, x4, aa, tt;
    // first iteration with j = 0
    aa = loaddup2(A);
    tt = loaddup2(X);
    store1(X, mul1(aa, tt));
    x0 = fnmadd2(tt, mul2(aa, loadu2(A+1)), loadu2(X+1));
    x2 = fnmadd2(tt, mul2(aa, loadu2(A+3)), loadu2(X+3));
    x4 = fnmadd2(tt, mul2(aa, loadu2(A+5)), loadu2(X+5)); // may load garbage if N <= KD
    A += lda;
    X += 1;
    // general case:
    for ( ; X < end; ++X ) //for ( int j = 1; j < nkd; ++j )
    {
        /*
         By expanding the multiplication, we reduce the dependency chain...
         tt = buf[0];
         aa = A[0];
         X[0] = buf[0] * aa;
         buf[0] = buf[1] - tt * (aa * A[1]);
         buf[1] = buf[2] - tt * (aa * A[2]);
         buf[2] = buf[3] - tt * (aa * A[3]);
         buf[3] = buf[4] - tt * (aa * A[4]);
         buf[4] = buf[5] - tt * (aa * A[5]);
         buf[5] =  X[KD] - tt * (aa * A[6]);
        */
#if DEVELOP_XTBSV
        aa = loaddup2(A);
        tt = unpacklo2(x0, x0);
        store1(X, mul1(aa, x0));
        x0 = catshift(x0, x2);
        x2 = catshift(x2, x4);
        x4 = catshift(x4, load1(X+KD));
        x0 = fnmadd2(tt, mul2(aa, loadu2(A+1)), x0);
        x2 = fnmadd2(tt, mul2(aa, loadu2(A+3)), x2);
        x4 = fnmadd2(tt, mul2(aa, loadu2(A+5)), x4);
#else
        /* This is a direct translation of the C-code above, but the dependency
        between each iteration is longer than above, which is limiting throughput:
        From x0 to x0, two sequential multiplications are needed */
        vec2 t = mul1(load1(A), x0);
        store1(X, t);
        x0 = catshift(x0, x2);
        x2 = catshift(x2, x4);
        x4 = catshift(x4, load1(X+KD));
        tt = unpacklo2(t, t);
        x0 = fnmadd2(tt, loadu2(A+1), x0);
        x2 = fnmadd2(tt, loadu2(A+3), x2);
        x4 = fnmadd2(tt, loadu2(A+5), x4);
#endif
        A += lda;
    }
    /*
     The ending sequence avoids loading elements of X beyond X+N, and elements
     of A in the last column that should normally be equal to zero, but otherwise
     it performs the same calculation than the normal iteration six times.
     */
    // process truncated case: j = N - KD
    aa = loaddup2(A);
    tt = unpacklo2(x0, x0);
    vec2 yy = mul1(aa, x0);
    x0 = catshift(x0, x2);
    x2 = catshift(x2, x4);
    x4 = unpackhi2(x4, setzero2());
    x0 = fnmadd2(tt, mul2(aa, loadu2(A+1)), x0);
    x2 = fnmadd2(tt, mul2(aa, loadu2(A+3)), x2);
    x4 = fnmadd1(tt, mul1(aa, load1(A+5)), x4);
    A += lda;
    // process truncated case: j = N - 5
    aa = loaddup2(A);
    tt = unpacklo2(x0, x0);
    storeu2(X, unpacklo2(yy, mul1(aa, x0)));
    x0 = catshift(x0, x2);
    x2 = catshift(x2, x4);
    x0 = fnmadd2(tt, mul2(aa, loadu2(A+1)), x0);
    x2 = fnmadd2(tt, mul2(aa, loadu2(A+3)), x2);
    A += lda;
    // process truncated case: j = N - 4
    aa = loaddup2(A);
    tt = unpacklo2(x0, x0);
    yy = mul1(aa, x0);
    x0 = catshift(x0, x2);
    x2 = catshift(x2, setzero2());
    x0 = fnmadd2(tt, mul2(aa, loadu2(A+1)), x0);
    x2 = fnmadd1(tt, mul1(aa, load1(A+3)), x2);
    A += lda;
    // process truncated case: j = N - 3
    aa = loaddup2(A);
    tt = unpacklo2(x0, x0);
    storeu2(X+2, unpacklo2(yy, mul1(aa, x0)));
    x0 = catshift(x0, x2);
    x0 = fnmadd2(tt, mul2(aa, loadu2(A+1)), x0);
    A += lda;
    // process last two cases
    aa = load1(A);
    tt = x0;
    yy = mul1(aa, x0);
    if ( N <= KD ) {
        store1(X, yy);
        return;
    }
    x0 = unpackhi2(x0, x0);
    x0 = fnmadd1(tt, mul1(aa, load1(A+1)), x0);
    x0 = mul1(load1(A+lda), x0);
    storeu2(X+4, unpacklo2(yy, x0));
}


/**
 Optmized version for KD == 6
 Beware: this works assuming that N >= KD, and it will in particular write
 to X[i] for i = 0 ... 6 for any value of N.
*/
void alsatian_xtbsvLTN6SSE(const int N, const double* A, const int lda, double* X)
{
    const double * end = X;
    //constexpr int KD = 6;
    assert_true( lda > 6 );
    assert_true( N >= 6 );
    A += ( N - 1 ) * lda;
    X += N - 1;
    vec2 aa, tt;
    vec2 x0 = setzero2();
    vec2 x2 = setzero2();
    vec2 x4 = setzero2();
#if 1
    /*
     The starting sequence avoids loading elements of X beyond X+N, and elements
     of A in the last column that should normally be equal to zero, but otherwise
     it performs the same calculation than the normal iteration six times.
     */
    // truncated round at j = N-1
    x0 = mul1(load1(A), load1Z(X));
    A -= lda;
    // truncated round at j = N-2
    tt = fnmadd1(x0, load1(A+1), load1(X-1));
    tt = mul1(load1(A), tt);
    x0 = unpacklo2(tt, x0);
    storeu2(X-1, x0);
    A -= lda;
    // truncated round at j = N-3
    aa = loaddup2(A);
    tt = mul2(aa, load1Z(X-2));
    x2 = catshift(x0, x2);
    tt = fnmadd2(x0, mul2(aa, loadu2(A+1)), tt);
    tt = add1(tt, swap2(tt));
    x0 = unpacklo2(tt, x0);
    A -= lda;
    // truncated round at j = N-4
    aa = loaddup2(A);
    tt = mul2(aa, load1Z(X-3));
    tt = fnmadd2(x2, mul2(aa, load1Z(A+3)), tt);
    x2 = catshift(x0, x2);
    tt = fnmadd2(x0, mul2(aa, loadu2(A+1)), tt);
    tt = add1(tt, swap2(tt));
    x0 = unpacklo2(tt, x0);
    storeu2(X-3, x0);
    A -= lda;
    // truncated round at j = N-5
    aa = loaddup2(A);
    tt = mul2(aa, load1Z(X-4));
    tt = fnmadd2(x2, mul2(aa, loadu2(A+3)), tt);
    x4 = catshift(x2, x4);
    x2 = catshift(x0, x2);
    tt = fnmadd2(x0, mul2(aa, loadu2(A+1)), tt);
    tt = add1(tt, swap2(tt));
    x0 = unpacklo2(tt, x0);
    A -= lda;
    // truncated round at j = N-6
    aa = loaddup2(A);
    tt = mul2(aa, fnmadd2(x4, load1Z(A+5), load1Z(X-5)));
    tt = fnmadd2(x2, mul2(aa, loadu2(A+3)), tt);
    x4 = catshift(x2, x4);
    x2 = catshift(x0, x2);
    tt = fnmadd2(x0, mul2(aa, loadu2(A+1)), tt);
    tt = add1(tt, swap2(tt));
    x0 = unpacklo2(tt, x0);
    storeu2(X-5, x0);
    A -= lda;
    X -= 6;
#endif
#if 0
    // unrolling does not change the dependency chain...
    for ( A -= lda, --X; X > end; X -= 2 )
    {
        const real* B = A + lda;
        vec2 xx, bb, ab = unpacklo2(load1(A), load1(B));
        aa = loadu2(A+5);
        bb = loadu2(B+5);
        tt = mul2(ab, fnmadd2(x4, unpackhi2(aa, bb), loadu2(X)));
        tt = fnmadd2(catshift(x2, x4), mul2(ab, unpacklo2(aa, bb)), tt);
        aa = loadu2(A+3);
        bb = loadu2(B+3);
        tt = fnmadd2(x2, mul2(ab, unpackhi2(aa, bb)), tt);
        tt = fnmadd2(catshift(x0, x2), mul2(ab, unpacklo2(aa, bb)), tt);
        aa = loadu2(A+1);
        bb = loadu2(B+1);
        tt = fnmadd2(x0, mul2(ab, unpackhi2(aa, bb)), tt);
        vec2 AB = mul2(ab, unpacklo2(aa, bb));
        // set correct value in high position:
        xx = fnmadd2(unpacklo2(x0, x0), AB, tt);
        // set correct values in both positions:
        tt = fnmadd2(catshift(xx, x0), AB, tt);
        x4 = x2;
        x2 = x0;
        x0 = tt;
        storeu2(X, tt);
        A -= lda*2;
    }
    A += lda;
    ++X;
#endif
#if 1
    // this does not use SIMD parallelism, and it is as fast as anything else...
    vec2 x1 = unpackhi2(x0, x0);
    vec2 x3 = unpackhi2(x2, x2);
    vec2 x5 = unpackhi2(x4, x4);
    for ( ; X >= end; --X )
    {
        aa = load1(A);
        vec2 ab = loadu2(A+5);
        tt = mul1(aa, fnmadd1(x5, unpackhi2(ab, ab), load1(X)));
        tt = fnmadd1(x4, mul1(aa, unpacklo2(ab, ab)), tt);
        x5 = x4;
        x4 = x3;
        ab = loadu2(A+3);
        tt = fnmadd1(x3, mul1(aa, unpackhi2(ab, ab)), tt);
        tt = fnmadd1(x2, mul1(aa, unpacklo2(ab, ab)), tt);
        x3 = x2;
        x2 = x1;
        ab = loadu2(A+1);
        tt = fnmadd1(x1, mul1(aa, unpackhi2(ab, ab)), tt);
        tt = fnmadd1(x0, mul1(aa, unpacklo2(ab, ab)), tt);
        x1 = x0;
        x0 = tt;
        store1(X, tt);
        A -= lda;
    }
    return;
#endif
    // general case, downward!
    for ( ; X >= end; --X ) //for ( int j = nkd-1; j >= 0; --j )
    {
        /*
         // distributing the multiplication to reduce the dependency chain...
         real tmp = A[0] * ( X[0] - A[6] * buf[5] );
         tmp -= A[0] * ( A[5] * buf[4] );
         tmp -= (A[0] * A[4]) * buf[3];
         tmp -= (A[0] * A[3]) * buf[2];
         tmp -= (A[0] * A[2]) * buf[1];
         tmp -= (A[0] * A[1]) * buf[0];
         for ( int i = KD-1; i > 0; --i )
             buf[i] = buf[i-1];
         buf[0] = tmp;
         X[0] = buf[0];
         */
#if DEVELOP_XTBSV
        aa = loaddup2(A);
        tt = mul2(aa, fnmadd2(x4, loadu2(A+5), load1Z(X)));
        //tt = fnmadd2(x4, mul2(aa, loadu2(A+5)), mul2(aa, load1Z(X)));
        tt = fnmadd2(x2, mul2(aa, loadu2(A+3)), tt);
        x4 = catshift(x2, x4);
        x2 = catshift(x0, x2);
        tt = fnmadd2(x0, mul2(aa, loadu2(A+1)), tt);
        tt = add1(tt, swap2(tt));
#else
        /* This is a direct translation of the C-code above, but the dependency
        between each iteration is longer than above, which is limiting throughput:
        From x0 to x0, four sequential multiplications are needed */
        tt = fnmadd2(x4, loadu2(A+5), load1Z(X));
        tt = fnmadd2(x2, loadu2(A+3), tt);
        x4 = catshift(x2, x4);
        x2 = catshift(x0, x2);
        tt = fnmadd2(x0, loadu2(A+1), tt);
        tt = add1(tt, swap2(tt));
        tt = mul1(load1(A), tt);
#endif
        x0 = unpacklo2(tt, x0);
        store1(X, tt);
        A -= lda;
    }
}
#endif

//------------------------------------------------------------------------------
#pragma mark - Isotropic Alsatian Solve xTBSV with KD==2


#if ( 1 )

/// this version is quite fast...
template < int ORD >
void alsatian_xtbsvLNN(const int N, const int KD, const real* A, const int lda, real* X)
{
    int kx = 0;
    int jx = 0;
    for ( int j = 0; j < N; ++j )
    {
        kx += ORD;
        //if ( X[jx] != 0. )
        {
            real* pX = X + kx; // kx == ORD*(j+1);
            const real * pA = A + j * lda;
            real buf[ORD];
            for ( int d = 0; d < ORD; ++d )
            {
                buf[d] = X[jx+d] * pA[0]; // X[jx] *= pA[0];
                X[jx+d] = buf[d]; //real buf = X[jx];
            }
            const int sup = std::min(N-1-j, KD); // ( N-1-j < KD ) if ( j >= N-KD )
            for ( int ij = 1; ij <= sup; ++ij )
            {
                for ( int d = 0; d < ORD; ++d )
                     pX[d] -= buf[d] * pA[ij];  // X[ix] -= buf * pA[i-j];
                pX += ORD;
            }
        }
        jx += ORD;
    }
}

template < int ORD >
void alsatian_xtbsvLTN(const int N, const int KD, const real* A, const int lda, real* X)
{
    int kx = (N-1) * ORD;
    int jx = kx;
    for ( int j = N-1; j >= 0; --j )
    {
        real buf[ORD];
        for ( int d = 0; d < ORD; ++d )
            buf[d] = X[jx+d]; //real buf = X[jx];
        real* pX = X + kx;
        const real * pA = A + j * lda;
        const int sup = std::min(N-1-j, KD); // ( N-1-j < KD ) if ( j >= N-KD )
        for ( int ij = sup; ij > 0; --ij )
        {
            for ( int d = 0; d < ORD; ++d )
                buf[d] -= pA[ij] * pX[d]; // buf -= pA[i-j] * X[ix];
            pX -= ORD;
        }
        for ( int d = 0; d < ORD; ++d )
            X[jx+d] = buf[d] * pA[0]; //X[jx] = buf * pA[0];
        jx -= ORD;
        if ( j < N-KD )
            kx -= ORD;
    }
}

#else

/*
 */
template < int ORD >
void alsatian_xtbsvLNN(const int N, const int KD, const real* A, const int lda, real* X)
{
    for ( int j = 0; j < N; ++j )
    {
        const real * pA = A + j * lda - j;
        real buf[ORD];
        for ( int d = 0; d < ORD; ++d )
        {
            /// X[j] /= pA[0]; // buf = X[j];
            buf[d] = X[ORD*j+d] * pA[j];
            X[ORD*j+d] = buf[d];
        }
        const int sup = std::min(N-1, j+KD);
        for ( int i = j + 1; i <= sup; ++i )
        {
            for ( int d = 0; d < ORD; ++d )
                X[ORD*i+d] -= buf[d] * pA[i]; // X[i] -= buf * pA[i-j];
        }
    }
}

template < int ORD >
void alsatian_xtbsvLTN(const int N, const int KD, const real* A, const int lda, real* X)
{
    for ( int j = N-1; j >= 0; --j )
    {
        real buf[ORD];
        for ( int d = 0; d < ORD; ++d )
            buf[d] = X[ORD*j+d]; // real buf = X[j];
        const real * pA = A + j * lda - j;
        const int sup = std::min(N-1, j+KD);
        for ( int i = sup; i > j; --i )
        {
            for ( int d = 0; d < ORD; ++d )
                buf[d] -= pA[i] * X[ORD*i+d]; //buf -= pA[i-j] * X[i];
        }
        for ( int d = 0; d < ORD; ++d )
            X[ORD*j+d] = buf[d] * pA[j]; // buf /= pA[0]; X[j] = buf;
    }
}

#endif

//------------------------------------------------------------------------------
#pragma mark - Isotropic Dimension-specific Solve xTBSV for KD==2

#if defined(__AVX__)
/// specialized version for KD==2 and ORD==3
void alsatian_xtbsvLNN3(const int N, const double* pA, const int lda, double* pX)
{
    constexpr int ORD = 3;
    const double*const end = pA + (N-2) * lda;
    vec4 a1 = loadu4(pX); //may load garbage
    vec4 a2 = loadu4(pX+ORD); //may load garbage @ pX+7
    while ( pA < end ) // for ( int j = 0; j < N-2; ++j )
    {
#if DEVELOP_XTBSV
        vec4 a0 = a1;
        vec4 aa = broadcast1(pA);
        a1 = fnmadd4(mul4(aa, broadcast1(pA+1)), a0, a2);
        a2 = fnmadd4(mul4(aa, broadcast1(pA+2)), a0, loadu4(pX+2*ORD));
        storeu4(pX, mul4(aa, a0));
#else
        vec4 a0 = mul4(broadcast1(pA), a1);      // a1 = loadu4(pX);
        a1 = fnmadd4(broadcast1(pA+1), a0, a2);  // a2 = loadu4(pX+ORD);
        a2 = fnmadd4(broadcast1(pA+2), a0, loadu4(pX+2*ORD));
        storeu4(pX, a0);
#endif
        pA += lda;
        pX += ORD;
    }
    if ( N >= 2 ) // j = N-2
    {
        vec4 a0 = mul4(broadcast1(pA), a1);
        a1 = fnmadd4(broadcast1(pA+1), a0, a2);
        storeu4(pX, a0);
        pA += lda;
        pX += ORD;
    }
    if ( N >= 1 ) // j = N-1
    {
        vec4 a0 = mul4(broadcast1(pA), a1);
        store3(pX, a0);
        //pA += lda;
        //pX += ORD;
    }
}
#endif

#if defined(__AVX__)
/// specialized version for KD==2 and ORD==3
void alsatian_xtbsvLTN3(const int N, const double* pA, const int lda, double* pX)
{
    constexpr int ORD = 3;
    const double*const end = pA + lda;
    pX += ( N - 1 ) * ORD;
    pA += ( N - 1 ) * lda;
    vec4 a0, a1, a2;
    if ( N <= 2 )
    {
        a1 = setzero4();
        if ( N > 0 ) // j = N-1
        {
            a0 = load3(pX);   // would load garbage in 4th position
            a1 = mul4(broadcast1(pA), a0);
            store3(pX, a1); //storeu4(pX, blend31(a1, a0));
            pA -= lda;
            pX -= ORD;
        }
        a2 = a1;
        if ( N > 1 ) // j = N-2
        {
            a0 = loadu4(pX);
            a0 = fnmadd4(broadcast1(pA+1), a1, a0);
            a1 = mul4(broadcast1(pA), a0);
            storeu4(pX, blend31(a1, a0));
            pA -= lda;
            pX -= ORD;
        }
        return;
    }
    // j = N-1
    a0 = load3(pX);   // would load garbage in 4th position
    a1 = mul4(broadcast1(pA), a0);
    store3(pX, a1); //storeu4(pX, blend31(a1, a0));
    pA -= lda;
    pX -= ORD;
    a2 = a1;
    // j = N-2
    a0 = loadu4(pX);
    a0 = fnmadd4(broadcast1(pA+1), a1, a0);
    a1 = mul4(broadcast1(pA), a0);
    storeu4(pX, blend31(a1, a0));
    pA -= lda;
    pX -= ORD;
    vec4 tt = broadcastX(a1);
    vec4 af = loadu4(pX);
    while ( pA >= end ) // for ( int j = N-3; j > 0; --j )
    {
#if DEVELOP_XTBSV
        vec4 aa = broadcast1(pA);
        a0 = fnmadd4(mul4(aa, broadcast1(pA+2)), a2, mul4(aa, af));
        af = loadu4(pX-ORD);
        a2 = a1;
        a1 = fnmadd4(mul4(aa, broadcast1(pA+1)), a1, a0);
#else
        a0 = fnmadd4(broadcast1(pA+2), a2, af);
        af = loadu4(pX-ORD);
        a2 = a1;
        a0 = fnmadd4(broadcast1(pA+1), a1, a0);
        a1 = mul4(broadcast1(pA), a0);
#endif
        // restore 4th position that was saved in 'tt'
        storeu4(pX, blend31(a1, tt));
        tt = broadcastX(a1); // save 4th position for next round
        pA -= lda;
        pX -= ORD;
    }
    // j = 0
    a0 = fnmadd4(broadcast1(pA+2), a2, af);
    a0 = fnmadd4(broadcast1(pA+1), a1, a0);
    // restore 4th position that was saved in 'tt'
    a0 = blend31(mul4(broadcast1(pA), a0), tt);
    storeu4(pX, a0);
}
#endif


#if defined(__SSE3__)
/// specialized version for single precision, KD==2 and ORD==3
void alsatian_xtbsvLNN3(const int N, const float* pA, const int lda, float* pX)
{
    constexpr int ORD = 3;
    const float*const end = pA + (N-2) * lda;
    vec4f a1 = loadu4f(pX);     //may load garbage
    vec4f a2 = loadu4f(pX+ORD); //may load garbage @ pX+7
    while ( pA < end ) // for ( int j = 0; j < N-2; ++j )
    {
#if DEVELOP_XTBSV
        vec4f a0 = a1;
        vec4f aa = broadcast1f(pA);
        a1 = fnmadd4f(mul4f(aa, broadcast1f(pA+1)), a0, a2);
        a2 = fnmadd4f(mul4f(aa, broadcast1f(pA+2)), a0, loadu4f(pX+2*ORD));
        storeu4f(pX, mul4f(aa, a0));
#else
        vec4f a0 = mul4f(broadcast1f(pA), a1);      // a1 = loadu4(pX);
        a1 = fnmadd4f(broadcast1f(pA+1), a0, a2);  // a2 = loadu4(pX+ORD);
        a2 = fnmadd4f(broadcast1f(pA+2), a0, loadu4f(pX+2*ORD));
        storeu4f(pX, a0);
#endif
        pA += lda;
        pX += ORD;
    }
    if ( N >= 2 ) // j = N-2
    {
        vec4f a0 = mul4f(broadcast1f(pA), a1);
        a1 = fnmadd4f(broadcast1f(pA+1), a0, a2);
        storeu4f(pX, a0);
        pA += lda;
        pX += ORD;
    }
    if ( N >= 1 ) // j = N-1
    {
        vec4f a0 = mul4f(broadcast1f(pA), a1);
        store3f(pX, a0);
        //pA += lda;
        //pX += ORD;
    }
}
#endif

#if defined(__SSE3__)
/// specialized version for single precision, KD==2 and ORD==3
void alsatian_xtbsvLTN3(const int N, const float* pA, const int lda, float* pX)
{
    constexpr int ORD = 3;
    const float*const end = pA + lda;
    pX += ( N - 1 ) * ORD;
    pA += ( N - 1 ) * lda;
    vec4f a0, a1, a2;
    if ( N <= 2 )
    {
        a1 = setzero4f();
        if ( N > 0 ) // j = N-1
        {
            a0 = load3f(pX);   // would load garbage in 4th position
            a1 = mul4f(broadcast1f(pA), a0);
            store3f(pX, a1); // storeu4f(pX, blend31f(a1, a0));
            pA -= lda;
            pX -= ORD;
        }
        if ( N > 1 ) // j = N-2
        {
            a0 = loadu4f(pX);
            a0 = fnmadd4f(broadcast1f(pA+1), a1, a0);
            a1 = mul4f(broadcast1f(pA), a0);
            storeu4f(pX, blend31f(a1, a0));
        }
        return;
    }
    // j = N-1
    a0 = load3f(pX);   // would load garbage in 4th position
    a1 = mul4f(broadcast1f(pA), a0);
    store3f(pX, a1); // storeu4f(pX, blend31f(a1, a0));
    pA -= lda;
    pX -= ORD;
    a2 = a1;
    // j = N-2
    a0 = loadu4f(pX);
    a0 = fnmadd4f(broadcast1f(pA+1), a1, a0);
    a1 = mul4f(broadcast1f(pA), a0);
    storeu4f(pX, blend31f(a1, a0));
    pA -= lda;
    pX -= ORD;
    vec4f tt = broadcastXf(a1);
    vec4f af = loadu4f(pX);
    while ( pA >= end ) // for ( int j = N-3; j > 0; --j )
    {
#if DEVELOP_XTBSV
        vec4f aa = broadcast1f(pA);
        a0 = fnmadd4f(mul4f(aa, broadcast1f(pA+2)), a2, mul4f(aa, af));
        af = loadu4f(pX-ORD);
        a2 = a1;
        a1 = fnmadd4f(mul4f(aa, broadcast1f(pA+1)), a1, a0);
#else
        a0 = fnmadd4f(broadcast1f(pA+2), a2, af);
        af = loadu4f(pX-ORD);
        a2 = a1;
        a0 = fnmadd4f(broadcast1f(pA+1), a1, a0);
        a1 = mul4f(broadcast1f(pA), a0);
#endif
        // restore 4th position that was saved in 'tt'
        storeu4f(pX, blend31f(a1, tt));
        tt = broadcastXf(a1); // save 4th position for next round
        pA -= lda;
        pX -= ORD;
    }
    // j = 0
    a0 = fnmadd4f(broadcast1f(pA+2), a2, af);
    a0 = fnmadd4f(broadcast1f(pA+1), a1, a0);
    // restore 4th position that was saved in 'tt'
    a0 = blend31f(mul4f(broadcast1f(pA), a0), tt);
    storeu4f(pX, a0);
}
#endif


#if defined(__SSE3__)
/// specialized version for KD==2 and ORD==2
void alsatian_xtbsvLNN2(const int N, const double* pA, const int lda, double* pX)
{
    constexpr int ORD = 2;
    vec2 a1 = load2(pX);     //may load garbage if N == 0
    vec2 a2 = load2(pX+ORD); //may load garbage if N < 1
    for ( int j = 0; j < N-2; ++j )
    {
#if DEVELOP_XTBSV
        // this reduces the dependency path
        vec2 aa = loaddup2(pA);
        vec2 a0 = a1;
        store2(pX, mul2(aa, a1));
        a1 = fnmadd2(mul2(aa, loaddup2(pA+1)), a0, a2);
        a2 = fnmadd2(mul2(aa, loaddup2(pA+2)), a0, load2(pX+2*ORD));
#else
        vec2 a0 = mul2(loaddup2(pA), a1);      // a1 = loadu4(pX);
        a1 = fnmadd2(loaddup2(pA+1), a0, a2);  // a2 = loadu4(pX+ORD);
        a2 = fnmadd2(loaddup2(pA+2), a0, load2(pX+2*ORD));
        store2(pX, a0);
#endif
        pA += lda;
        pX += ORD;
    }
    if ( N >= 2 ) // j = N-2
    {
        vec2 a0 = mul2(loaddup2(pA), a1);
        a1 = fnmadd2(loaddup2(pA+1), a0, a2);
        store2(pX, a0);
        pA += lda;
        pX += ORD;
    }
    if ( N >= 1 ) // j = N-1
    {
        vec2 a0 = mul2(loaddup2(pA), a1);
        store2(pX, a0);
        //pA += lda;
        //pX += ORD;
    }
}


/// specialized version for KD==2 and ORD==2
void alsatian_xtbsvLTN2(const int N, const double* pA, const int lda, double* pX)
{
    constexpr int ORD = 2;
    pX += ( N - 1 ) * ORD;
    pA += ( N - 1 ) * lda;
    vec2 a1 = setzero2();
    if ( N >= 1 ) // j = N-1
    {
        a1 = mul2(loaddup2(pA), load2(pX));
        store2(pX, a1);
        pA -= lda;
        pX -= ORD;
    }
    vec2 a2 = a1;
    if ( N >= 2 ) // j = N-2
    {
        vec2 a0 = fnmadd2(loaddup2(pA+1), a1, load2(pX));
        a1 = mul2(loaddup2(pA), a0);
        store2(pX, a1);
        pA -= lda;
        pX -= ORD;
    }
    for ( int j = N-3; j > 0; --j )
    {
#if DEVELOP_XTBSV
        // this reduces the dependency path
        vec2 aa = loaddup2(pA);
        vec2 a0 = fnmadd2(mul2(aa, loaddup2(pA+2)), a2, mul2(aa, load2(pX)));
        a2 = a1;
        a1 = fnmadd2(mul2(aa, loaddup2(pA+1)), a1, a0);
#else
        vec2 a0 = fnmadd2(loaddup2(pA+2), a2, load2(pX));
        a0 = fnmadd2(loaddup2(pA+1), a1, a0);
        a2 = a1;
        a1 = mul2(loaddup2(pA), a0);
#endif
        store2(pX, a1);
        pA -= lda;
        pX -= ORD;
    }
    if ( N >= 3 ) // j = 0
    {
        vec2 a0 = fnmadd2(loaddup2(pA+2), a2, load2(pX));
        a0 = fnmadd2(loaddup2(pA+1), a1, a0);
        store2(pX, mul2(loaddup2(pA), a0));
    }
}

#endif


/// specialized version for KD==2 and ORD==1
void alsatian_xtbsvLNN1(const int N, const real* pA, const int lda, real* pX)
{
    real a1 = pX[0]; //may load garbage
    real a2 = pX[1]; //may load garbage
    for ( int j = 0; j < N-2; ++j )
    {
#if DEVELOP_XTBSV
        // this reduces the dependency path
        real aa = pA[0];
        real a0 = a1;
        a1 = a2 - ( aa * pA[1] ) * a0;
        a2 = pX[2] - ( aa * pA[2] ) * a0;
        pX[0] = aa * a0;
#else
        real a0 = pA[0] * a1;
        a1 = a2 - pA[1] * a0;
        a2 = pX[2] - pA[2] * a0;
        pX[0] = a0;
#endif
        pA += lda;
        pX += 1;
    }
    if ( N >= 2 ) // j = N-2
    {
        real a0 = pA[0] * a1;
        a1 = a2 - pA[1] * a0;
        pX[0] = a0;
        pA += lda;
        pX += 1;
    }
    if ( N >= 1 ) // j = N-1
    {
        real a0 = pA[0] * a1;
        pX[0] = a0;
        //pA += lda;
        //pX += 1;
    }
}


/// specialized version for KD==2 and ORD==1
void alsatian_xtbsvLTN1(const int N, const real* pA, const int lda, real* pX)
{
    pX += ( N - 1 );
    pA += ( N - 1 ) * lda;
    real a1 = 0;
    if ( N >= 1 ) // j = N-1
    {
        real a0 = pX[0];
        a1 = pA[0] * a0;
        pX[0] = a1;
        pA -= lda;
        pX -= 1;
    }
    real a2 = a1;
    if ( N >= 2 ) // j = N-2
    {
        real a0 = pX[0] - pA[1] * a1;
        a1 = pA[0] * a0;
        pX[0] = a1;
        pA -= lda;
        pX -= 1;
    }
    for ( int j = N-3; j > 0; --j )
    {
#if DEVELOP_XTBSV
        // this reduces the dependency path
        real aa = pA[0];
        real a0 = aa * pX[0] - ( aa * pA[2] ) * a2;
        a2 = a1;
        a1 = a0 - ( aa * pA[1] ) * a1;
#else
        real a0 = pX[0] - pA[2] * a2;
        a0 = a0 - pA[1] * a1;
        a2 = a1;
        a1 = pA[0] * a0;
#endif
        pX[0] = a1;
        pA -= lda;
        pX -= 1;
    }
    if ( N >= 3 ) // j = 0
    {
        real a0 = pX[0] - pA[2] * a2;
        a0 = a0 - pA[1] * a1;
        pX[0] = pA[0] * a0;
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
void alsatian_xpbtrs(char UPLO, int N, int KD, real const* AB, int LDAB, real* B, int, int* INFO)
{
    *INFO = 0;
    if ( UPLO == 'U' )
    {
        ABORT_NOW("unfinished alsatian_xpbtrs('U', ...)");
        //alsatian_xtbsvUTN<ORD>(N, KD, AB, LDAB, B);
        //alsatian_xtbsvUNN<ORD>(N, KD, AB, LDAB, B);
    }
    else if ( UPLO == 'L' )
    {
        alsatian_xtbsvLNN<ORD>(N, KD, AB, LDAB, B);
        alsatian_xtbsvLTN<ORD>(N, KD, AB, LDAB, B);
    }
    else
        *INFO = 1;
}


template < int ORD >
void alsatian_xpbtrsL(const int N, real const* AB, int LDAB, real* B)
{
#if defined(__AVX__) && REAL_IS_DOUBLE
    /* use routines for KD=2, and interleaved vectors of size `ORD*N` */
    if ( ORD == 3 )
    {
        alsatian_xtbsvLNN3(N, AB, LDAB, B);
        alsatian_xtbsvLTN3(N, AB, LDAB, B);
    }
    else if ( ORD == 2 )
    {
        alsatian_xtbsvLNN2(N, AB, LDAB, B);
        alsatian_xtbsvLTN2(N, AB, LDAB, B);
    }
    else if ( ORD == 1 )
    {
        alsatian_xtbsvLNN1(N, AB, LDAB, B);
        alsatian_xtbsvLTN1(N, AB, LDAB, B);
    }
    else
        ABORT_NOW("unexpected dimension!");
#elif defined(__SSE3__) && !REAL_IS_DOUBLE
    if ( ORD == 3 )
    {
        alsatian_xtbsvLNN3(N, AB, LDAB, B);
        alsatian_xtbsvLTN3(N, AB, LDAB, B);
    }
    else
    {
        alsatian_xtbsvLNN<ORD>(N, 2, AB, LDAB, B);
        alsatian_xtbsvLTN<ORD>(N, 2, AB, LDAB, B);
    }
#else
    alsatian_xtbsvLNN<ORD>(N, 2, AB, LDAB, B);
    alsatian_xtbsvLTN<ORD>(N, 2, AB, LDAB, B);
#endif
}

#include "vecprint.h"

template < int KD >
void alsatian_xpbtrsLK(const int N, real const* AB, int LDAB, real* B)
{
#if 0
    // Comparing two implementations
    int S = std::min(N, 16);
    real* tmp = new_real(N);
    copy_real(N, B, tmp);
    alsatian_xtbsvLNNK<KD>(N, AB, LDAB, B);
    alsatian_xtbsvLNN(N, KD, AB, LDAB, tmp);
    printf("\n  t "); VecPrint::print(S, tmp, 5);
    printf("\n  L "); VecPrint::print(S, B, 5);
    copy_real(N, B, tmp);
    alsatian_xtbsvLTNK<KD>(N, AB, LDAB, B);
    alsatian_xtbsvLTN(N, KD, AB, LDAB, tmp);
    printf("\n  - "); VecPrint::print(S, tmp, 5);
    printf("\n  L "); VecPrint::print(S, B, 5);
    printf("\n");
    free_real(tmp);
#else
    if ( KD == 6 )
    {
#if REAL_IS_DOUBLE && defined(__SSE3__)
        alsatian_xtbsvLNN6SSE(N, AB, LDAB, B);
        alsatian_xtbsvLTN6SSE(N, AB, LDAB, B);
#else
        alsatian_xtbsvLNN6(N, AB, LDAB, B);
        alsatian_xtbsvLTN6(N, AB, LDAB, B);
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
