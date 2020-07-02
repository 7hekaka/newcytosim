// Cytosim was created by Francois Nedelec.
// Copyright 2020 Sainsbury Laboratory, Cambridge University

#ifndef XTBSV_H
#define XTBSV_H

#ifdef __AVX__
#  include "simd.h"
#  include "simd_float.h"
#elif defined(__SSE3__)
#  include "simd_float.h"
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
    if ( incX == 1 )
    {
        for (int j = N-1; j >= 0; --j)
        {
            if (X[j] != 0.)
            {
                const real * pA = A + KD + j * lda;
                if ( diag == 'N' ) X[j] /= pA[0];
                else if ( diag == 'I' ) X[j] *= pA[0];
                real temp = X[j];
                const int inf = std::max(0, j-KD);
                for (int i = j - 1; i >= inf; --i)
                    X[i] -= temp * pA[i-j];
            }
        }
    }
    else
    {
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
                real temp = X[jx];
                const int inf = std::max(0, j-KD);
                for (int i = j - 1; i >= inf; --i)
                {
                    X[ix] -= temp * pA[i-j];
                    ix -= incX;
                }
            }
            jx -= incX;
        }
    }
}


template < char diag >
void blas_xtbsvLN(const int N, const int KD, const real* A, const int lda, real* X, const int incX)
{
    if ( incX == 1 )
    {
        for (int j = 0; j < N; ++j)
        {
            if (X[j] != 0.)
            {
                const real * pA = A + j * lda;
                if ( diag == 'N' ) X[j] /= pA[0];
                else if ( diag == 'I' ) X[j] *= pA[0];
                real temp = X[j];
                const int sup = std::min(N-1, j+KD);
                for (int i = j + 1; i <= sup; ++i)
                    X[i] -= temp * pA[i-j];
            }
        }
    }
    else
    {
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
                real temp = X[jx];
                const int sup = std::min(N-1, j+KD);
                for (int i = j + 1; i <= sup; ++i)
                {
                    X[ix] -= temp * pA[i-j];
                    ix += incX;
                }
            }
            jx += incX;
        }
    }
}


template < char diag >
void blas_xtbsvUT(const int N, const int KD, const real* A, const int lda, real* X, const int incX)
{
    if ( incX == 1 )
    {
        for (int j = 0; j < N; ++j)
        {
            real temp = X[j];
            const real * pA = A + KD + j * lda;
            for (int i = std::max(0, j-KD); i < j; ++i)
                temp -= pA[i-j] * X[i];
            if ( diag == 'N' ) temp /= pA[0];
            else if ( diag == 'I' ) temp *= pA[0];
            X[j] = temp;
        }
    }
    else
    {
        int kx = 0;
        if ( incX <= 0 )
            kx = - (N-1) * incX;
        int jx = kx;
        for (int j = 0; j < N; ++j)
        {
            real temp = X[jx];
            int ix = kx;
            const real * pA = A + KD + j * lda;
            for (int i = std::max(0, j-KD); i < j; ++i)
            {
                temp -= pA[i-j] * X[ix];
                ix += incX;
            }
            if ( diag == 'N' ) temp /= pA[0];
            else if ( diag == 'I' ) temp *= pA[0];
            X[jx] = temp;
            jx += incX;
            if (j >= KD)
                kx += incX;
        }
    }
}


template < char diag >
void blas_xtbsvLT(const int N, const int KD, const real* A, const int lda, real* X, const int incX)
{
    if ( incX == 1 )
    {
        for (int j = N-1; j >= 0; --j)
        {
            real temp = X[j];
            const real * pA = A + j * lda;
            const int sup = std::min(N-1, j+KD);
            for (int i = sup; i > j; --i)
                temp -= pA[i-j] * X[i];
            if ( diag == 'N' ) temp /= pA[0];
            else if ( diag == 'I' ) temp *= pA[0];
            X[j] = temp;
        }
    }
    else
    {
        int kx = 0;
        if ( incX > 0 )
            kx = (N-1) * incX;
        int jx = kx;
        for (int j = N-1; j >= 0; --j)
        {
            real temp = X[jx];
            int ix = kx;
            const real * pA = A + j * lda;
            const int sup = std::min(N-1, j+KD);
            for (int i = sup; i > j; --i)
            {
                temp -= pA[i-j] * X[ix];
                ix -= incX;
            }
            if ( diag == 'N' ) temp /= pA[0];
            else if ( diag == 'I' ) temp *= pA[0];
            X[jx] = temp;
            jx -= incX;
            if ( j < N-KD )
                kx -= incX;
        }
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
#pragma mark - ALSATIAN versions that inverse the diagonal elements

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
void  blas_xsyrL(int N, real ALPHA, const real* X, real* A, int LDA)
{
    for ( int J = 0; J < N; ++J )
    {
        if ( X[J] != 0 )
        {
            real temp = ALPHA * X[J];
            for ( int I = J; I < N; ++I )
                A[I] = A[I] + X[I] * temp;
        }
        A += LDA;
    }
}


/**
 This calls the standard lapack::pbtf2()
 and then *** inverts *** the diagonal terms
 
 SUBROUTINE DPBTF2( UPLO, N, KD, AB, LDAB, INFO )
 
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
     dia = 1.0 / sqrt(dia);
     AB[0] = dia;
     // Compute elements J+1:J+KN of column J
     int KN = std::min(KD, N-1-J);  // N-1-J < KD iff J >= N-KD
     // scale off diagonal terms in column:
     for ( int K = 1; K <= KN; ++K )
         AB[K] *= dia;
     // update the trailing submatrix within the band.
     real* A = AB + LDAB; // next column of AB
     real const* X = AB + 1;
     for ( int K = 0; K < KN; ++K )
     {
         real temp = X[K];
         for ( int I = K; I < KN; ++I )
             A[I] -= X[I] * temp;
         A += KLD;
     }
     AB += LDAB;
 }
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
    int KLD = std::max(1, LDAB-1);
    //Compute the Cholesky factorization A = L*L**T.
    for ( int J = 0; J < N-KD; ++J )
    {
        //Compute L(J,J) and test for non-positive-definiteness.
        real dia = AB[0];
        if ( dia <= 0 )
        {
            *INFO = J;
            return;
        }
        dia = 1.0 / sqrt(dia); // inverse the diagonal term!!!
        AB[0] = dia;
        /* Compute elements J+1:J+KN of column J and update the
         trailing submatrix within the band.*/
        {
            for ( int K = 1; K <= KD; ++K )
                AB[K] *= dia;
            real* A = AB + LDAB; // next column of AB
            real const* X = AB + 1;
            for ( int K = 0; K < KD; ++K )
            {
                real temp = X[K];
                for ( int I = K; I < KD; ++I )
                    A[I] -= X[I] * temp;
                A += KLD;
            }
        }
        AB += LDAB;
    }
    // process remaining columns with < KD terms
    for ( int J = N-KD; J < N; ++J )
    {
        //Compute L(J,J) and test for non-positive-definiteness.
        real dia = AB[0];
        if ( dia <= 0 )
        {
            *INFO = J;
            return;
        }
        dia = 1.0 / sqrt(dia);
        AB[0] = dia;
        /* Compute elements J+1:J+KN of column J and update the
         trailing submatrix within the band.*/
        int KN = std::min(KD, N-1-J); // always N-1-J
        {
            for ( int K = 1; K <= KN; ++K )
                AB[K] *= dia;
            real* A = AB + LDAB; // next column of AB
            real const* X = AB + 1;
            for ( int K = 0; K < KN; ++K )
            {
                real temp = X[K];
                for ( int I = K; I < KN; ++I )
                    A[I] -= X[I] * temp;
                A += KLD;
            }
        }
        AB += LDAB;
    }
    *INFO = 0;
}


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
            real temp[ORD];
            for ( int d = 0; d < ORD; ++d )
            {
                temp[d] = X[jx+d] * pA[0]; // X[jx] *= pA[0];
                X[jx+d] = temp[d]; //real temp = X[jx];
            }
            const int sup = std::min(N-1-j, KD); // ( N-1-j < KD ) if ( j >= N-KD )
            for ( int ij = 1; ij <= sup; ++ij )
            {
                for ( int d = 0; d < ORD; ++d )
                     pX[d] -= temp[d] * pA[ij];  // X[ix] -= temp * pA[i-j];
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
        real temp[ORD];
        for ( int d = 0; d < ORD; ++d )
            temp[d] = X[jx+d]; //real temp = X[jx];
        real* pX = X + kx;
        const real * pA = A + j * lda;
        const int sup = std::min(N-1-j, KD); // ( N-1-j < KD ) if ( j >= N-KD )
        for ( int ij = sup; ij > 0; --ij )
        {
            for ( int d = 0; d < ORD; ++d )
                temp[d] -= pA[ij] * pX[d]; // temp -= pA[i-j] * X[ix];
            pX -= ORD;
        }
        for ( int d = 0; d < ORD; ++d )
            X[jx+d] = temp[d] * pA[0]; //X[jx] = temp * pA[0];
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
        real temp[ORD];
        for ( int d = 0; d < ORD; ++d )
        {
            /// X[j] /= pA[0]; // temp = X[j];
            temp[d] = X[ORD*j+d] * pA[j];
            X[ORD*j+d] = temp[d];
        }
        const int sup = std::min(N-1, j+KD);
        for ( int i = j + 1; i <= sup; ++i )
        {
            for ( int d = 0; d < ORD; ++d )
                X[ORD*i+d] -= temp[d] * pA[i]; // X[i] -= temp * pA[i-j];
        }
    }
}

template < int ORD >
void alsatian_xtbsvLTN(const int N, const int KD, const real* A, const int lda, real* X)
{
    for ( int j = N-1; j >= 0; --j )
    {
        real temp[ORD];
        for ( int d = 0; d < ORD; ++d )
            temp[d] = X[ORD*j+d]; // real temp = X[j];
        const real * pA = A + j * lda - j;
        const int sup = std::min(N-1, j+KD);
        for ( int i = sup; i > j; --i )
        {
            for ( int d = 0; d < ORD; ++d )
                temp[d] -= pA[i] * X[ORD*i+d]; //temp -= pA[i-j] * X[i];
        }
        for ( int d = 0; d < ORD; ++d )
            X[ORD*j+d] = temp[d] * pA[j]; // temp /= pA[0]; X[j] = temp;
    }
}

#endif

//------------------------------------------------------------------------------
#pragma mark - DIMENSION-SPECIFIC ALSATIAN DPBTF2

#if defined(__AVX__)

/// specialized version for KD==2 and ORD==3
void alsatian_xtbsvLNN3(const int N, const double* pA, const int lda, double* pX)
{
    const double*const end = pA + (N-2) * lda;
    constexpr int ORD = 3;
    vec4 a1 = loadu4(pX); //may load garbage
    vec4 a2 = loadu4(pX+ORD); //may load garbage
#if ( 0 )
    while ( pA < end-lda )
    {
        vec4 a0 = mul4(broadcast1(pA), a1);      // a1 = loadu4(pX);
        vec4 b1 = fnmadd4(broadcast1(pA+1), a0, a2);  // a2 = loadu4(pX+ORD);
        vec4 b2 = fnmadd4(broadcast1(pA+2), a0, loadu4(pX+2*ORD));
        storeu4(pX, a0);
        vec4 b0 = mul4(broadcast1(pA+lda), b1);
        a1 = fnmadd4(broadcast1(pA+lda+1), b0, b2);
        a2 = fnmadd4(broadcast1(pA+lda+2), b0, loadu4(pX+3*ORD));
        storeu4(pX+ORD, b0);
        pA += 2*lda;
        pX += 2*ORD;
    }
#endif
    while ( pA < end ) // for ( int j = 0; j < N-2; ++j )
    {
        vec4 a0 = mul4(broadcast1(pA), a1);      // a1 = loadu4(pX);
        a1 = fnmadd4(broadcast1(pA+1), a0, a2);  // a2 = loadu4(pX+ORD);
        a2 = fnmadd4(broadcast1(pA+2), a0, loadu4(pX+2*ORD));
        storeu4(pX, a0);
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


/// specialized version for KD==2 and ORD==3
void alsatian_xtbsvLTN3(const int N, const double* pA, const int lda, double* pX)
{
    const double*const end = pA;
    constexpr int ORD = 3;
    const vec4 zero = setzero4();
    pX += ( N - 1 ) * ORD;
    pA += ( N - 1 ) * lda;
    vec4 a1 = zero;
    if ( N >= 1 ) // j = N-1
    {
        vec4 a0 = loadu4(pX);
        a1 = mul4(broadcast1(pA), a0);
        storeu4(pX, blend4(a1, a0, 0b1000));
        a1 = blend4(a1, zero, 0b1000);
        pA -= lda;
        pX -= ORD;
    }
    vec4 a2 = a1;
    if ( N >= 2 ) // j = N-2
    {
        vec4 a0 = loadu4(pX);
        a0 = fnmadd4(broadcast1(pA+1), a1, a0);
        a1 = mul4(broadcast1(pA), a0);
        storeu4(pX, blend4(a1, a0, 0b1000));
        a1 = blend4(a1, zero, 0b1000);
        pA -= lda;
        pX -= ORD;
    }
    vec4 tt = broadcast1(a1);
    vec4 af = loadu4(pX);
    while ( pA >= end ) // for ( int j = N-3; j >= 0; --j )
    {
        vec4 a0 = fnmadd4(broadcast1(pA+2), a2, af);  // a2 = load3(pX+6);
        af = loadu4(pX-ORD);
        a2 = a1;
        a0 = fnmadd4(broadcast1(pA+1), a1, a0);  // a1 = load3(pX+3);
        // restore 4th position that was saved in 'tt'
        a0 = blend4(mul4(broadcast1(pA), a0), tt, 0b1000);
        tt = broadcast1(a0); // save 4th position for next round
        storeu4(pX, a0);
        a1 = blend4(a0, zero, 0b1000);
        pA -= lda;
        pX -= ORD;
    }
}

/*
 void alsatian_xtbsvLTN3(const int N, const real* pA, const int lda, real* pX)
 {
     pX += ( N - 1 ) * 3;
     pA += ( N - 1 ) * lda;
     if ( N >= 1 ) // j = N-1
     {
         store3(pX, mul4(loadu4(pX), broadcast1(pA)));
         pA -= lda;
         pX -= 3;
     }
     if ( N >= 2 ) // j = N-2
     {
         vec4 temp = load3(pX); //real temp = X[jx];
         temp = fnmadd4(broadcast1(pA+1), loadu4(pX+3), temp);
         store3(pX, mul4(temp, broadcast1(pA)));
         pA -= lda;
         pX -= 3;
     }
     for ( int j = N-3; j >= 0; --j )
     {
         vec4 temp = load3(pX); //real temp = X[jx];
         temp = fnmadd4(broadcast1(pA+2), loadu4(pX+6), temp);
         temp = fnmadd4(broadcast1(pA+1), loadu4(pX+3), temp);
         store3(pX, mul4(temp, broadcast1(pA)));
         pA -= lda;
         pX -= 3;
     }
 }
*/
#endif


#if defined(__SSE3__)
/// specialized version for KD==2 and ORD==3
void alsatian_xtbsvLNN3(const int N, const float* pA, const int lda, float* pX)
{
    const float*const end = pA + (N-2) * lda;
    constexpr int ORD = 3;
    vec4f a1 = loadu4f(pX);     //may load garbage
    vec4f a2 = loadu4f(pX+ORD); //may load garbage
#if ( 0 )
    while ( pA+lda < end )
    {
        vec4f a0 = mul4f(broadcast1f(pA), a1);      // a1 = loadu4(pX);
        vec4f b1 = fnmadd4f(broadcast1f(pA+1), a0, a2);  // a2 = loadu4(pX+ORD);
        vec4f b2 = fnmadd4f(broadcast1f(pA+2), a0, loadu4f(pX+2*ORD));
        storeu4f(pX, a0);
        vec4f b0 = mul4f(broadcast1f(pA+lda), b1);
        a1 = fnmadd4f(broadcast1f(pA+lda+1), b0, b2);
        a2 = fnmadd4f(broadcast1f(pA+lda+2), b0, loadu4f(pX+3*ORD));
        storeu4f(pX+ORD, b0);
        pA += 2*lda;
        pX += 2*ORD;
    }
#endif
    while ( pA < end ) // for ( int j = 0; j < N-2; ++j )
    {
        vec4f a0 = mul4f(broadcast1f(pA), a1);     // a1 = loadu4(pX);
        a1 = fnmadd4f(broadcast1f(pA+1), a0, a2);  // a2 = loadu4(pX+ORD);
        a2 = fnmadd4f(broadcast1f(pA+2), a0, loadu4f(pX+2*ORD));
        storeu4f(pX, a0);
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
/// specialized version for KD==2 and ORD==3
void alsatian_xtbsvLTN3(const int N, const float* pA, const int lda, float* pX)
{
    const float*const end = pA;
    constexpr int ORD = 3;
    const vec4f zero = setzero4f();
    pX += ( N - 1 ) * ORD;
    pA += ( N - 1 ) * lda;
    vec4f a1 = zero;
    if ( N >= 1 ) // j = N-1
    {
        vec4f a0 = loadu4f(pX);
        a1 = mul4f(broadcast1f(pA), a0);
        storeu4f(pX, blend4f(a1, a0, 0b1000));
        a1 = blend4f(a1, zero, 0b1000);
        pA -= lda;
        pX -= ORD;
    }
    vec4f a2 = a1;
    if ( N >= 2 ) // j = N-2
    {
        vec4f a0 = loadu4f(pX);
        a0 = fnmadd4f(broadcast1f(pA+1), a1, a0);
        a1 = mul4f(broadcast1f(pA), a0);
        storeu4f(pX, blend4f(a1, a0, 0b1000));
        a1 = blend4f(a1, zero, 0b1000);
        pA -= lda;
        pX -= ORD;
    }
    vec4f tt = broadcast1f(a1);
    vec4f af = loadu4f(pX);
    while ( pA >= end ) // for ( int j = N-3; j >= 0; --j )
    {
        vec4f a0 = fnmadd4f(broadcast1f(pA+2), a2, af);  // a2 = load3(pX+6);
        af = loadu4f(pX-ORD);
        a2 = a1;
        a0 = fnmadd4f(broadcast1f(pA+1), a1, a0);  // a1 = load3(pX+3);
        // restore 4th position that was saved in 'tt'
        a0 = blend4f(mul4f(broadcast1f(pA), a0), tt, 0b1000);
        tt = broadcast1f(a0); // save 4th position for next round
        storeu4f(pX, a0);
        a1 = blend4f(a0, zero, 0b1000);
        pA -= lda;
        pX -= ORD;
    }
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
        vec2 a0 = mul2(loaddup2(pA), a1);      // a1 = loadu4(pX);
        a1 = fnmadd2(loaddup2(pA+1), a0, a2);  // a2 = loadu4(pX+ORD);
        a2 = fnmadd2(loaddup2(pA+2), a0, load2(pX+2*ORD));
        store2(pX, a0);
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
    for ( int j = N-3; j >= 0; --j )
    {
        vec2 a0 = fnmadd2(loaddup2(pA+2), a2, load2(pX));
        a0 = fnmadd2(loaddup2(pA+1), a1, a0);
        a2 = a1;
        a1 = mul2(loaddup2(pA), a0);
        store2(pX, a1);
        pA -= lda;
        pX -= ORD;
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
        real a0 = pA[0] * a1;
        a1 = a2 - pA[1] * a0;
        a2 = pX[2] - pA[2] * a0;
        pX[0] = a0;
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
        real a0 = pX[0];
        a0 = a0 - pA[1] * a1;
        a1 = pA[0] * a0;
        pX[0] = a1;
        pA -= lda;
        pX -= 1;
    }
    for ( int j = N-3; j >= 0; --j )
    {
        real a0 = pX[0];
        a0 = a0 - pA[2] * a2;
        a0 = a0 - pA[1] * a1;
        a2 = a1;
        a1 = pA[0] * a0;
        pX[0] = a1;
        pA -= lda;
        pX -= 1;
    }
}

//------------------------------------------------------------------------------
#pragma mark - LAPACK-STYLE DPBTRS


inline void lapack_xpbtrs(char UPLO, int N, int KD, int NRHS, real const* AB, int LDAB, real* B, int LDB, int* INFO)
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
            blas_xtbsvLN<'N'>(N, KD, AB, LDAB, B+i*LDB, 1); //blas::xtbsv('L', 'N', 'N', N, KD, AB, LDAB, B+i*LDB, 1);
            blas_xtbsvLT<'N'>(N, KD, AB, LDAB, B+i*LDB, 1); //blas::xtbsv('L', 'T', 'N', N, KD, AB, LDAB, B+i*LDB, 1);
        }
    }
    else
        *INFO = 1;
}


template < int ORD >
inline void iso_xpbtrs(char UPLO, int N, int KD, real const* AB, int LDAB, real* B, int LDB, int* INFO)
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
inline void alsatian_xpbtrs(char UPLO, int N, int KD, real const* AB, int LDAB, real* B, int, int* INFO)
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
inline void alsatian_xpbtrsL(const int N, real const* AB, int LDAB, real* B)
{
#if defined(__AVX__) && REAL_IS_DOUBLE
    /* use routines for KD=2, and interleaved vectors of size `DIM*nbp` */
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
        ABORT_NOW("unexpected DIM!");
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

#endif
