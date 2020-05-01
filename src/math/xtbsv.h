// Cytosim was created by Francois Nedelec.
// Copyright 2020 Sainsbury Laboratory, Cambridge University

#ifndef XTBSV_H
#define XTBSV_H

/**
 This is a C-translation of the BLAS reference implementation of DTBSV

 
 */
void blas_xtbsvUN(char Diag, int N, int K, const real* A, int lda, real* X, int incX)
{
    bool nounit = ( Diag == 'N' );
    /* pointer adjustments for indices starting from 1 */
    A -= 1 + lda;
    --X;
    int kplus1 = K + 1;
    
    if (incX == 1)
    {
        for (int j = N; j >= 1; --j)
        {
            if (X[j] != 0.)
            {
                int l = kplus1 - j;
                if (nounit)
                    X[j] /= A[kplus1 + j * lda];
                real temp = X[j];
                /* Computing MAX */
                int i__2 = 1, i__3 = j - K;
                int i__1 = std::max(i__2,i__3);
                for (int i__ = j - 1; i__ >= i__1; --i__)
                    X[i__] -= temp * A[l + i__ + j * lda];
            }
        }
    }
    else
    {
        int kx = 1;
        if ( incX > 0 )
            kx = 1 + (N - 1) * incX;
        int jx = kx;
        for (int j = N; j >= 1; --j)
        {
            kx -= incX;
            if (X[jx] != 0.)
            {
                int ix = kx;
                int l = kplus1 - j;
                if (nounit)
                    X[jx] /= A[kplus1 + j * lda];
                real temp = X[jx];
                /* Computing MAX */
                int i__2 = 1, i__3 = j - K;
                int i__1 = std::max(i__2,i__3);
                for (int i__ = j - 1; i__ >= i__1; --i__)
                {
                    X[ix] -= temp * A[l + i__ + j * lda];
                    ix -= incX;
                }
            }
            jx -= incX;
        }
    }
}

void blas_xtbsvLN(char Diag, int N, int K, const real* A, int lda, real* X, int incX)
{
    bool nounit = ( Diag == 'N' );
    /* pointer adjustments for indices starting from 1 */
    A -= 1 + lda;
    --X;

    if (incX == 1)
    {
        int i__1 = N;
        for (int j = 1; j <= i__1; ++j)
        {
            if (X[j] != 0.) {
                int l = 1 - j;
                if (nounit)
                    X[j] /= A[j * lda + 1];
                real temp = X[j];
                /* Computing MIN */
                int i__3 = N, i__4 = j + K;
                int i__2 = std::min(i__3,i__4);
                for (int i__ = j + 1; i__ <= i__2; ++i__)
                    X[i__] -= temp * A[l + i__ + j * lda];
            }
        }
    }
    else
    {
        int kx = 1;
        if ( incX <= 0 )
            kx = 1 - (N - 1) * incX;
        int jx = kx;
        int i__1 = N;
        for (int j = 1; j <= i__1; ++j)
        {
            kx += incX;
            if (X[jx] != 0.)
            {
                int ix = kx;
                int l = 1 - j;
                if (nounit)
                    X[jx] /= A[j * lda + 1];
                real temp = X[jx];
                /* Computing MIN */
                int i__3 = N, i__4 = j + K;
                int i__2 = std::min(i__3,i__4);
                for (int i__ = j + 1; i__ <= i__2; ++i__)
                {
                    X[ix] -= temp * A[l + i__ + j * lda];
                    ix += incX;
                }
            }
            jx += incX;
        }
    }
}

void blas_xtbsvUT(char Diag, int N, int K, const real* A, int lda, real* X, int incX)
{
    bool nounit = ( Diag == 'N' );
    /* pointer adjustments for indices starting from 1 */
    A -= 1 + lda;
    --X;
    int kplus1 = K + 1;
    
    if (incX == 1)
    {
        int i__1 = N;
        for (int j = 1; j <= i__1; ++j)
        {
            real temp = X[j];
            int l = kplus1 - j;
            /* Computing MAX */
            int i__2 = 1, i__3 = j - K;
            int i__4 = j - 1;
            for (int i__ = std::max(i__2,i__3); i__ <= i__4; ++i__)
                temp -= A[l + i__ + j * lda] * X[i__];
            if (nounit)
                temp /= A[kplus1 + j * lda];
            X[j] = temp;
        }
    }
    else
    {
        int kx = 1;
        if ( incX <= 0 )
            kx = 1 - (N - 1) * incX;
        int jx = kx;
        int i__1 = N;
        for (int j = 1; j <= i__1; ++j)
        {
            real temp = X[jx];
            int ix = kx;
            int l = kplus1 - j;
            /* Computing MAX */
            int i__4 = 1, i__2 = j - K;
            int i__3 = j - 1;
            for (int i__ = std::max(i__4,i__2); i__ <= i__3; ++i__)
            {
                temp -= A[l + i__ + j * lda] * X[ix];
                ix += incX;
            }
            if (nounit)
                temp /= A[kplus1 + j * lda];
            X[jx] = temp;
            jx += incX;
            if (j > K)
                kx += incX;
        }
    }
}

void blas_xtbsvLT(char Diag, int N, int K, const real* A, int lda, real* X, int incX)
{
    bool nounit = ( Diag == 'N' );
    /* pointer adjustments for indices starting from 1 */
    A -= 1 + lda;
    --X;
    if (incX == 1)
    {
        for (int j = N; j >= 1; --j)
        {
            real temp = X[j];
            int l = 1 - j;
            /* Computing MIN */
            int i__1 = N, i__3 = j + K;
            int i__4 = j + 1;
            for (int i__ = std::min(i__1,i__3); i__ >= i__4; --i__)
                temp -= A[l + i__ + j * lda] * X[i__];
            if (nounit)
                temp /= A[j * lda + 1];
            X[j] = temp;
        }
    }
    else
    {
        int kx = 1;
        if ( incX > 0 )
            kx = 1 + (N - 1) * incX;
        int jx = kx;
        for (int j = N; j >= 1; --j)
        {
            real temp = X[jx];
            int ix = kx;
            int l = 1 - j;
            /* Computing MIN */
            int i__4 = N, i__1 = j + K;
            int i__3 = j + 1;
            for (int i__ = std::min(i__4,i__1); i__ >= i__3; --i__)
            {
                temp -= A[l + i__ + j * lda] * X[ix];
                ix -= incX;
            }
            if (nounit)
                temp /= A[j * lda + 1];
            X[jx] = temp;
            jx -= incX;
            if (N - j >= K)
                kx -= incX;
        }
    }
}


void blas_xtbsv(char Uplo, char Trans, char Diag, int N, int K, const real* A, int lda, real*X, int incX)
{
    if ( Uplo == 'U' )
    {
        if ( Trans == 'N' )
            blas_xtbsvUN(Diag, N, K, A, lda, X, incX);
        else
            blas_xtbsvUT(Diag, N, K, A, lda, X, incX);
    }
    else if ( Uplo == 'L' )
    {
        if ( Trans == 'N' )
            blas_xtbsvLN(Diag, N, K, A, lda, X, incX);
        else
            blas_xtbsvLT(Diag, N, K, A, lda, X, incX);
    }
}

//------------------------------------------------------------------------------
#pragma mark - LAPACK DPBTF2


//------------------------------------------------------------------------------
#pragma mark - LAPACK DPBTRS


inline void lapack_xpbtrs(char UPLO, int N, int KD, int NRHS, real const* AB, int LDAB, real* B, int LDB, int* INFO)
{
    if ( UPLO == 'U' )
    {
        for ( int i = 0; i < NRHS; ++i )
        {
            blas::xtbsv('U', 'T', 'N', N, KD, AB, LDAB, B+i*LDB, 1);
            blas::xtbsv('U', 'N', 'N', N, KD, AB, LDAB, B+i*LDB, 1);
        }
    }
    else if ( UPLO == 'L' )
    {
        for ( int i = 0; i < NRHS; ++i )
        {
            blas::xtbsv('L', 'N', 'N', N, KD, AB, LDAB, B+i*LDB, 1);
            blas::xtbsv('L', 'T', 'N', N, KD, AB, LDAB, B+i*LDB, 1);
        }
    }
}

#endif
