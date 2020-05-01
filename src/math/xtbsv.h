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
    if ( incX == 1 )
    {
        for (int j = N-1; j >= 0; --j)
        {
            if (X[j] != 0.)
            {
                const real * pA = A + K + j * lda;
                if ( nounit ) X[j] /= pA[0];
                real temp = X[j];
                const int inf = std::max(0, j-K);
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
                const real * pA = A + K + j * lda;
                if ( nounit ) X[jx] /= pA[0];
                real temp = X[jx];
                const int inf = std::max(0, j-K);
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

void blas_xtbsvLN(char Diag, int N, int K, const real* A, int lda, real* X, int incX)
{
    bool nounit = ( Diag == 'N' );
    if ( incX == 1 )
    {
        for (int j = 0; j < N; ++j)
        {
            if (X[j] != 0.)
            {
                const real * pA = A + j * lda;
                if ( nounit ) X[j] /= pA[0];
                real temp = X[j];
                const int sup = std::min(N-1, j+K);
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
                if ( nounit ) X[jx] /= pA[0];
                real temp = X[jx];
                const int sup = std::min(N-1, j+K);
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

void blas_xtbsvUT(char Diag, int N, int K, const real* A, int lda, real* X, int incX)
{
    bool nounit = ( Diag == 'N' );
    if ( incX == 1 )
    {
        for (int j = 0; j < N; ++j)
        {
            real temp = X[j];
            const real * pA = A + K + j * lda;
            for (int i = std::max(0, j-K); i < j; ++i)
                temp -= pA[i-j] * X[i];
            if ( nounit ) temp /= pA[0];
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
            const real * pA = A + K + j * lda;
            for (int i = std::max(0, j-K); i < j; ++i)
            {
                temp -= pA[i-j] * X[ix];
                ix += incX;
            }
            if ( nounit ) temp /= pA[0];
            X[jx] = temp;
            jx += incX;
            if (j >= K)
                kx += incX;
        }
    }
}

void blas_xtbsvLT(char Diag, int N, int K, const real* A, int lda, real* X, int incX)
{
    bool nounit = ( Diag == 'N' );
    if ( incX == 1 )
    {
        for (int j = N-1; j >= 0; --j)
        {
            real temp = X[j];
            const real * pA = A + j * lda;
            for (int i = std::min(N-1, j+K); i > j; --i)
                temp -= pA[i-j] * X[i];
            if ( nounit ) temp /= pA[0];
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
            for (int i = std::min(N-1, j+K); i > j; --i)
            {
                temp -= pA[i-j] * X[ix];
                ix -= incX;
            }
            if ( nounit ) temp /= pA[0];
            X[jx] = temp;
            jx -= incX;
            if ( j < N-K )
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
