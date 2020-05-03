// Cytosim was created by Francois Nedelec.
// Copyright 2020 Sainsbury Laboratory, Cambridge University

#ifndef XTRSM_H
#define XTRSM_H

/**
 This is a C-translation of the BLAS reference implementation of DTRSM
 FJN 03.05.2020
 
 SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
 
 SIDE   = 'L' : solves   op(A) * X = alpha * B.
        = 'R' : solves   X * op(A) = alpha * B.
 UPLO   = 'U' or 'L'
 TRANSA = 'N' or 'T' :  op(A) = A or op(A) = transpose(A)
 DIAG   = 'U' or 'N'   A is assumed to be unit triangular or not.
 M      = number of rows of B
 N      = number of columns of B.
 A      = DOUBLE PRECISION array of DIMENSION ( LDA, k )
 LDA    = leading dimension of A
 B      = DOUBLE PRECISION array of DIMENSION ( LDB, n )
 LDB    = leading dimension of B
 */


/**
 Solve A*X = alpha*B, overwriting B with X.
 DTRSM('L', 'L', 'N', 'N', N, 1, 1.0, A, LDA, tmp, N);
 
 DO J = 1,N
    IF (ALPHA.NE.ONE) THEN
        DO I = 1,M
            B(I,J) = ALPHA*B(I,J)
        CONTINUE
    END IF
    DO K = 1,M
        IF (B(K,J).NE.ZERO) THEN
            IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
            DO I = K + 1,M
                B(I,J) = B(I,J) - B(K,J)*A(I,K)
            CONTINUE
        END IF
    CONTINUE
CONTINUE
*/
void blas_xtrsmLLN(char Diag, int M, int N, real ALPHA, const real* A, int lda, real* B, int ldb)
{
    bool nounit = ( Diag == 'N' );
    for ( int J = 0; J < N; ++J )
    {
        if ( ALPHA != 1.0 )
        {
            for ( int I = 0; I < M; ++I )
                B[I+ldb*J] *= ALPHA;
        }
        for ( int K = 0; K < M; ++K )
        {
            if ( B[K+ldb*J] != 0.0 )
            {
                if (nounit)
                    B[K+ldb*J] /= A[K+lda*K];
                for ( int I = K + 1; I < M; ++I )
                    B[I+ldb*J] -= B[K+ldb*J] * A[I+lda*K];
            }
        }
    }
}


/**
 Solve transposed(A)*X = alpha*B, overwriting B with X.
 DTRSM('L', 'L', 'T', 'N', N, 1, 1.0, A, LDA, tmp, N);
 
 DO J = 1,N
     DO I = M,1,-1
         TEMP = ALPHA*B(I,J)
         DO K = I + 1,M
             TEMP = TEMP - A(K,I)*B(K,J)
         CONTINUE
         IF (NOUNIT) TEMP = TEMP/A(I,I)
         B(I,J) = TEMP
     CONTINUE
 CONTINUE
*/
void blas_xtrsmLLT(char Diag, int M, int N, real ALPHA, const real* A, int lda, real* B, int ldb)
{
    bool nounit = ( Diag == 'N' );
    for ( int J = 0; J < N; ++J )
    for ( int I = M-1; I >= 0; --I )
    {
        real temp = ALPHA * B[I+ldb*J];
        for ( int K = I + 1; K < M; ++K )
            temp -= A[K+lda*I] * B[K+ldb*J];
        if (nounit) temp /= A[I+lda*I];
        B[I+ldb*J] = temp;
    }
}


/**
 Solve A*X = alpha*B, overwriting B with X.
 DTRSM('L', 'U', 'N', 'N', N, 1, 1.0, A, LDA, tmp, N);
 
 DO J = 1,N
     IF (ALPHA.NE.ONE) THEN
         DO I = 1,M
             B(I,J) = ALPHA*B(I,J)
         CONTINUE
     END IF
     DO K = M,1,-1
         IF (B(K,J).NE.ZERO) THEN
             IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
             DO  I = 1,K - 1
                 B(I,J) = B(I,J) - B(K,J)*A(I,K)
             CONTINUE
         END IF
     CONTINUE
 CONTINUE
*/
void blas_xtrsmLUN(char Diag, int M, int N, real ALPHA, const real* A, int lda, real* B, int ldb)
{
    bool nounit = ( Diag == 'N' );
    for ( int J = 0; J < N; ++J )
    {
         if ( ALPHA != 1.0 )
         {
             for ( int I = 0; I < M; ++I )
                 B[I+ldb*J] *= ALPHA;
         }
        for ( int K = M-1; K >= 0; --K )
        {
             if ( B[K+ldb*J] != 0 )
             {
                 if (nounit)
                     B[K+ldb*J] /= A[K+lda*K];
                 for ( int I = 0; I < K; ++I )
                     B[I+ldb*J] -= B[K+ldb*J] * A[I+lda*K];
             }
        }
    }
}

/**
 Solve transposed(A)*X = alpha*B, overwriting B with X.
 DTRSM('L', 'U', 'T', Diag', N, 1, 1.0, A, LDA, tmp, N);

 DO J = 1,N
     DO I = 1,M
         TEMP = ALPHA*B(I,J)
         DO K = 1,I - 1
             TEMP = TEMP - A(K,I)*B(K,J)
         CONTINUE
         IF (NOUNIT) TEMP = TEMP/A(I,I)
         B(I,J) = TEMP
     CONTINUE
 CONTINUE
*/
void blas_xtrsmLUT(char Diag, int M, int N, real ALPHA, const real* A, int lda, real* B, int ldb)
{
    bool nounit = ( Diag == 'N' );
    for ( int J = 0; J < N; ++J )
    for ( int I = 0; I < M; ++I )
    {
        real temp = ALPHA * B[I+ldb*J];
        for ( int K = 0; K < I; ++K )
            temp -= A[K+lda*I] * B[K+ldb*J];
        if (nounit)
            temp /= A[I+lda*I];
        B[I+ldb*J] = temp;
    }
}

/**
 Patch to BLAS' DTRSM
 */
void blas_xtrsm(char Side, char Uplo, char Trans, char Diag, int M, int N, real ALPHA, const real* A, int LDA, real* B, int LDB)
{
    if ( Side == 'L' )
    {
        if ( Uplo == 'U' )
        {
            if ( Trans == 'N' )
                blas_xtrsmLUN(Diag, M, N, ALPHA, A, LDA, B, LDB);
            else
                blas_xtrsmLUT(Diag, M, N, ALPHA, A, LDA, B, LDB);
        }
        else if ( Uplo == 'L' )
        {
            if ( Trans == 'N' )
                blas_xtrsmLLN(Diag, M, N, ALPHA, A, LDA, B, LDB);
            else
                blas_xtrsmLLT(Diag, M, N, ALPHA, A, LDA, B, LDB);
        }
    }
    else
    {
        blas::xtrsm(Side, Uplo, Trans, Diag, M, N, ALPHA, A, LDA, B, LDB);
    }
}


//------------------------------------------------------------------------------
#pragma mark - ALSATIAN DTRSM


//------------------------------------------------------------------------------
#pragma mark - DIMENSION-SPECIFIC ALSATIAN DTRSM

#ifdef __AVX__

/// specialized version for ORD==3
void alsatian_xtrsmLLNN_3D(int M, int N, real ALPHA, const real* A, int lda, real* B, int ldb)
{
    constexpr int ORD = 3;
}


/// specialized version for ORD==3
void alsatian_xtrsmLLTN_3D(int M, int N, real ALPHA, const real* A, int lda, real* B, int ldb)
{
    constexpr int ORD = 3;
}

#endif


#ifdef __SSE3__

/// specialized version for ORD==2
void alsatian_xtrsmLLNN_2D(int M, int N, real ALPHA, const real* A, int lda, real* B, int ldb)
{
    constexpr int ORD = 2;
}


/// specialized version for ORD==2
void alsatian_xtrsmLLTN_2D(int M, int N, real ALPHA, const real* A, int lda, real* B, int ldb)
{
    constexpr int ORD = 2;
}

#endif


/// specialized version for ORD==1
void alsatian_xtrsmLLNN(int M, int N, real ALPHA, const real* A, int lda, real* B, int ldb)
{
}


/// specialized version for ORD==1
void alsatian_xtrsmLLTN(int M, int N, real ALPHA, const real* A, int lda, real* B, int ldb)
{
}

//------------------------------------------------------------------------------
#pragma mark - LAPACK-STYLE ROUTINES

inline void lapack_xgetrs(char TRANS, int N, int NRHS, const real* A, int LDA, const int* IPIV, real* B, int LDB, int& INFO)
{
    INFO = 0;
    if ( TRANS == 'N' )
    {
        // Apply row interchanges to the right hand sides.
        lapack::xlaswp(NRHS, B, LDB, 1, N, IPIV, 1);
        // Solve L*X = B, overwriting B with X.
        blas::xtrsm('L', 'L', 'N', 'U', N, NRHS, 1.0, A, LDA, B, LDB);
        // Solve U*X = B, overwriting B with X.
        blas::xtrsm('L', 'U', 'N', 'N', N, NRHS, 1.0, A, LDA, B, LDB);
    }
    else
    {
        // Solve U**T *X = B, overwriting B with X.
        blas::xtrsm('L', 'U', 'T', 'N', N, NRHS, 1.0, A, LDA, B, LDB );
        // Solve L**T *X = B, overwriting B with X.
        blas::xtrsm('L', 'L', 'T', 'U', N, NRHS, 1.0, A, LDA, B, LDB );
        // Apply row interchanges to the solution vectors.
        lapack::xlaswp(NRHS, B, LDB, 1, N, IPIV, -1);
    }
}


inline void lapack_xpotrs(char UPLO, int N, int NRHS, const real* A, int LDA, real* B, int LDB, int& INFO)
{
    INFO = 0;
    if ( UPLO == 'U' )
    {
        // Solve U**T *X = B, overwriting B with X.
        blas::xtrsm('L', 'U', 'T', 'N', N, NRHS, 1.0, A, LDA, B, LDB);
        // Solve U*X = B, overwriting B with X.
        blas::xtrsm('L', 'U', 'N', 'N', N, NRHS, 1.0, A, LDA, B, LDB);
    }
    else
    {
        // Solve L*X = B, overwriting B with X. ALPHA = 1.0
        blas::xtrsm('L', 'L', 'N', 'N', N, NRHS, 1.0, A, LDA, B, LDB);
        // Solve U*X = B, overwriting B with X. ALPHA = 1.0
        blas::xtrsm('L', 'L', 'T', 'N', N, NRHS, 1.0, A, LDA, B, LDB);
    }
}



inline void iso_xgetrsL(int N, const real* A, int LDA, const int* IPIV, real* B)
{
    /*
     we cannot call lapack::DGETRS('N', bks, 1, mec->block(), bks, mec->pivot(), Y, bks, &info);
     because the coordinates of the vector 'Y' are not contiguous but offset by 'DIM'.
     But calling DTBSV gets the required work done.
     */
    real * tmp = new_real(N);
    for ( int d = 0; d < DIM; ++d )
    {
        for ( int u = 0; u < N; ++u )
            tmp[u] = B[d+DIM*u];
#if 0
        int info = 0;
        lapack::xgetrs('N', N, 1, A, N, IPIV, tmp, N, &info);
#else
        // Apply row interchanges to the right hand sides.
        lapack::xlaswp(1, tmp, N, 1, N, IPIV, 1);
        // Solve L*X = B, overwriting B with X.
        blas::xtrsm('L', 'L', 'N', 'U', N, 1, 1.0, A, LDA, tmp, N);
        // Solve U*X = B, overwriting B with X.
        blas::xtrsm('L', 'U', 'N', 'N', N, 1, 1.0, A, LDA, tmp, N);
#endif
        for ( int u = 0; u < N; ++u )
            B[d+DIM*u] = tmp[u];
    }
    free_real(tmp);
}


inline void iso_xpotrsL(int N, const real* A, int LDA, real* B)
{
    /*
     we cannot call lapack::DPOTRS('L', N, 1, A, LDA, B, N, &info);
     because the coordinates of the vector 'Y' are not contiguous but offset by 'DIM'.
     But calling DTBSV gets the required work done.
     */
    real * tmp = new_real(N);
    for ( int d = 0; d < DIM; ++d )
    {
        for ( int u = 0; u < N; ++u )
            tmp[u] = B[d+DIM*u];
#if 1
        int info = 0;
        lapack::xpotrs('L', N, 1, A, LDA, tmp, N, &info);
#else
        // Solve L*X = B, overwriting B with X. ALPHA = 1.0
        blas::xtrsm('L', 'L', 'N', 'N', N, 1, 1.0, A, LDA, tmp, N);
        // Solve U*X = B, overwriting B with X. ALPHA = 1.0
        blas::xtrsm('L', 'L', 'T', 'N', N, 1, 1.0, A, LDA, tmp, N);
#endif
        for ( int u = 0; u < N; ++u )
            B[d+DIM*u] = tmp[u];
    }
    free_real(tmp);
}


#endif
