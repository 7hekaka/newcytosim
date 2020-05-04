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
 DTRSM('L', 'L', 'N', Diag, M, N, ALPHA, A, LDA, B, LDB);
 
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
    if ( ALPHA == 0.0 )
    {
        for ( int U = 0; U < N*M; ++U )
            B[U] = 0.0;
    }
    const bool nounit = ( Diag == 'N' );
    for ( int J = 0; J < N; ++J )
    {
        if ( ALPHA != 1.0 )
        {
            for ( int I = 0; I < M; ++I )
                B[I] *= ALPHA;
        }
        for ( int K = 0; K < M; ++K )
        {
            if ( B[K] != 0.0 )
            {
                if (nounit)
                    B[K] /= A[K+lda*K];
                real temp = B[K];
                for ( int I = K + 1; I < M; ++I )
                    B[I] -= temp * A[I+lda*K];
            }
        }
        B += ldb;
    }
}


/**
 Solve transposed(A)*X = alpha*B, overwriting B with X.
 DTRSM('L', 'L', 'T', Diag, M, N, ALPHA, A, LDA, B, LDB);
 
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
    if ( ALPHA == 0.0 )
    {
        for ( int U = 0; U < N*M; ++U )
            B[U] = 0.0;
    }
    const bool nounit = ( Diag == 'N' );
    for ( int J = 0; J < N; ++J )
    {
        for ( int I = M-1; I >= 0; --I )
        {
            real temp = ALPHA * B[I];
            for ( int K = I + 1; K < M; ++K )
                temp -= A[K+lda*I] * B[K];
            if (nounit)
                temp /= A[I+lda*I];
            B[I] = temp;
        }
        B += ldb;
    }
}


/**
 Solve A*X = alpha*B, overwriting B with X.
 DTRSM('L', 'U', 'N', Diag, M, N, ALPHA, A, LDA, B, LDB);
 
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
    if ( ALPHA == 0.0 )
    {
        for ( int U = 0; U < N*M; ++U )
            B[U] = 0.0;
    }
    const bool nounit = ( Diag == 'N' );
    for ( int J = 0; J < N; ++J )
    {
        if ( ALPHA != 1.0 )
        {
            for ( int I = 0; I < M; ++I )
                B[I] *= ALPHA;
        }
        for ( int K = M-1; K >= 0; --K )
        {
             if ( B[K] != 0 )
             {
                 if (nounit)
                     B[K] /= A[K+lda*K];
                 real temp = B[K];
                 for ( int I = 0; I < K; ++I )
                     B[I] -= temp * A[I+lda*K];
             }
        }
        B += ldb;
    }
}

/**
 Solve transposed(A)*X = alpha*B, overwriting B with X.
 DTRSM('L', 'U', 'T', Diag, M, N, ALPHA, A, LDA, B, LDB);

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
    if ( ALPHA == 0.0 )
    {
        for ( int U = 0; U < N*M; ++U )
            B[U] = 0.0;
    }
    const bool nounit = ( Diag == 'N' );
    for ( int J = 0; J < N; ++J )
    {
        for ( int I = 0; I < M; ++I )
        {
            real temp = ALPHA * B[I];
            for ( int K = 0; K < I; ++K )
                temp -= A[K+lda*I] * B[K];
            if (nounit)
                temp /= A[I+lda*I];
            B[I] = temp;
        }
        B += ldb;
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
#pragma mark - DTRSM-STYLE for interleaved vectors

/**
 Solve A*X = B, for 'ORD' interleaved vectors B
 like DTRSM('L', 'L', 'N', Diag, M, 1, 1.0, A, LDA, B, LDB);
 N = 1
 ALPHA = 1.0
 but for the ORD vectors
 
 for ( int K = 0; K < M; ++K )
 {
     if (nounit)
         B[K] /= A[K+lda*K];
     real temp = B[K];
     for ( int I = K + 1; I < M; ++I )
         B[I] -= temp * A[I+lda*K];
 }

 */
template < int ORD >
void iso_xtrsmLLN(char Diag, int M, const real* A, int lda, real* B)
{
    const bool nounit = ( Diag == 'N' );
    for ( int K = 0; K < M; ++K )
    {
        real temp[ORD];
        if (nounit)
        {
            for ( int d = 0; d < ORD; ++d )
            {
                temp[d] = B[ORD*K+d] / A[K+lda*K];
                B[ORD*K+d] = temp[d];
            }
        }
        else
        {
            for ( int d = 0; d < ORD; ++d )
                temp[d] = B[ORD*K+d];
        }
        for ( int I = K + 1; I < M; ++I )
        {
            for ( int d = 0; d < ORD; ++d )
                B[ORD*I+d] -= temp[d] * A[I+lda*K];
        }
    }
}


/**
 Solve transposed(A)*X = B, for 'ORD' interleaved vectors B
 DTRSM('L', 'L', 'T', Diag, M, 1, 1.0, A, LDA, B, LDB);
 N = 1
 ALPHA = 1.0
 
 for ( int I = M-1; I >= 0; --I )
 {
     real temp = ALPHA * B[I];
     for ( int K = I + 1; K < M; ++K )
         temp -= A[K+lda*I] * B[K];
     if (nounit)
         temp /= A[I+lda*I];
     B[I] = temp;
 }
*/
template < int ORD >
void iso_xtrsmLLT(char Diag, int M, const real* A, int lda, real* B)
{
    const bool nounit = ( Diag == 'N' );
    for ( int I = M-1; I >= 0; --I )
    {
        real temp[ORD];
        for ( int d = 0; d < ORD; ++d )
            temp[d] = B[ORD*I+d];
        for ( int K = I + 1; K < M; ++K )
        {
            for ( int d = 0; d < ORD; ++d )
                temp[d] -= A[K+lda*I] * B[ORD*K+d];
        }
        if (nounit)
        {
            for ( int d = 0; d < ORD; ++d )
                B[ORD*I+d] = temp[d] / A[I+lda*I];
        }
        else
        {
            for ( int d = 0; d < ORD; ++d )
                B[ORD*I+d] = temp[d];
        }
    }
}

/**
 Solve transposed(A)*X = B, for 'ORD' interleaved vectors B
 DTRSM('L', 'U', 'N', Diag, M, 1, 1.0, A, LDA, B, LDB);
 N = 1
 ALPHA = 1.0
 
 for ( int K = M-1; K >= 0; --K )
 {
     if (nounit)
         B[K] /= A[K+lda*K];
     real temp = B[K];
     for ( int I = 0; I < K; ++I )
         B[I] -= temp * A[I+lda*K];
 }
*/
template < int ORD >
void iso_xtrsmLUN(char Diag, int M, const real* A, int lda, real* B)
{
    const bool nounit = ( Diag == 'N' );
    for ( int K = M-1; K >= 0; --K )
    {
        if ( B[K] != 0 )
        {
            real temp[ORD];
            if (nounit)
            {
                for ( int d = 0; d < ORD; ++d )
                {
                    temp[d] = B[ORD*K+d] / A[K+lda*K];
                    B[ORD*K+d] = temp[d];
                }
            }
            else
            {
                for ( int d = 0; d < ORD; ++d )
                    temp[d] = B[ORD*K+d];
            }
            for ( int I = 0; I < K; ++I )
            {
                for ( int d = 0; d < ORD; ++d )
                    B[ORD*I+d] -= temp[d] * A[I+lda*K];
            }
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - DLASWP for interleaved vector

/*
 DLASWP performs a series of row interchanges on the matrix A.
 INCX = 1

 IF( INCX.GT.0 ) THEN
     IX0 = K1
     I1 = K1
     I2 = K2
     INC = 1
  ELSE IF( INCX.LT.0 ) THEN
     IX0 = K1 + ( K1-K2 )*INCX
     I1 = K2
     I2 = K1
     INC = -1
  ELSE
     RETURN
  END IF

   IX = IX0
   DO I = I1, I2, INC
      IP = IPIV( IX )
      IF( IP.NE.I ) THEN
         DO K = 1, N
            TEMP = A( I, K )
            A( I, K ) = A( IP, K )
            A( IP, K ) = TEMP
         CONTINUE
      END IF
      IX = IX + INCX
   CONTINUE
 */
template < int ORD >
void iso_xlaswp(real* A, int N, int K1, int K2, const int* IPIV, int INC)
{
    for ( int d = 0; d < ORD; ++d )
        lapack::xlaswp(1, A+d, N, K1, K2, IPIV, INC);
    ABORT_NOW("unfinished 'iso_xlaswp' code");
}

//------------------------------------------------------------------------------
#pragma mark - DIMENSION-SPECIFIC ALSATIAN DTRSM

/**
 This calls the standard lapack::xpotf2()
 and then inverts the diagonal terms
*/
void alsatian_xpotf2L(int N, real* A, int LDA, int* INFO)
{
    lapack::xpotf2('L', N, A, LDA, INFO);
    if ( 0 == *INFO )
    {
        const int S = LDA+1;
        for ( int u = 0; u < N; ++u )
            A[S*u] = 1.0 / A[S*u];
    }
}


#ifdef __AVX__

/// specialized version for ORD==3
/*
 for ( int K = 0; K < M; ++K )
 {
     real temp = B[K] * A[K]; // DIV
     B[K] = temp;
     for ( int I = K + 1; I < M; ++I )
         B[I] -= temp * A[I];
     A += lda;
 }
 */
void alsatian_xtrsmLLNN_3D(int M, const real* A, int lda, real* B)
{
    for ( int K = 0; K < M; ++K )
    {
        real * pB = B + 3 * K;
        vec4 temp0, temp1, temp2;
        const vec4 s = mul4(loadu4(pB), broadcast1(A+K)); // DIV
        store3(pB, s);
        {
            vec4 p = permute2f128(s, s, 0x01);
            vec4 h = shuffle4(s, p, 0b0001);
            temp0 = blend4(s, h, 0b1000);
            temp1 = blend4(h, p, 0b1100);
            temp2 = shuffle4(p, s, 0b0100);
        }
        pB += 3;
        int I = K + 1;
        for ( ; I+3 < M; I += 4 )
        {
            vec4 a0 = broadcast1(A+I  );
            vec4 a1 = broadcast1(A+I+1);
            vec4 a2 = broadcast1(A+I+2);
            vec4 a3 = broadcast1(A+I+3);
            storeu4(pB  , fnmadd4(blend4(a0, a1, 0b1000), temp0, loadu4(pB  )));
            storeu4(pB+4, fnmadd4(blend4(a1, a2, 0b1100), temp1, loadu4(pB+4)));
            storeu4(pB+8, fnmadd4(blend4(a2, a3, 0b1110), temp2, loadu4(pB+8)));
            pB += 12;
        }
        vec4 n = loadu4(pB);
        for ( ; I < M; ++I )
        {
            vec4 t = fnmadd4(s, broadcast1(A+I), n);
            n = loadu4(pB+3);
            storeu4(pB, t);
            pB += 3;
        }
        A += lda;
    }
}


/// specialized version for ORD==3
/*
 A += M * lda;
 for ( int I = M-1; I >= 0; --I )
 {
     A -= lda;
     real temp = B[I];
     for ( int K = I + 1; K < M; ++K )
         temp -= A[K] * B[K];
     B[I] = temp * A[I]; // DIV
 }

 */
void alsatian_xtrsmLLTN_3D(int M, const real* A, int lda, real* B)
{
    const vec4 zero = setzero4();
    A += M * lda;
    for ( int I = M-1; I >= 0; --I )
    {
        A -= lda;
#if 0
        vec4 temp = loadu4(B+3*I);
        for ( int K = I + 1; K < M; ++K )
            temp = fnmadd4(broadcast1(A+K), loadu4(B+3*K), temp);
        store3(B+3*I, mul4(temp, broadcast1(A+I)));
#else
        vec4 temp = loadu4(B+3*I);
        vec4 s0 = setzero4();
        vec4 s1 = setzero4();
        vec4 s2 = setzero4();
        // can unroll
        int K = I + 1;
        for ( ; K+3 < M; K += 4 )
        {
            vec4 a0 = broadcast1(A+K  );
            vec4 a1 = broadcast1(A+K+1);
            vec4 a2 = broadcast1(A+K+2);
            vec4 a3 = broadcast1(A+K+3);
            s0 = fnmadd4(blend4(a0, a1, 0b1000), loadu4(B+3*K  ), s0);
            s1 = fnmadd4(blend4(a1, a2, 0b1100), loadu4(B+3*K+4), s1);
            s2 = fnmadd4(blend4(a2, a3, 0b1110), loadu4(B+3*K+8), s2);
        }
        vec4 d = permute2f128(s0, s1, 0x21);
        d = shuffle4(d, s1, 0b0101);
        s0 = add4(s0, d);
        d = blend4(s1, s2, 0b0011);
        d = permute2f128(d, d, 0b0001);
        s0 = add4(s0, d);
        d = permute2f128(s2, s2, 0x01);
        d = shuffle4(s2, d, 0b0101);
        s0 = blend4(add4(s0, d), zero, 0b1000);
        temp = add4(temp, s0);
        for ( ; K < M; ++K )
            temp = fnmadd4(broadcast1(A+K), loadu4(B+3*K), temp);
        store3(B+3*I, mul4(temp, broadcast1(A+I)));
#endif
    }
}

#endif


#ifdef __SSE3__

/// specialized version for ORD==2
void alsatian_xtrsmLLNN_2D(int M, const real* A, int lda, real* B)
{
    for ( int K = 0; K < M; ++K )
    {
        vec2 temp = mul2(load2(B+2*K), loaddup2(A+K)); // DIV
        store2(B+2*K, temp);
        for ( int I = K + 1; I < M; ++I )
            store2(B+2*I, fnmadd2(temp, loaddup2(A+I), load2(B+2*I)));
        A += lda;
    }
}


/// specialized version for ORD==2
void alsatian_xtrsmLLTN_2D(int M, const real* A, int lda, real* B)
{
    A += M * lda;
    for ( int I = M-1; I >= 0; --I )
    {
        A -= lda;
        vec2 temp = load2(B+2*I);
        // can unroll
        for ( int K = I + 1; K < M; ++K )
            temp = fnmadd2(loaddup2(A+K), load2(B+2*K), temp);
        temp = mul2(temp, loaddup2(A+I)); // DIV
        store2(B+2*I, temp);
    }
}

#endif


/// specialized version for ORD==1
void alsatian_xtrsmLLNN(int M, const real* A, int lda, real* B)
{
    for ( int K = 0; K < M; ++K )
    {
        real temp = B[K] * A[K]; // DIV
        B[K] = temp;
        for ( int I = K + 1; I < M; ++I )
            B[I] -= temp * A[I];
        A += lda;
    }
}


/// specialized version for ORD==1
void alsatian_xtrsmLLTN(int M, const real* A, int lda, real* B)
{
    A += M * lda;
    for ( int I = M-1; I >= 0; --I )
    {
        A -= lda;
        real temp = B[I];
        for ( int K = I + 1; K < M; ++K )
            temp -= A[K] * B[K];
        B[I] = temp * A[I]; // DIV
    }
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
        blas::xtrsm('L', 'U', 'T', 'N', N, NRHS, 1.0, A, LDA, B, LDB);
        // Solve L**T *X = B, overwriting B with X.
        blas::xtrsm('L', 'L', 'T', 'U', N, NRHS, 1.0, A, LDA, B, LDB);
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


template < int ORD >
inline void iso_xgetrsL(int N, const real* A, int LDA, const int* IPIV, real* B)
{
    /*
     we cannot call lapack::DGETRS('N', bks, 1, mec->block(), bks, mec->pivot(), Y, bks, &info);
     because the coordinates of the vector 'Y' are not contiguous but offset by 'ORD'.
     But calling DTBSV gets the required work done.
     */
#if 1
    real * tmp = new_real(N);
    for ( int d = 0; d < ORD; ++d )
    {
        for ( int u = 0; u < N; ++u )
            tmp[u] = B[d+ORD*u];
        //int info = 0;
        //lapack::xgetrs('N', N, 1, A, N, IPIV, tmp, N, &info);
        // Apply row interchanges to the right hand sides.
        lapack::xlaswp(1, tmp, N, 1, N, IPIV, 1);
        // Solve L*X = B, overwriting B with X.
        blas::xtrsm('L', 'L', 'N', 'U', N, 1, 1.0, A, LDA, tmp, N);
        // Solve U*X = B, overwriting B with X.
        blas::xtrsm('L', 'U', 'N', 'N', N, 1, 1.0, A, LDA, tmp, N);
        for ( int u = 0; u < N; ++u )
            B[d+ORD*u] = tmp[u];
    }
    free_real(tmp);
#else
    // Apply row interchanges to the right hand sides.
    iso_xlaswp<ORD>(B, N, 1, N, IPIV, 1);
    // Solve L*X = B, overwriting B with X.
    iso_xtrsmLLN<ORD>('U', N, A, LDA, B);
    // Solve U*X = B, overwriting B with X.
    iso_xtrsmLUN<ORD>('N', N, A, LDA, B);
#endif
}


template < int ORD >
inline void iso_xpotrsL(int N, const real* A, int LDA, real* B)
{
    /*
     we cannot call lapack::DPOTRS('L', N, 1, A, LDA, B, N, &info);
     because the coordinates of the vector 'Y' are not contiguous but offset by 'ORD'.
     But calling DTBSV gets the required work done.
     */
#if 1
    real * tmp = new_real(N);
    for ( int d = 0; d < ORD; ++d )
    {
        for ( int u = 0; u < N; ++u )
            tmp[u] = B[d+ORD*u];
        //int info = 0;
        //lapack::xpotrs('L', N, 1, A, LDA, tmp, N, &info);
        // Solve L*X = B, overwriting B with X. ALPHA = 1.0
        blas::xtrsm('L', 'L', 'N', 'N', N, 1, 1.0, A, LDA, tmp, N);
        // Solve U*X = B, overwriting B with X. ALPHA = 1.0
        blas::xtrsm('L', 'L', 'T', 'N', N, 1, 1.0, A, LDA, tmp, N);
        for ( int u = 0; u < N; ++u )
            B[d+ORD*u] = tmp[u];
    }
    free_real(tmp);
#else
    // Solve L*X = B, overwriting B with X. ALPHA = 1.0
    iso_xtrsmLLN<ORD>('N', N, A, LDA, B);
    // Solve U*X = B, overwriting B with X. ALPHA = 1.0
    iso_xtrsmLLT<ORD>('N', N, A, LDA, B);
#endif
}


#endif
