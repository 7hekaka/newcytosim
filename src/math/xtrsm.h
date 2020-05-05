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
template < int ORD, bool nounit >
void iso_xtrsmLLN(int M, const real* A, int lda, real* B)
{
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
template < int ORD, bool nounit >
void iso_xtrsmLLT(int M, const real* A, int lda, real* B)
{
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
template < int ORD, bool nounit >
void iso_xtrsmLUN(int M, const real* A, int lda, real* B)
{
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
template < bool nounit >
void alsatian_xtrsmLLN_3D(int M, const real* A, int lda, real* B)
{
    for ( int K = 0; K < M; ++K )
    {
        vec4 T = loadu4(B+3*K);
        vec4 n = T;
        if ( nounit )
            T = mul4(n, broadcast1(A+K)); // DIV
        storeu4(B+3*K, blend4(T, n, 0b1000)); // blend to keep 4th value!
        int I = K + 1;
        {
            vec4 temp0, temp1, temp2;
            {
                /*
                 Convert temp = { XYZ? }
                 into temp0 = { XYZX } temp1 = { YZXY } temp2 = { ZXYZ }
                 */
                vec4 p = permute2f128(T, T, 0x01);
                vec4 h = shuffle4(T, p, 0b0001);
                temp0 = blend4(T, h, 0b1000);
                temp1 = blend4(h, p, 0b1100);
                temp2 = shuffle4(p, T, 0b0100);
            }
            // this loop could be unrolled further
            for ( ; I+3 < M; I += 4 )
            {
                /*
                 broadcast values of A:
                 a0 = { AAAA } a1 = { BBBB } a2 = { CCCC } a3 = { DDDD }
                 */
#if 1
                vec4 a0 = broadcast1(A+I  );
                vec4 a1 = broadcast1(A+I+1);
                vec4 a2 = broadcast1(A+I+2);
                vec4 a3 = broadcast1(A+I+3);
#else
                vec4 a1 = broadcast2(A+I);
                vec4 a3 = broadcast2(A+I+2);
                vec4 a0 = unpacklo4(a1, a1);
                vec4 a2 = unpacklo4(a3, a3);
                a1 = unpackhi4(a1, a1);
                a3 = unpackhi4(a3, a3);
#endif
                /*
                 blend broadcasted values of A to generate the required vec4:
                  { AAAB } { BBCC } { CDDD }
                 */
                real * pB = B+3*I;
                storeu4(pB  , fnmadd4(blend4(a0, a1, 0b1000), temp0, loadu4(pB  )));
                storeu4(pB+4, fnmadd4(blend4(a1, a2, 0b1100), temp1, loadu4(pB+4)));
                storeu4(pB+8, fnmadd4(blend4(a2, a3, 0b1110), temp2, loadu4(pB+8)));
            }
        }
        // load the next vector, before store4() will change it
        n = loadu4(B+3*I);
        for ( ; I < M; ++I )
        {
            vec4 a = fnmadd4(T, broadcast1(A+I), n);
            n = loadu4(B+3*I+3);
            storeu4(B+3*I, a);
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
     if ( nounit )
         temp *= A[I]; // DIV
     B[I] = temp;
}

 */
template < bool nounit >
void alsatian_xtrsmLLT_3D(int M, const real* A, int lda, real* B)
{
    A += M * lda;
    for ( int I = M-1; I >= 0; --I )
    {
        A -= lda;
        real * pB = B + 3 * I;
        const vec4 ori = loadu4(pB); //(B+3*I);
        vec4 s0 = setzero4();  // temp
        vec4 s1 = setzero4();
        vec4 s2 = setzero4();
        // can unroll
        pB += 3;
        int K = I + 1;
        for ( ; K+3 < M; K += 4 )
        {
            /*
             broadcast values of A:
             a0 = { AAAA } a1 = { BBBB } a2 = { CCCC } a3 = { DDDD }
             */
#if 1
            vec4 a0 = broadcast1(A+K  );
            vec4 a1 = broadcast1(A+K+1);
            vec4 a2 = broadcast1(A+K+2);
            vec4 a3 = broadcast1(A+K+3);
#else
            vec4 a1 = broadcast2(A+K);
            vec4 a3 = broadcast2(A+K+2);
            vec4 a0 = unpacklo4(a1, a1);
            vec4 a2 = unpacklo4(a3, a3);
            a1 = unpackhi4(a1, a1);
            a3 = unpackhi4(a3, a3);
#endif
            /*
             blend broadcasted values of A to generate the required vec4:
              { AAAB } { BBCC } { CDDD }
             */
            s0 = fnmadd4(blend4(a0, a1, 0b1000), loadu4(pB  ), s0); //(B+3*K  )
            s1 = fnmadd4(blend4(a1, a2, 0b1100), loadu4(pB+4), s1);
            s2 = fnmadd4(blend4(a2, a3, 0b1110), loadu4(pB+8), s2);
            pB += 12;
        }
        {
            /*
             Sum the X, Y and Z components:
             from s0 = { XYZX } s1 = { YZXY } s2 = { ZXYZ }
             into s0 = { X+X+X, Y+Y+Y, Z+Z+Z, ? }
             */
            vec4 h = shuffle4(blend4(s1, s0, 0b1000), s2, 0b0101);
            vec4 d3 = permute2f128(s1, s2, 0x21);
            vec4 d2 = shuffle4(s2, s1, 0b0101);
            vec4 d1 = permute2f128(h, h, 0x01);
            s0 = add4(add4(s0, d2), add4(d3, d1));
        }
        for ( ; K < M; ++K )
        {
            s0 = fnmadd4(broadcast1(A+K), loadu4(pB), s0);
            pB += 3;
        }
        if ( nounit )
            s0 = fmadd4(s0, broadcast1(A+I), ori);
        else
            s0 = add4(s0, ori);
        storeu4(B+3*I, blend4(s0, ori, 0b1000));
    }
}

#endif


#ifdef __SSE3__

/// specialized version for ORD==2
template < bool nounit >
void alsatian_xtrsmLLN_2D(int M, const real* A, int lda, real* B)
{
    for ( int K = 0; K < M; ++K )
    {
        vec2 temp = load2(B+2*K);
        if ( nounit )
            temp = mul2(temp, loaddup2(A+K)); // DIV
        store2(B+2*K, temp);
        for ( int I = K + 1; I < M; ++I )
            store2(B+2*I, fnmadd2(temp, loaddup2(A+I), load2(B+2*I)));
        A += lda;
    }
}


/// specialized version for ORD==2
template < bool nounit >
void alsatian_xtrsmLLT_2D(int M, const real* A, int lda, real* B)
{
    A += M * lda;
    for ( int I = M-1; I >= 0; --I )
    {
        A -= lda;
        vec2 temp = load2(B+2*I);
        // can unroll
        for ( int K = I + 1; K < M; ++K )
            temp = fnmadd2(loaddup2(A+K), load2(B+2*K), temp);
        if ( nounit )
            temp = mul2(temp, loaddup2(A+I)); // DIV
        store2(B+2*I, temp);
    }
}

#endif


/// specialized version for ORD==1
template < bool nounit >
void alsatian_xtrsmLLN_1D(int M, const real* A, int lda, real* B)
{
    for ( int K = 0; K < M; ++K )
    {
        real temp = B[K];
        if ( nounit )
            B[K] *= A[K]; // DIV
        B[K] = temp;
        for ( int I = K + 1; I < M; ++I )
            B[I] -= temp * A[I];
        A += lda;
    }
}


/// specialized version for ORD==1
template < bool nounit >
void alsatian_xtrsmLLT_1D(int M, const real* A, int lda, real* B)
{
    A += M * lda;
    for ( int I = M-1; I >= 0; --I )
    {
        A -= lda;
        real temp = B[I];
        for ( int K = I + 1; K < M; ++K )
            temp -= A[K] * B[K];
        if ( nounit )
            temp *= A[I]; // DIV
        B[I] = temp;
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
    iso_xtrsmLLN<ORD, 0>(N, A, LDA, B);
    // Solve U*X = B, overwriting B with X.
    iso_xtrsmLUN<ORD, 1>(N, A, LDA, B);
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
#if 0
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
    iso_xtrsmLLN<ORD, 1>(N, A, LDA, B);
    // Solve U*X = B, overwriting B with X. ALPHA = 1.0
    iso_xtrsmLLT<ORD, 1>(N, A, LDA, B);
#endif
}


#endif
