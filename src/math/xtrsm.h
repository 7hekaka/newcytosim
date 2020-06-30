// Cytosim was created by Francois Nedelec.
// Copyright 2020 Sainsbury Laboratory, Cambridge University

#ifndef XTRSM_H
#define XTRSM_H

#ifdef __SSE3__
#  define XTRSM_USES_SSE3 REAL_IS_DOUBLE
#else
#  define XTRSM_USES_SSE3 0
#endif

#ifdef __AVX__
#  define XTRSM_USES_AVX REAL_IS_DOUBLE
#else
#  define XTRSM_USES_AVX 0
#endif

/**
 This is a C-translation of the BLAS reference implementation of DTRSM
 FJN 03.05.2020
 
 SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
 
 DTRSM  solves one of the matrix equations
 
    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
 
 where alpha is a scalar, X and B are m by n matrices, A is a unit, or
 non-unit, upper or lower triangular matrix and op(A) is one of
 
     op(A) = A   or   op(A) = A**T.
 
 The matrix X is overwritten on B.
 
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
template < char diag >
void blas_xtrsmLLN(const int M, const int N, real ALPHA, const real* A, const int lda, real* B, const int ldb)
{
    if ( ALPHA == 0.0 )
    {
        for ( int U = 0; U < N*M; ++U )
            B[U] = 0.0;
    }
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
                if ( diag == 'N' )
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
template < char diag >
void blas_xtrsmLLT(const int M, const int N, real ALPHA, const real* A, const int lda, real* B, const int ldb)
{
    if ( ALPHA == 0.0 )
    {
        for ( int U = 0; U < N*M; ++U )
            B[U] = 0.0;
    }
    for ( int J = 0; J < N; ++J )
    {
        for ( int I = M-1; I >= 0; --I )
        {
            real temp = ALPHA * B[I];
            for ( int K = I + 1; K < M; ++K )
                temp -= A[K+lda*I] * B[K];
            if ( diag == 'N' )
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
template < char diag >
void blas_xtrsmLUN(const int M, const int N, real ALPHA, const real* A, const int lda, real* B, const int ldb)
{
    if ( ALPHA == 0.0 )
    {
        for ( int U = 0; U < N*M; ++U )
            B[U] = 0.0;
    }
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
                 if ( diag == 'N' )
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
template < char diag >
void blas_xtrsmLUT(const int M, const int N, real ALPHA, const real* A, const int lda, real* B, const int ldb)
{
    if ( ALPHA == 0.0 )
    {
        for ( int U = 0; U < N*M; ++U )
            B[U] = 0.0;
    }
    for ( int J = 0; J < N; ++J )
    {
        for ( int I = 0; I < M; ++I )
        {
            real temp = ALPHA * B[I];
            for ( int K = 0; K < I; ++K )
                temp -= A[K+lda*I] * B[K];
            if ( diag == 'N' )
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
            if ( Trans == 'N' ) {
                if ( Diag == 'N' )
                    blas_xtrsmLUN<'N'>(M, N, ALPHA, A, LDA, B, LDB);
                else
                    blas_xtrsmLUN<'U'>(M, N, ALPHA, A, LDA, B, LDB);
            } else {
                if ( Diag == 'N' )
                    blas_xtrsmLUT<'N'>(M, N, ALPHA, A, LDA, B, LDB);
                else
                    blas_xtrsmLUT<'U'>(M, N, ALPHA, A, LDA, B, LDB);
            }
        }
        else if ( Uplo == 'L' )
        {
            if ( Trans == 'N' ) {
                if ( Diag == 'N' )
                    blas_xtrsmLLN<'N'>(M, N, ALPHA, A, LDA, B, LDB);
                else
                    blas_xtrsmLLN<'U'>(M, N, ALPHA, A, LDA, B, LDB);
            } else {
                if ( Diag == 'N' )
                    blas_xtrsmLLT<'N'>(M, N, ALPHA, A, LDA, B, LDB);
                else
                    blas_xtrsmLLT<'U'>(M, N, ALPHA, A, LDA, B, LDB);
            }
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
template < int ORD, char diag >
void iso_xtrsmLLN(const int M, const real* A, const int lda, real* B)
{
    for ( int K = 0; K < M; ++K )
    {
        real temp[ORD];
        if ( diag == 'N' ) {
            for ( int d = 0; d < ORD; ++d )
            {
                temp[d] = B[ORD*K+d] / A[K+lda*K];
                B[ORD*K+d] = temp[d];
            }
        } else if ( diag == 'I' ) {
            for ( int d = 0; d < ORD; ++d )
            {
                temp[d] = B[ORD*K+d] * A[K+lda*K];
                B[ORD*K+d] = temp[d];
            }
        } else {
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
     real temp = B[I];
     for ( int K = I + 1; K < M; ++K )
         temp -= A[K+lda*I] * B[K];
     if (nounit)
         temp /= A[I+lda*I];
     B[I] = temp;
 }
*/
template < int ORD, char diag >
void iso_xtrsmLLT(const int M, const real* A, const int lda, real* B)
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
        if ( diag == 'N' ) {
            for ( int d = 0; d < ORD; ++d )
                B[ORD*I+d] = temp[d] / A[I+lda*I];
        } else if ( diag == 'I' ) {
            for ( int d = 0; d < ORD; ++d )
                B[ORD*I+d] = temp[d] * A[I+lda*I];
        } else {
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
template < int ORD, char diag >
void iso_xtrsmLUN(const int M, const real* A, const int lda, real* B)
{
    for ( int K = M-1; K >= 0; --K )
    {
        if ( B[K] != 0 )
        {
            real temp[ORD];
            if ( diag == 'N' ) {
                for ( int d = 0; d < ORD; ++d )
                {
                    temp[d] = B[ORD*K+d] / A[K+lda*K];
                    B[ORD*K+d] = temp[d];
                }
            } else if ( diag == 'I' ) {
                for ( int d = 0; d < ORD; ++d )
                {
                    temp[d] = B[ORD*K+d] * A[K+lda*K];
                    B[ORD*K+d] = temp[d];
                }
            } else {
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

/**
 Solve transposed(A)*X = alpha*B, overwriting B with X.
 DTRSM('L', 'U', 'T', Diag, M, N, ALPHA, A, LDA, B, LDB);

 for ( int I = 0; I < M; ++I )
 {
     real temp = B[I];
     for ( int K = 0; K < I; ++K )
         temp -= A[K+lda*I] * B[K];
     if (nounit)
         temp /= A[I+lda*I];
     B[I] = temp;
 }
*/
template < int ORD, char diag >
void iso_xtrsmLUT(const int M, const real* A, const int lda, real* B)
{
    for ( int I = 0; I < M; ++I )
    {
        real temp[ORD];
        for ( int d = 0; d < ORD; ++d )
            temp[d] = B[I];
        for ( int K = 0; K < I; ++K )
        {
            for ( int d = 0; d < ORD; ++d )
                temp[d] -= A[K+lda*I] * B[ORD*K+d];
        }
        if ( diag == 'N' ) {
            for ( int d = 0; d < ORD; ++d )
                B[ORD*I+d] = temp[d] / A[I+lda*I];
        } else if ( diag == 'I' ) {
            for ( int d = 0; d < ORD; ++d )
                B[ORD*I+d] = temp[d] * A[I+lda*I];
        } else {
            for ( int d = 0; d < ORD; ++d )
                B[ORD*I+d] = temp[d];
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark - DIMENSION-SPECIFIC ALSATIAN DTRSM

#if XTRSM_USES_AVX

/// specialized version for ORD==3, FJN 4.5.2020
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
template < char diag >
void alsatian_xtrsmLLN3(const int M, const real* A, const int lda, real* B)
{
    for ( int K = 0; K < M; ++K )
    {
        vec4 T = loadu4(B+3*K);
        if ( diag == 'N' ) {
            vec4 n = T;
            T = div4(T, broadcast1(A+K));
            storeu4(B+3*K, blend4(T, n, 0b1000)); // blend to keep 4th value!
        } else if ( diag == 'I' ) {
            vec4 n = T;
            T = mul4(T, broadcast1(A+K)); // DIV
            storeu4(B+3*K, blend4(T, n, 0b1000)); // blend to keep 4th value!
        }
        int I = K + 1;
        real* pB = B+3*I;
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
#if 0
            // this unrolling may not work so well
            for ( ; I < M-7; I += 8 )
            {
                real const* pA = A + I;
#if 1
                vec4 a0 = broadcast1(pA  );
                vec4 a1 = broadcast1(pA+1);
                vec4 a2 = broadcast1(pA+2);
                vec4 a3 = broadcast1(pA+3);
                vec4 a4 = broadcast1(pA+4);
                vec4 a5 = broadcast1(pA+5);
                vec4 a6 = broadcast1(pA+6);
                vec4 a7 = broadcast1(pA+7);
#else
                vec4 a1 = broadcast2(pA);
                vec4 a3 = broadcast2(pA+2);
                vec4 a0 = unpacklo4(a1, a1);
                vec4 a2 = unpacklo4(a3, a3);
                a1 = unpackhi4(a1, a1);
                a3 = unpackhi4(a3, a3);
                vec4 a5 = broadcast2(pA+4);
                vec4 a7 = broadcast2(pA+6);
                vec4 a4 = unpacklo4(a5, a5);
                vec4 a6 = unpacklo4(a7, a7);
                a5 = unpackhi4(a5, a5);
                a7 = unpackhi4(a7, a7);
#endif
                vec4 t0 = fnmadd4(blend4(a0, a1, 0b1000), temp0, loadu4(pB  ));
                vec4 t1 = fnmadd4(blend4(a1, a2, 0b1100), temp1, loadu4(pB+4));
                vec4 t2 = fnmadd4(blend4(a2, a3, 0b1110), temp2, loadu4(pB+8));
                vec4 t3 = fnmadd4(blend4(a4, a5, 0b1000), temp0, loadu4(pB+12));
                vec4 t4 = fnmadd4(blend4(a5, a6, 0b1100), temp1, loadu4(pB+16));
                vec4 t5 = fnmadd4(blend4(a6, a7, 0b1110), temp2, loadu4(pB+20));
                storeu4(pB  , t0);
                storeu4(pB+4, t1);
                storeu4(pB+8, t2);
                storeu4(pB+12, t3);
                storeu4(pB+16, t4);
                storeu4(pB+20, t5);
                pB += 24;
            }
#endif
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
                vec4 t0 = fnmadd4(blend4(a0, a1, 0b1000), temp0, loadu4(pB  ));
                vec4 t1 = fnmadd4(blend4(a1, a2, 0b1100), temp1, loadu4(pB+4));
                vec4 t2 = fnmadd4(blend4(a2, a3, 0b1110), temp2, loadu4(pB+8));
                storeu4(pB  , t0);
                storeu4(pB+4, t1);
                storeu4(pB+8, t2);
                pB += 12;
            }
        }
        if ( I < M )
        {
            // load the next vector, before store4() will change it
            vec4 n = loadu4(pB);
            #pragma nounroll
            for ( ; I < M-1; ++I )
            {
                vec4 a = fnmadd4(T, broadcast1(A+I), n);
                n = loadu4(B+3*I+3);
                storeu4(B+3*I, a);
                pB += 3;
            }
            storeu4(pB, fnmadd4(T, broadcast1(A+I), n));
            assert_true(pB==B+3*I);
        }
        A += lda;
    }
}


/// specialized version for ORD==3, FJN 4.5.2020
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
template < char diag >
void alsatian_xtrsmLLT3(const int M, const real* A, const int lda, real* B)
{
    A += M * lda;
    for ( int I = M-1; I >= 0; --I )
    {
        A -= lda;
        real * pB = B+3*I;
        const vec4 ori = loadu4(pB); //(B+3*I);
        vec4 s0 = setzero4();  // temp
        vec4 s1 = setzero4();
        vec4 s2 = setzero4();
        pB += 3;
        int K = I + 1;
#if ( 0 )
        for ( ; K < M-7; K += 8 )
        {
            const real* pA = A + K;
            /*
             broadcast values of A:
             a0 = { AAAA } a1 = { BBBB } a2 = { CCCC } a3 = { DDDD }
             */
#if 1
            vec4 a0 = broadcast1(pA  );
            vec4 a1 = broadcast1(pA+1);
            vec4 a2 = broadcast1(pA+2);
            vec4 a3 = broadcast1(pA+3);
            vec4 a4 = broadcast1(pA+4);
            vec4 a5 = broadcast1(pA+5);
            vec4 a6 = broadcast1(pA+6);
            vec4 a7 = broadcast1(pA+7);
#else
            vec4 a1 = broadcast2(pA);
            vec4 a3 = broadcast2(pA+2);
            vec4 a0 = unpacklo4(a1, a1);
            vec4 a2 = unpacklo4(a3, a3);
            a1 = unpackhi4(a1, a1);
            a3 = unpackhi4(a3, a3);
            vec4 a5 = broadcast2(pA+4);
            vec4 a7 = broadcast2(pA+6);
            vec4 a4 = unpacklo4(a5, a5);
            vec4 a6 = unpacklo4(a7, a7);
            a5 = unpackhi4(a5, a5);
            a7 = unpackhi4(a7, a7);
#endif
            /*
             blend broadcasted values of A to generate the required vec4:
             { AAAB } { BBCC } { CDDD }
             */
            s0 = fnmadd4(blend4(a0, a1, 0b1000), loadu4(pB  ), s0);
            s1 = fnmadd4(blend4(a1, a2, 0b1100), loadu4(pB+4), s1);
            s2 = fnmadd4(blend4(a2, a3, 0b1110), loadu4(pB+8), s2);
            s0 = fnmadd4(blend4(a4, a5, 0b1000), loadu4(pB+12), s0);
            s1 = fnmadd4(blend4(a5, a6, 0b1100), loadu4(pB+16), s1);
            s2 = fnmadd4(blend4(a6, a7, 0b1110), loadu4(pB+20), s2);
            pB += 24;
        }
#endif
        for ( ; K < M-3; K += 4 )
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
        #pragma nounroll
        for ( ; K < M; ++K )
        {
            s0 = fnmadd4(broadcast1(A+K), loadu4(pB), s0);
            pB += 3;
        }
        if ( diag == 'N' )
            s0 = sub4(div4(s0, broadcast1(A+I)), ori);
        else if ( diag == 'I' )
            s0 = fmadd4(s0, broadcast1(A+I), ori);
        else
            s0 = add4(s0, ori);
        storeu4(B+3*I, blend4(s0, ori, 0b1000));
    }
}



/// specialized version for ORD==3, FJN 4.5.2020
/*
 A += M * lda;
 for ( int K = M-1; K >= 0; --K )
 {
     A -= lda;
     if (nounit)
         B[K] /= A[K];
     real temp = B[K];
     for ( int I = 0; I < K; ++I )
         B[I] -= temp * A[I];
 }
 */
template < char diag >
void alsatian_xtrsmLUN3(const int M, const real* A, const int lda, real* B)
{
    A += M * lda;
    for ( int K = M-1; K >= 0; --K )
    {
        A -= lda;
        vec4 T = loadu4(B+3*K);
        vec4 ori = T;
        if ( diag == 'N' ) {
            T = div4(ori, broadcast1(A+K));
            ori = blend4(T, ori, 0b1000);  // blend to keep 4th value!
        } else if ( diag == 'I' ) {
            T = mul4(ori, broadcast1(A+K)); // DIV
            ori = blend4(T, ori, 0b1000);
        }
        int I = 0;
        real * pB = B;
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
            for ( ; I+3 < K; I += 4 )
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
                vec4 t0 = fnmadd4(blend4(a0, a1, 0b1000), temp0, loadu4(pB  ));
                vec4 t1 = fnmadd4(blend4(a1, a2, 0b1100), temp1, loadu4(pB+4));
                vec4 t2 = fnmadd4(blend4(a2, a3, 0b1110), temp2, loadu4(pB+8));
                storeu4(pB  , t0);
                storeu4(pB+4, t1);
                storeu4(pB+8, t2);
                pB += 12;
            }
        }
        if ( I < K )
        {
            // load the next vector, before store4() will change it
            vec4 n = loadu4(B+3*I);
            #pragma nounroll
            for ( ; I < K-1; ++I )
            {
                vec4 a = fnmadd4(T, broadcast1(A+I), n);
                n = loadu4(B+3*I+3);
                storeu4(B+3*I, a);
            }
            // last is I = K-1
            storeu4(B+3*I, fnmadd4(T, broadcast1(A+I), n));
        }
        storeu4(B+3*K, ori);
    }
}



#endif


#if XTRSM_USES_SSE3

/// specialized version for ORD==2
template < char diag >
void alsatian_xtrsmLLN2(const int M, const real* A, const int lda, real* B)
{
    for ( int K = 0; K < M; ++K )
    {
        vec2 temp = load2(B+2*K);
        if ( diag == 'N' ) {
            temp = div2(temp, loaddup2(A+K));
            store2(B+2*K, temp);
        } else if ( diag == 'I' ) {
            temp = mul2(temp, loaddup2(A+K)); // DIV
            store2(B+2*K, temp);
        }
        // could unroll, using AVX
        for ( int I = K + 1; I < M; ++I )
            store2(B+2*I, fnmadd2(temp, loaddup2(A+I), load2(B+2*I)));
        A += lda;
    }
}


/// specialized version for ORD==2
template < char diag >
void alsatian_xtrsmLLT2(const int M, const real* A, const int lda, real* B)
{
    A += M * lda;
    for ( int I = M-1; I >= 0; --I )
    {
        A -= lda;
        vec2 temp = load2(B+2*I);
        // could unroll, using AVX
        for ( int K = I + 1; K < M; ++K )
            temp = fnmadd2(loaddup2(A+K), load2(B+2*K), temp);
        if ( diag == 'N' )
            temp = div2(temp, loaddup2(A+I));
        else if ( diag == 'I' )
            temp = mul2(temp, loaddup2(A+I)); // DIV
        store2(B+2*I, temp);
    }
}


/// specialized version for ORD==2
template < char diag >
void alsatian_xtrsmLUN2(const int M, const real* A, const int lda, real* B)
{
    A += M * lda;
    for ( int K = M-1; K >= 0; --K )
    {
        A -= lda;
        vec2 temp = load2(B+2*K);
        if ( diag == 'N' ) {
            temp = div2(temp, loaddup2(A+K));
            store2(B+2*K, temp);
        } else if ( diag == 'I' ) {
            temp = mul2(temp, loaddup2(A+K)); // DIV
            store2(B+2*K, temp);
        }
        for ( int I = 0; I < K; ++I )
            store2(B+2*I, fnmadd2(temp, loaddup2(A+I), load2(B+2*I)));
    }
}

#endif


/// specialized version for ORD==1
template < char diag >
void alsatian_xtrsmLLN1(const int M, const real* A, const int lda, real* B)
{
    for ( int K = 0; K < M; ++K )
    {
        real temp = B[K];
        if ( diag == 'N' ) {
            temp /= A[K];
            B[K] = temp;
        } else if ( diag == 'I' ) {
            temp *= A[K]; // DIV
            B[K] = temp;
        }
        # pragma ivdep
        for ( int I = K + 1; I < M; ++I )
            B[I] -= temp * A[I];
        A += lda;
    }
}


/// specialized version for ORD==1
template < char diag >
void alsatian_xtrsmLLT1(const int M, const real* A, const int lda, real* B)
{
    A += M * lda;
    for ( int I = M-1; I >= 0; --I )
    {
        A -= lda;
        real temp = B[I];
        # pragma ivdep
        for ( int K = I + 1; K < M; ++K )
            temp -= A[K] * B[K];
        if ( diag == 'N' )
            temp /= A[I];
        else if ( diag == 'I' )
            temp *= A[I]; // DIV
        B[I] = temp;
    }
}


/// specialized version for ORD==1
template < char diag >
void alsatian_xtrsmLUN1(const int M, const real* A, const int lda, real* B)
{
    A += M * lda;
    for ( int K = M-1; K >= 0; --K )
    {
        A -= lda;
        real temp = B[K];
        if ( diag == 'N' ) {
            temp /= A[K];
            B[K] = temp;
        } else if ( diag == 'I' ) {
            temp *= A[K]; // DIV
            B[K] = temp;
        }
        # pragma ivdep
        for ( int I = 0; I < K; ++I )
            B[I] -= temp * A[I];
    }
}

//------------------------------------------------------------------------------
#pragma mark - LAPACK-STYLE ROUTINES for Positive Symmetric Matrices

inline void lapack_xpotrs(char UPLO, int N, int NRHS, const real* A, int LDA, real* B, int LDB, int* INFO)
{
    *INFO = 0;
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
void iso_xpotrsL(const int N, const real* A, const int LDA, real* B)
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
    iso_xtrsmLLN<ORD, 'N'>(N, A, LDA, B);
    // Solve U*X = B, overwriting B with X. ALPHA = 1.0
    iso_xtrsmLLT<ORD, 'N'>(N, A, LDA, B);
#endif
}



/**
 This calls the standard lapack::xpotf2()
 and then inverts the diagonal terms
*/
void alsatian_xpotf2L(const int N, real* A, const int LDA, int* INFO)
{
    lapack::xpotf2('L', N, A, LDA, INFO);
    if ( 0 == *INFO )
    {
        const int S = LDA+1;
        for ( int u = 0; u < N; ++u )
            A[S*u] = 1.0 / A[S*u];
    }
}


template < int ORD >
void alsatian_xpotrsLtest(const int N, const real* A, const int LDA, real* B)
{
    real * tmp = new_real(N);
    for ( int d = 0; d < ORD; ++d )
    {
        for ( int u = 0; u < N; ++u )
            tmp[u] = B[d+ORD*u];
        // Solve L*X = B, overwriting B with X. ALPHA = 1.0
        alsatian_xtrsmLLN1<'I'>(N, A, LDA, tmp);
        // Solve U*X = B, overwriting B with X. ALPHA = 1.0
        alsatian_xtrsmLLT1<'I'>(N, A, LDA, tmp);
        for ( int u = 0; u < N; ++u )
            B[d+ORD*u] = tmp[u];
    }
    free_real(tmp);
}


template < int ORD >
void alsatian_xpotrsL(const int N, const real* A, const int LDA, real* B)
{
#if XTRSM_USES_AVX
    if ( ORD == 3 )
    {
        alsatian_xtrsmLLN3<'I'>(N, A, LDA, B);
        alsatian_xtrsmLLT3<'I'>(N, A, LDA, B);
    }
    else if ( ORD == 2 )
    {
        alsatian_xtrsmLLN2<'I'>(N, A, LDA, B);
        alsatian_xtrsmLLT2<'I'>(N, A, LDA, B);
    }
    else if ( ORD == 1 )
    {
        alsatian_xtrsmLLN1<'I'>(N, A, LDA, B);
        alsatian_xtrsmLLT1<'I'>(N, A, LDA, B);
    }
    else
        ABORT_NOW("unexpected DIM!");
#else
    iso_xtrsmLLN<ORD, 'I'>(N, A, LDA, B);
    iso_xtrsmLLT<ORD, 'I'>(N, A, LDA, B);
#endif
}

//------------------------------------------------------------------------------
#pragma mark - LAPACK-STYLE ROUTINES for General LU factorization

/*
 DLASWP performs a series of row interchanges on the matrix A.

  if ( INCX > 0 )
  {
      IX0 = K1
      I1 = K1
      I2 = K2
      INC = 1
  }
  else if ( INCX < 0 )
  {
      IX0 = K1 + ( K1-K2 )*INCX
      I1 = K2
      I2 = K1
      INC = -1
  }
  int IX = IX0;
  for ( int I = I1; I <= I2; I += INC )
  {
      int IP = IPIV[IX];
      if ( IP != I )
          swap<ORD>(A+ORD*I, A+ORD*IP);
      IX = IX + INCX;
  }
 */

template < int ORD >
inline void swap(real* A, real* B)
{
    for ( int d = 0; d < ORD; ++d )
    {
        real t = A[d];
        A[d] = B[d];
        B[d] = t;
    }
}


template < int ORD >
void iso_xlaswp(real* A, int K1, int K2, const int* IPIV, int INCX)
{
    if ( INCX == 1 )
    {
        for ( int I = K1-1; I < K2; ++I )
        {
            int IP = IPIV[I] - 1;  // zero-based array indexing!
            if ( IP != I )
                swap<ORD>(A+ORD*I, A+ORD*IP);
        }
    }
    else if ( INCX > 0 )
    {
        for ( int I = K1-1; I < K2; ++I )
        {
            int IP = IPIV[INCX*I] - 1;  // zero-based array indexing!
            if ( IP != I )
                swap<ORD>(A+ORD*I, A+ORD*IP);
        }
    }
    else if ( INCX < 0 )
    {
        int IX = K1 + ( K1-K2 )*INCX - 1;
        for ( int I = K2-1; I < K1; --I )
        {
            int IP = IPIV[IX] - 1;  // zero-based array indexing!
            if ( IP != I )
                swap<ORD>(A+ORD*I, A+ORD*IP);
            IX += INCX;
        }
    }
}


inline void lapack_xgetrs(char TRANS, int N, int NRHS, const real* A, int LDA, const int* IPIV, real* B, int LDB, int* INFO)
{
    *INFO = 0;
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


/// following the standard interface
void lapack_xgetrsN(int N, int NRHS, const real* A, int LDA, const int* IPIV, real* B, int LDB, int* INFO)
{
    *INFO = 0;
    for ( int i = 0; i < NRHS; ++i )
    {
        // Apply row interchanges to the right hand side.
        iso_xlaswp<1>(B, 1, N, IPIV, 1);
        // Solve L*X = B, overwriting B with X.
        alsatian_xtrsmLLN1<'U'>(N, A, LDA, B);
        // Solve U*X = B, overwriting B with X.
        alsatian_xtrsmLUN1<'N'>(N, A, LDA, B);
        B += LDB;
    }
}

/// version for just one RHS
void lapack_xgetrsN(int N, const real* A, int LDA, const int* IPIV, real* B)
{
    // Apply row interchanges to the right hand side.
    iso_xlaswp<1>(B, 1, N, IPIV, 1);
    // Solve L*X = B, overwriting B with X.
    alsatian_xtrsmLLN1<'U'>(N, A, LDA, B);
    // Solve U*X = B, overwriting B with X.
    alsatian_xtrsmLUN1<'N'>(N, A, LDA, B);
}


template < int ORD >
void iso_xgetrsN(const int N, const real* A, const int LDA, const int* IPIV, real* B)
{
    /*
     we cannot call lapack::DGETRS('N', bks, 1, mec->block(), bks, mec->pivot(), Y, bks, &info);
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
        //lapack::xgetrs('N', N, 1, A, N, IPIV, tmp, N, &info);
        // Apply row interchanges to the right hand side.
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
    // Apply row interchanges to the right hand side.
    iso_xlaswp<ORD>(B, 1, N, IPIV, 1);
    // Solve L*X = B, overwriting B with X.
    iso_xtrsmLLN<ORD, 'U'>(N, A, LDA, B);
    // Solve U*X = B, overwriting B with X.
    iso_xtrsmLUN<ORD, 'N'>(N, A, LDA, B);
#endif
}


template < int ORD >
void alsatian_xgetrsN(const int N, const real* A, const int LDA, const int* IPIV, real* B)
{
    // Apply row interchanges to the right hand side.
    iso_xlaswp<ORD>(B, 1, N, IPIV, 1);
#if XTRSM_USES_AVX
    if ( ORD == 3 )
    {
        // Solve L*X = B, overwriting B with X.
        alsatian_xtrsmLLN3<'U'>(N, A, LDA, B);
        // Solve U*X = B, overwriting B with X.
        alsatian_xtrsmLUN3<'N'>(N, A, LDA, B);
    }
    else if ( ORD == 2 )
    {
        // Solve L*X = B, overwriting B with X.
        alsatian_xtrsmLLN2<'U'>(N, A, LDA, B);
        // Solve U*X = B, overwriting B with X.
        alsatian_xtrsmLUN2<'N'>(N, A, LDA, B);
    }
    else if ( ORD == 1 )
    {
        // Solve L*X = B, overwriting B with X.
        alsatian_xtrsmLLN1<'U'>(N, A, LDA, B);
        // Solve U*X = B, overwriting B with X.
        alsatian_xtrsmLUN1<'N'>(N, A, LDA, B);
    }
    else
        ABORT_NOW("unexpected DIM!");
#else
    iso_xtrsmLLN<ORD, 'U'>(N, A, LDA, B);
    iso_xtrsmLLT<ORD, 'N'>(N, A, LDA, B);
#endif
}

#endif
