// Cytosim was created by Francois Nedelec.
// Copyright 2020 Sainsbury Laboratory, Cambridge University

#ifndef XTRSM_H
#define XTRSM_H

#if USE_SIMD
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
 
    op(A)*X = alpha*B,   or   X*op(A) = alpha*B,
 
 where alpha is a scalar, X and B are M by N matrices, A is a unit, or
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
 NOUNIT = LSAME(DIAG,'N')
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
template < char diag, typename REAL >
void blas_xtrsmLLN(const int M, const int N, REAL ALPHA, const REAL* A, const int lda, REAL* B, const int ldb)
{
    if ( ALPHA == 0 )
    {
        for ( int U = 0; U < N*M; ++U )
            B[U] = 0;
        return;
    }
    for ( int J = 0; J < N; ++J )
    {
        if ( ALPHA != 1 )
        {
            for ( int I = 0; I < M; ++I )
                B[I] *= ALPHA;
        }
        for ( int K = 0; K < M; ++K )
        {
            if ( B[K] != 0 )
            {
                if ( diag == 'N' )
                    B[K] /= A[K+lda*K];
                else if ( diag == 'I' )
                    B[K] *= A[K+lda*K]; // DIV
                REAL tmp = B[K];
                for ( int I = K + 1; I < M; ++I )
                    B[I] -= tmp * A[I+lda*K];
            }
        }
        B += ldb;
    }
}


/**
 Solve transposed(A)*X = alpha*B, overwriting B with X.
 DTRSM('L', 'L', 'T', Diag, M, N, ALPHA, A, LDA, B, LDB);
 NOUNIT = LSAME(DIAG,'N')
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
template < char diag, typename REAL >
void blas_xtrsmLLT(const int M, const int N, REAL ALPHA, const REAL* A, const int lda, REAL* B, const int ldb)
{
    if ( ALPHA == 0 )
    {
        for ( int U = 0; U < N*M; ++U )
            B[U] = 0;
        return;
    }
    for ( int J = 0; J < N; ++J )
    {
        for ( int I = M-1; I >= 0; --I )
        {
            REAL tmp = ALPHA * B[I];
            for ( int K = I + 1; K < M; ++K )
                tmp -= A[K+lda*I] * B[K];
            if ( diag == 'N' )
                tmp /= A[I+lda*I];
            else if ( diag == 'I' )
                tmp *= A[I+lda*I]; // DIV
            B[I] = tmp;
        }
        B += ldb;
    }
}


/**
 Solve A*X = alpha*B, overwriting B with X.
 DTRSM('L', 'U', 'N', Diag, M, N, ALPHA, A, LDA, B, LDB);
 NOUNIT = LSAME(DIAG,'N')
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
template < char diag, typename REAL >
void blas_xtrsmLUN(const int M, const int N, REAL ALPHA, const REAL* A, const int lda, REAL* B, const int ldb)
{
    if ( ALPHA == 0 )
    {
        for ( int U = 0; U < N*M; ++U )
            B[U] = 0;
        return;
    }
    for ( int J = 0; J < N; ++J )
    {
        if ( ALPHA != 1 )
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
                 else if ( diag == 'I' )
                     B[K] *= A[K+lda*K]; // DIV
                 REAL tmp = B[K];
                 for ( int I = 0; I < K; ++I )
                     B[I] -= tmp * A[I+lda*K];
             }
        }
        B += ldb;
    }
}

/**
 Solve transposed(A)*X = alpha*B, overwriting B with X.
 DTRSM('L', 'U', 'T', Diag, M, N, ALPHA, A, LDA, B, LDB);
 NOUNIT = LSAME(DIAG,'N')
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
template < char diag, typename REAL >
void blas_xtrsmLUT(const int M, const int N, REAL ALPHA, const REAL* A, const int lda, REAL* B, const int ldb)
{
    if ( ALPHA == 0 )
    {
        for ( int U = 0; U < N*M; ++U )
            B[U] = 0;
        return;
    }
    for ( int J = 0; J < N; ++J )
    {
        for ( int I = 0; I < M; ++I )
        {
            REAL tmp = ALPHA * B[I];
            for ( int K = 0; K < I; ++K )
                tmp -= A[K+lda*I] * B[K];
            if ( diag == 'N' )
                tmp /= A[I+lda*I];
            else if ( diag == 'I' )
                tmp *= A[I+lda*I]; // DIV
            B[I] = tmp;
        }
        B += ldb;
    }
}

/**
 Patch to a BLAS' style interface to DTRSM
 */
template < typename REAL >
void blas_xtrsm(char Side, char Uplo, char Trans, char Diag, int M, int N, REAL ALPHA, const REAL* A, int LDA, REAL* B, int LDB)
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
 Solve A*X = B, when B contains 'ORD' interleaved vectors of size M
 A is a unit, or non-unit, lower triangular matrix of size M x M

 like DTRSM('L', 'L', 'N', Diag, M, 1, 1.0, A, LDA, B, LDB);
 with N = 1 and ALPHA = 1.0

 for ( int K = 0; K < M; ++K )
 {
     if (nounit)
         B[K] /= A[K+lda*K];
     real tmp = B[K];
     for ( int I = K + 1; I < M; ++I )
         B[I] -= tmp * A[I+lda*K];
 }

 */
template < int ORD, char diag, typename FLOAT, typename REAL >
void iso_xtrsmLLN(const int M, const FLOAT* A, const int lda, REAL* B)
{
    assert_true( M <= lda );
    for ( int K = 0; K < M; ++K )
    {
        REAL tmp[ORD];
        if ( diag == 'N' ) {
            for ( int d = 0; d < ORD; ++d )
            {
                tmp[d] = B[ORD*K+d] / A[K+lda*K];
                B[ORD*K+d] = tmp[d];
            }
        } else if ( diag == 'I' ) {
            for ( int d = 0; d < ORD; ++d )
            {
                tmp[d] = B[ORD*K+d] * A[K+lda*K];
                B[ORD*K+d] = tmp[d];
            }
        } else {
            for ( int d = 0; d < ORD; ++d )
                tmp[d] = B[ORD*K+d];
        }
        for ( int I = K + 1; I < M; ++I )
        {
            for ( int d = 0; d < ORD; ++d )
                B[ORD*I+d] -= tmp[d] * A[I+lda*K];
        }
    }
}


/**
 Solve transposed(A)*X = B, when B contains 'ORD' interleaved vectors of size M
 A is a unit, or non-unit, lower triangular matrix of size M x M
 
 like DTRSM('L', 'L', 'T', Diag, M, 1, 1.0, A, LDA, B, LDB);
 with N = 1 and ALPHA = 1.0
 
 for ( int I = M-1; I >= 0; --I )
 {
     real tmp = B[I];
     for ( int K = I + 1; K < M; ++K )
         tmp -= A[K+lda*I] * B[K];
     if (nounit)
         tmp /= A[I+lda*I];
     B[I] = tmp;
 }
*/
template < int ORD, char diag, typename FLOAT, typename REAL >
void iso_xtrsmLLT(const int M, const FLOAT* A, const int lda, REAL* B)
{
    assert_true( M <= lda );
    for ( int I = M-1; I >= 0; --I )
    {
        REAL tmp[ORD];
        for ( int d = 0; d < ORD; ++d )
            tmp[d] = B[ORD*I+d];
        for ( int K = I + 1; K < M; ++K )
        {
            for ( int d = 0; d < ORD; ++d )
                tmp[d] -= A[K+lda*I] * B[ORD*K+d];
        }
        if ( diag == 'N' ) {
            for ( int d = 0; d < ORD; ++d )
                B[ORD*I+d] = tmp[d] / A[I+lda*I];
        } else if ( diag == 'I' ) {
            for ( int d = 0; d < ORD; ++d )
                B[ORD*I+d] = tmp[d] * A[I+lda*I];
        } else {
            for ( int d = 0; d < ORD; ++d )
                B[ORD*I+d] = tmp[d];
        }
    }
}

/**
 Solve A*X = B, when B contains 'ORD' interleaved vectors of size M
 A is a unit, or non-unit, upper triangular matrix of size M x M
 
 like DTRSM('L', 'U', 'N', Diag, M, 1, 1.0, A, LDA, B, LDB);
 with N = 1 and ALPHA = 1.0

 for ( int K = M-1; K >= 0; --K )
 {
     if (nounit)
         B[K] /= A[K+lda*K];
     real tmp = B[K];
     for ( int I = 0; I < K; ++I )
         B[I] -= tmp * A[I+lda*K];
 }
*/
template < int ORD, char diag, typename FLOAT, typename REAL >
void iso_xtrsmLUN(const int M, const FLOAT* A, const int lda, REAL* B)
{
    assert_true( M <= lda );
    for ( int K = M-1; K >= 0; --K )
    {
        if ( B[K] != 0 )
        {
            REAL tmp[ORD];
            if ( diag == 'N' ) {
                for ( int d = 0; d < ORD; ++d )
                {
                    tmp[d] = B[ORD*K+d] / A[K+lda*K];
                    B[ORD*K+d] = tmp[d];
                }
            } else if ( diag == 'I' ) {
                for ( int d = 0; d < ORD; ++d )
                {
                    tmp[d] = B[ORD*K+d] * A[K+lda*K];
                    B[ORD*K+d] = tmp[d];
                }
            } else if ( diag == 'C' ) {
                // multiply vector component by diagonal term:
                for ( int d = 0; d < ORD; ++d )
                {
                    tmp[d] = B[ORD*K+d];
                    B[ORD*K+d] = tmp[d] * A[K+lda*K];
                }
            } else if ( diag == 'U' ) {
                for ( int d = 0; d < ORD; ++d )
                    tmp[d] = B[ORD*K+d];
            } else {
                printf("Error: unknown `diag` in alsatian_xtrsmLUN<>\n");
            }

            for ( int I = 0; I < K; ++I )
            {
                for ( int d = 0; d < ORD; ++d )
                    B[ORD*I+d] -= tmp[d] * A[I+lda*K];
            }
        }
    }
}

/**
 Solve transposed(A)*X = B, when B contains 'ORD' interleaved vectors of size M
 A is a unit, or non-unit, upper triangular matrix of size M x M

 like DTRSM('L', 'U', 'T', Diag, M, N, ALPHA, A, LDA, B, LDB);
 with N = 1 and ALPHA = 1.0

 for ( int I = 0; I < M; ++I )
 {
     real tmp = B[I];
     for ( int K = 0; K < I; ++K )
         tmp -= A[K+lda*I] * B[K];
     if (nounit)
         tmp /= A[I+lda*I];
     B[I] = tmp;
 }
*/
template < int ORD, char diag, typename FLOAT, typename REAL >
void iso_xtrsmLUT(const int M, const FLOAT* A, const int lda, REAL* B)
{
    assert_true( M <= lda );
    for ( int I = 0; I < M; ++I )
    {
        REAL tmp[ORD];
        for ( int d = 0; d < ORD; ++d )
            tmp[d] = B[I];
        for ( int K = 0; K < I; ++K )
        {
            for ( int d = 0; d < ORD; ++d )
                tmp[d] -= A[K+lda*I] * B[ORD*K+d];
        }
        if ( diag == 'N' ) {
            for ( int d = 0; d < ORD; ++d )
                B[ORD*I+d] = tmp[d] / A[I+lda*I];
        } else if ( diag == 'I' ) {
            for ( int d = 0; d < ORD; ++d )
                B[ORD*I+d] = tmp[d] * A[I+lda*I];
        } else {
            for ( int d = 0; d < ORD; ++d )
                B[ORD*I+d] = tmp[d];
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark - 3D AVX ALSATIAN DTRSM

#if defined(__AVX__)

/// specialized version for ORD==3, FJN 4.5.2020
/*
 for ( int K = 0; K < M; ++K )
 {
     real tmp = B[K] * A[K]; // DIV
     B[K] = tmp;
     for ( int I = K + 1; I < M; ++I )
         B[I] -= tmp * A[I];
     A += lda;
 }
 */
template < char diag >
void alsatian_iso3_xtrsmLLN_AVX(const int M, const double* A, const int lda, double* B)
{
    assert_true( M <= lda );
    for ( int K = 0; K < M; ++K )
    {
        vec4 T = loadu4(B+3*K);
        if ( diag == 'N' ) {
            vec4 n = T;
            T = div4(T, broadcast1(A+K));
            storeu4(B+3*K, blend31(T, n)); // blend to keep 4th value!
        } else if ( diag == 'I' ) {
            vec4 n = T;
            T = mul4(T, broadcast1(A+K)); // DIV
            storeu4(B+3*K, blend31(T, n)); // blend to keep 4th value!
        }
        int I = K + 1;
        double* pB = B+3*I;
        {
            vec4 temp0, temp1, temp2;
            {
                /*
                 Convert temp = { XYZ? }
                 into temp0 = { XYZX } temp1 = { YZXY } temp2 = { ZXYZ }
                 */
                vec4 p = swap2f128(T);
                vec4 h = shuffle4(T, p, 0b0001);
                temp0 = blend31(T, h);
                temp1 = blend22(h, p);
                temp2 = shuffle4(p, T, 0b0100);
            }
#if 0
            // this unrolling may not work so well
            for ( ; I < M-7; I += 8 )
            {
                double const* pA = A + I;
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
                vec4 t0 = fnmadd4(blend31(a0, a1), temp0, loadu4(pB  ));
                vec4 t1 = fnmadd4(blend22(a1, a2), temp1, loadu4(pB+4));
                vec4 t2 = fnmadd4(blend13(a2, a3), temp2, loadu4(pB+8));
                vec4 t3 = fnmadd4(blend31(a4, a5), temp0, loadu4(pB+12));
                vec4 t4 = fnmadd4(blend22(a5, a6), temp1, loadu4(pB+16));
                vec4 t5 = fnmadd4(blend13(a6, a7), temp2, loadu4(pB+20));
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
                vec4 t0 = fnmadd4(blend31(a0, a1), temp0, loadu4(pB  ));
                vec4 t1 = fnmadd4(blend22(a1, a2), temp1, loadu4(pB+4));
                vec4 t2 = fnmadd4(blend13(a2, a3), temp2, loadu4(pB+8));
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
     real tmp = B[I];
     for ( int K = I + 1; K < M; ++K )
         tmp -= A[K] * B[K];
     if ( nounit )
         tmp *= A[I]; // DIV
     B[I] = tmp;
}

 */
template < char diag >
void alsatian_iso3_xtrsmLLT_AVX(const int M, const double* A, const int lda, double* B)
{
    assert_true( M <= lda );
    A += M * lda;
    for ( int I = M-1; I >= 0; --I )
    {
        A -= lda;
        double * pB = B+3*I;
        const vec4 ori = loadu4(pB); //(B+3*I);
        vec4 s0 = setzero4();  // temp
        vec4 s1 = setzero4();
        vec4 s2 = setzero4();
        pB += 3;
        int K = I + 1;
#if ( 0 )
        for ( ; K < M-7; K += 8 )
        {
            const double* pA = A + K;
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
            s0 = fnmadd4(blend31(a0, a1), loadu4(pB  ), s0);
            s1 = fnmadd4(blend22(a1, a2), loadu4(pB+4), s1);
            s2 = fnmadd4(blend13(a2, a3), loadu4(pB+8), s2);
            s0 = fnmadd4(blend31(a4, a5), loadu4(pB+12), s0);
            s1 = fnmadd4(blend22(a5, a6), loadu4(pB+16), s1);
            s2 = fnmadd4(blend13(a6, a7), loadu4(pB+20), s2);
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
            s0 = fnmadd4(blend31(a0, a1), loadu4(pB  ), s0); //(B+3*K  )
            s1 = fnmadd4(blend22(a1, a2), loadu4(pB+4), s1);
            s2 = fnmadd4(blend13(a2, a3), loadu4(pB+8), s2);
            pB += 12;
        }
        {
            /*
             Sum the X, Y and Z components:
             from s0 = { XYZX } s1 = { YZXY } s2 = { ZXYZ }
             into s0 = { X+X+X, Y+Y+Y, Z+Z+Z, ? }
             */
            vec4 h = shuffle4(blend31(s1, s0), s2, 0b0101);
            vec4 d3 = catshift2(s1, s2);
            vec4 d2 = shuffle4(s2, s1, 0b0101);
            vec4 d1 = swap2f128(h);
            s0 = add4(add4(s0, d2), add4(d3, d1));
        }
        #pragma nounroll
        for ( ; K < M; ++K )
        {
            s0 = fnmadd4(broadcast1(A+K), loadu4(pB), s0);
            pB += 3;
        }
        s0 = add4(s0, ori);
        if ( diag == 'N' )
            s0 = div4(s0, broadcast1(A+I));
        else if ( diag == 'I' )
            s0 = mul4(s0, broadcast1(A+I));
        storeu4(B+3*I, blend31(s0, ori));
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
     real tmp = B[K];
     for ( int I = 0; I < K; ++I )
         B[I] -= tmp * A[I];
 }
 */
template < char diag >
void alsatian_iso3_xtrsmLUN_AVX(const int M, const double* A, const int lda, double* B)
{
    assert_true( M <= lda );
    A += M * lda;
    for ( int K = M-1; K >= 0; --K )
    {
        A -= lda;
        vec4 T = loadu4(B+3*K);
        vec4 ori = T;
        if ( diag == 'N' ) {
            T = div4(T, broadcast1(A+K));
            ori = blend31(T, ori);  // blend to keep 4th value!
        } else if ( diag == 'I' ) {
            T = mul4(T, broadcast1(A+K)); // DIV
            ori = blend31(T, ori);
        } else if ( diag == 'C' ) {
            // multiply vector component by diagonal term:
            ori = blend31(mul4(T, broadcast1(A+K)), T);
        } else if ( diag != 'U' ) {
            printf("Error: unknown `diag` in alsatian_xtrsmLUN3\n");
        }
        int I = 0;
        double * pB = B;
        {
            vec4 temp0, temp1, temp2;
            {
                /*
                 Convert temp = { XYZ? }
                 into temp0 = { XYZX } temp1 = { YZXY } temp2 = { ZXYZ }
                 */
                vec4 p = swap2f128(T);
                vec4 h = shuffle4(T, p, 0b0001);
                temp0 = blend31(T, h);
                temp1 = blend22(h, p);
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
                vec4 t0 = fnmadd4(blend31(a0, a1), temp0, loadu4(pB  ));
                vec4 t1 = fnmadd4(blend22(a1, a2), temp1, loadu4(pB+4));
                vec4 t2 = fnmadd4(blend13(a2, a3), temp2, loadu4(pB+8));
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

//------------------------------------------------------------------------------
#pragma mark - 2D SIMD ALSATIAN DTRSM


#if XTRSM_USES_SSE3

/// specialized version for ORD==2
template < char diag >
void alsatian_iso2_xtrsmLLN_SIMD(const int M, const double* A, const int lda, double* B)
{
    assert_true( M <= lda );
    for ( int K = 0; K < M; ++K )
    {
        vec2 tmp = load2(B+2*K);
        if ( diag == 'N' ) {
            tmp = div2(tmp, loaddup2(A+K));
            store2(B+2*K, tmp);
        } else if ( diag == 'I' ) {
            tmp = mul2(tmp, loaddup2(A+K)); // DIV
            store2(B+2*K, tmp);
        }
        // could unroll, using AVX
        for ( int I = K + 1; I < M; ++I )
            store2(B+2*I, fnmadd2(tmp, loaddup2(A+I), load2(B+2*I)));
        A += lda;
    }
}


/// specialized version for ORD==2
template < char diag >
void alsatian_iso2_xtrsmLLT_SIMD(const int M, const double* A, const int lda, double* B)
{
    assert_true( M <= lda );
    A += M * lda;
    for ( int I = M-1; I >= 0; --I )
    {
        A -= lda;
        vec2 tmp = load2(B+2*I);
        // could unroll, using AVX
        for ( int K = I + 1; K < M; ++K )
            tmp = fnmadd2(loaddup2(A+K), load2(B+2*K), tmp);
        if ( diag == 'N' )
            tmp = div2(tmp, loaddup2(A+I));
        else if ( diag == 'I' )
            tmp = mul2(tmp, loaddup2(A+I)); // DIV
        store2(B+2*I, tmp);
    }
}


/// specialized version for ORD==2
template < char diag >
void alsatian_iso2_xtrsmLUN_SIMD(const int M, const double* A, const int lda, double* B)
{
    assert_true( M <= lda );
    A += M * lda;
    for ( int K = M-1; K >= 0; --K )
    {
        A -= lda;
        vec2 tmp = load2(B+2*K);
        if ( diag == 'N' ) {
            tmp = div2(tmp, loaddup2(A+K));
            store2(B+2*K, tmp);
        } else if ( diag == 'I' ) {
            tmp = mul2(tmp, loaddup2(A+K)); // DIV
            store2(B+2*K, tmp);
        } else if ( diag == 'C' ) {
            // multiply vector component by diagonal term:
            store2(B+2*K, mul2(tmp, loaddup2(A+K)));
        } else if ( diag != 'U' ) {
            printf("Error: unknown `diag` in alsatian_xtrsmLUN2\n");
        }
        for ( int I = 0; I < K; ++I )
            store2(B+2*I, fnmadd2(tmp, loaddup2(A+I), load2(B+2*I)));
    }
}

#endif

//------------------------------------------------------------------------------
#pragma mark - 1D ALSATIAN DTRSM


/// specialized DTRSM('L','L','N', diag, M, N=1, ALPHA=1.0, A, LDA, B, LDB=UNUSED);
template < char diag, typename FLOAT >
void alsatian_xtrsmLLN1(const int M, const FLOAT* A, const int lda, real* B)
{
    assert_true( M <= lda );
    for ( int K = 0; K < M; ++K )
    {
        real tmp = B[K];
        if ( diag == 'N' ) {
            tmp /= A[K];
            B[K] = tmp;
        } else if ( diag == 'I' ) {
            tmp *= A[K]; // DIV
            B[K] = tmp;
        }
        #pragma omp simd
        for ( int I = K + 1; I < M; ++I )
            B[I] -= tmp * A[I];
        A += lda;
    }
}


/// specialized DTRSM('L','L','T', diag, M, N=1, ALPHA=1.0, A, LDA, B, LDB=UNUSED);
template < char diag, typename FLOAT >
void alsatian_xtrsmLLT1(const int M, const FLOAT* A, const int lda, real* B)
{
    assert_true( M <= lda );
    A += M * lda;
    for ( int I = M-1; I >= 0; --I )
    {
        A -= lda;
        real tmp = B[I];
        #pragma omp simd
        for ( int K = I + 1; K < M; ++K )
            tmp -= A[K] * B[K];
        if ( diag == 'N' )
            tmp /= A[I];
        else if ( diag == 'I' )
            tmp *= A[I]; // DIV
        B[I] = tmp;
    }
}


/// specialized version for ORD==1
template < char diag, typename FLOAT >
void alsatian_xtrsmLUN1(const int M, const FLOAT* A, const int lda, real* B)
{
    assert_true( M <= lda );
    A += M * lda;
    for ( int K = M-1; K >= 0; --K )
    {
        A -= lda;
        real tmp = B[K];
        if ( diag == 'N' ) {
            tmp /= A[K];
            B[K] = tmp;
        } else if ( diag == 'I' ) {
            tmp *= A[K]; // DIV
            B[K] = tmp;
        } else if ( diag == 'C' ) {
            // multiply vector component by diagonal term:
            B[K] = tmp * A[K];
        } else if ( diag != 'U' ) {
            printf("Error: unknown `diag` in alsatian_xtrsmLUN1\n");
        }
        #pragma omp simd
        for ( int I = 0; I < K; ++I )
            B[I] -= tmp * A[I];
    }
}

//------------------------------------------------------------------------------
#pragma mark - ALSATIAN DTRSM for single-precision matrix argument

/// specialized DTRSM('L','L','N','U', M, N=1, ALPHA=1.0, A, LDA, B, LDB=UNUSED);
void alsatian_xtrsmLLN1U(const int M, const float* A, const int lda, real* B)
{
    assert_true( M <= lda );
    for ( int K = 0; K < M; ++K )
    {
        const real tmp = B[K];
        #pragma omp simd
        for ( int I = K + 1; I < M; ++I )
            B[I] -= tmp * A[I];
        A += lda;
    }
}

/**
 Modified DTRSM('L','U','N','I', M, N=1, ALPHA=1.0, A, LDA, B, LDB=UNUSED),
 with multiplication by the diagonal terms, when BLAS normally uses divisions
 */
void alsatian_xtrsmLUN1C(const int M, const float* A, const int lda, real* B)
{
    assert_true( M <= lda );
    A += M * lda;
    for ( int K = M-1; K >= 0; --K )
    {
        A -= lda;
        const real tmp = B[K];
        // multiply vector component by diagonal term:
        B[K] = tmp * A[K];
        #pragma omp simd
        for ( int I = 0; I < K; ++I )
            B[I] -= tmp * A[I];
    }
}

//------------------------------------------------------------------------------
#pragma mark - Unrolled ALSATIAN DTRSM for dimension 3

#if ( DIM == 3 )
/// version for M = multiple of 3
void alsatian_xtrsmLLN1U_3D(const int M, const float* A, const int lda, real* B)
{
    assert_true( M <= lda );
    assert_true( M >= 3 );
    assert_true( M%3 == 0 );
    //process columns 3 by 3
    for ( int K = 0; K < M; K += 3 )
    {
        real T0 = B[K];
        real T1 = B[K+1] - T0 * A[K+1];
        real T2 = B[K+2] - T0 * A[K+2] - T1 * A[lda+K+2];
        B[K+1] = T1;
        B[K+2] = T2;
        #pragma omp simd
        for ( int I = K + 3; I < M; ++I )
            B[I] -= T0 * A[I] + T1 * A[I+lda] + T2 * A[I+lda*2];
        A += 3*lda;
    }
}

/// this version works for M multiple of 3
void alsatian_xtrsmLUN1C_3D(const int M, const float* A, const int lda, real* B)
{
    assert_true( M <= lda );
    assert_true( M >= 3 );
    assert_true( M%3 == 0 );
    A += M * lda;
    //process columns 3 by 3
    for ( int K = M - 3; K >= 0; K -= 3 )
    {
        A -= 3*lda;
        real T2 = B[K+2];
        real T1 = B[K+1] - T2 * A[lda*2+K+1];
        real T0 = B[K  ] - T2 * A[lda*2+K  ] - T1 * A[lda+K];
        // multiply vector component by diagonal term:
        B[K+2] = T2 * A[lda*2+K+2];
        B[K+1] = T1 * A[lda  +K+1];
        B[K  ] = T0 * A[      K  ];
        #pragma omp simd
        for ( int I = 0; I < K; ++I )
            B[I] -= T0 * A[I] + T1 * A[I+lda] + T2 * A[I+lda*2];
    }
}
#endif

//------------------------------------------------------------------------------
#pragma mark - SIMD ALSATIAN DTRSM for mixed double/single-precision matrix argument

#if USE_SIMD
/// version for M = multiple of 3
void alsatian_xtrsmLLN1U_3D_SSE(const int M, const float* A, const int lda, double* B)
{
    assert_true( M <= lda );
    assert_true( M >= 3 );
    assert_true( M%3 == 0 );
    double *const end = B + M;
    //process columns 3 by 3
    int K = 0;
    for ( ; K < M-6; K += 3 )
    {
        const float* pA = A + K;
        A += 3*lda;
        double * pB = B + K;
        const vec2 t0 = loaddup2(pB);  // { T0 T0 }
        vec2 t2 = fnmadd2(t0, load2d(pA+1), loadu2(pB+1));  // { T1 B2-T0*A[2] }
        vec2 t1 = load2d(pA+1+lda);  // using upper value
        t2 = fnmadd2(unpacklo2(setzero2(), t2), t1, t2);  // { T1 T2 }
        storeu2(pB+1, t2);
        pA += 3;
        pB += 3;
        //if ( pB < end )
        {
            t1 = unpacklo2(t2, t2);  // { T1 T1 }
            t2 = unpackhi2(t2, t2);  // { T2 T2 }
            if (( M & 1 ) == ( K & 1 ))
            {
                vec2 x = fnmadd1(t0, load1d(pA), load1(pB));
                x = fnmadd1(t1, load1d(pA+lda), x); // column K+1
                x = fnmadd1(t2, load1d(pA+lda*2), x); // column K+2
                store1(pB, x);
                pA += 1;
                pB += 1;
            }
            vec2 x2 = fnmadd2(t2, load2d(pA+lda*2), loadu2(pB));
            vec2 x1 = fnmadd2(t1, load2d(pA+lda), x2);
            x2 = fnmadd2(t2, load2d(pA+2+lda*2), loadu2(pB+2));
            for ( ; pB < end-4; pB += 2 )
            {
                storeu2(pB, fnmadd2(t0, load2d(pA), x1));
                x1 = fnmadd2(t1, load2d(pA+2+lda), x2);
                x2 = fnmadd2(t2, load2d(pA+4+lda*2), loadu2(pB+4));
                pA += 2;
            }
            storeu2(pB, fnmadd2(t0, load2d(pA), x1));
            x1 = fnmadd2(t1, load2d(pA+2+lda), x2);
            storeu2(pB+2, fnmadd2(t0, load2d(pA+2), x1));
        }
        assert_true( pB+4 == B+M );
    }
    for ( ; K < M-3; K += 3 )  // same as above, for case K = M-6
    {
        const float* pA = A + K;
        A += 3*lda;
        double * pB = B + K;
        const vec2 t0 = loaddup2(pB);  // { T0 T0 }
        vec2 t2 = fnmadd2(t0, load2d(pA+1), loadu2(pB+1));  // { T1 B2-T0*A[2] }
        vec2 t1 = load2d(pA+1+lda);  // using upper value
        t2 = fnmadd2(unpacklo2(setzero2(), t2), t1, t2);  // { T1 T2 }
        storeu2(pB+1, t2);
        pA += 3;
        pB += 3;
        //if ( pB < end )
        {
            t1 = unpacklo2(t2, t2);  // { T1 T1 }
            t2 = unpackhi2(t2, t2);  // { T2 T2 }
            if (( M & 1 ) == ( K & 1 ))
            {
                vec2 x = fnmadd1(t0, load1d(pA), load1(pB));
                x = fnmadd1(t1, load1d(pA+lda), x); // column K+1
                x = fnmadd1(t2, load1d(pA+lda*2), x); // column K+2
                store1(pB, x);
                pA += 1;
                pB += 1;
            }
            for ( ; pB < end; pB += 2 )
            {
                vec2 xx = fnmadd2(t2, load2d(pA+lda*2), loadu2(pB));
                xx = fnmadd2(t1, load2d(pA+lda), xx);
                xx = fnmadd2(t0, load2d(pA), xx);
                storeu2(pB, xx);
                pA += 2;
            }
        }
        assert_true( pB == B+M );
    }
    if ( K < M )  // same as above, for case K = M-3
    {
        const float* pA = A + K;
        A += 3*lda;
        double * pB = B + K;
        const vec2 t0 = loaddup2(pB);  // { T0 T0 }
        vec2 t2 = fnmadd2(t0, load2d(pA+1), loadu2(pB+1));  // { T1 B2-T0*A[2] }
        vec2 t1 = load2d(pA+1+lda);  // using upper value
        t2 = fnmadd2(unpacklo2(setzero2(), t2), t1, t2);  // { T1 T2 }
        storeu2(pB+1, t2);
    }
}
#endif

#if USE_SIMD
/// this version works for M multiple of 3
void alsatian_xtrsmLUN1C_3D_SSE(const int M, const float* A, const int lda, double* B)
{
    assert_true( M <= lda );
    assert_true( M >= 3 );
    assert_true( M%3 == 0 );
    A += M * lda;
    // process columns 3 by 3
    int K = M - 3;
    for ( ; K > 3; K -= 3 )
    {
        A -= 3*lda;
        double * pB = B + K;
        float const* pA = A + K;
        //real T2 = B[K+2];
        //real T1 = B[K+1] - T2 * A[lda*2+K+1];
        //real T0 = B[K  ] - T2 * A[lda*2+K  ] - T1 * A[lda+K];
        vec2 t2 = loaddup2(pB+2);
        vec2 t0 = fnmadd2(t2, load2d(pA+2*lda), loadu2(pB)); // { -, T1 }
        vec2 aa = load2d(pA+lda);
        vec2 t1 = duphi2(t0); // { T1, T1 }
        t0 = duplo2(fnmadd1(t1, aa, t0)); // { T0, T0 }
        // multiply vector component by diagonal term:
        //B[K+2] = T2 * A[K+lda*2+2];
        //B[K+1] = T1 * A[K+lda  +1];
        //B[K  ] = T0 * A[K        ];
        //B[i] = B[i] - T2 * A[i+2*lda] - T1 * A[i+lda] - T0 * A[i];
        aa = blend11(load1d(pA), aa);
        store1(pB+2, mul1(t2, load1d(pA+2*lda+2)));
        storeu2(pB, mul2(blend11(t0, t1), aa));
        //if ( K > 3 )
        {
            // there should be at least 3 lines remaining
            if ( K & 1 )
            {
                --pA;
                --pB;
                vec2 x = fnmadd1(t0, load1d(pA), load1(pB));
                x = fnmadd1(t1, load1d(pA+lda), x); // column K+1
                x = fnmadd1(t2, load1d(pA+lda*2), x); // column K+2
                store1(pB, x);
            }
            pA -= 2;
            pB -= 2;
            vec2 x2 = fnmadd2(t1, load2d(pA+lda), loadu2(pB));
            vec2 x1 = fnmadd2(t2, load2d(pA+lda*2), x2);
            x2 = fnmadd2(t2, load2d(pA-2+lda*2), loadu2(pB-2));
            while ( pB > B + 2 )
            {
                storeu2(pB, fnmadd2(t0, load2d(pA), x1));
                x1 = fnmadd2(t1, load2d(pA-2+lda), x2);
                x2 = fnmadd2(t2, load2d(pA-4+lda*2), loadu2(pB-4));
                pA -= 2;
                pB -= 2;
            }
            storeu2(pB, fnmadd2(t0, load2d(pA), x1));
            x1 = fnmadd2(t1, load2d(pA-2+lda), x2);
            storeu2(pB-2, fnmadd2(t0, load2d(pA-2), x1));
        }
    }
    for ( ; K > 0; K -= 3 ) // same as above, for case K = 3
    {
        A -= 3*lda;
        double * pB = B + K;
        float const* pA = A + K;
        vec2 t2 = loaddup2(pB+2);
        vec2 t0 = fnmadd2(t2, load2d(pA+2*lda), loadu2(pB)); // { -, T1 }
        vec2 aa = load2d(pA+lda);
        vec2 t1 = duphi2(t0); // { T1, T1 }
        t0 = duplo2(fnmadd1(t1, aa, t0)); // { T0, T0 }
        // multiply vector component by diagonal term:
        aa = blend11(load1d(pA), aa);
        store1(pB+2, mul1(t2, load1d(pA+2*lda+2)));
        storeu2(pB, mul2(blend11(t0, t1), aa));
        //if ( pB > B )
        {
            // there should be at least 3 lines remaining
            if ( K & 1 )
            {
                --pA;
                --pB;
                vec2 x = fnmadd1(t0, load1d(pA), load1(pB));
                x = fnmadd1(t1, load1d(pA+lda), x); // column K+1
                x = fnmadd1(t2, load1d(pA+lda*2), x); // column K+2
                store1(pB, x);
            }
            vec2 x2 = fnmadd2(t2, load2d(pA-2+lda*2), loadu2(pB-2));
            vec2 x1 = fnmadd2(t1, load2d(pA-2+lda), x2);
            storeu2(pB-2, fnmadd2(t0, load2d(pA-2), x1));
        }
    }
    if ( M > 0 )  // same as above, for case K = 0
    {
        vec2 t2 = loaddup2(B+2);
        vec2 t0 = fnmadd2(t2, load2d(A-lda), loadu2(B)); // { -, T1 }
        vec2 aa = load2d(A-2*lda);
        vec2 t1 = duphi2(t0); // { T1, T1 }
        t0 = duplo2(fnmadd2(t1, aa, t0)); // { T0, T0 }
        aa = blend11(load1d(A-3*lda), aa);
        // multiply vector component by diagonal term:
        store1(B+2, mul1(t2, load1d(A-lda+2)));
        storeu2(B, mul2(blend11(t0, t1), aa));
    }
}
#endif


#if USE_SIMD
/// this version works for any M
void alsatian_xtrsmLLN1U_4U_SSE(const int M, const float* A, const int lda, double* B)
{
    assert_true( M <= lda );
    double *const end = B + M - 1;
    //process columns 4 by 4
    int K = 0;
    for ( ; K < M-3; K += 4 )
    {
        const float * pA = A + K;
        A += 4 * lda;
        double * pB = B + K;
#if 0
        real T0 = pB[0];
        real T1 = pB[1] - T0 * pA[1];
        real T2 = pB[2] - T0 * pA[2] - T1 * pA[lda+2];
        real T3 = pB[3] - T0 * pA[3] - T1 * pA[lda+3] - T2 * pA[2*lda+3];
        pB[1] = T1;
        pB[2] = T2;
        pB[3] = T3;
        vec2 t0{T0, T0}, t1{T1, T1}, t2{T2, T2}, t3{T3, T3};
#else
        vec2 tt = loadu2(pB);
        vec2 t0 = duplo2(tt);
        vec2 t1 = duphi2(fnmadd2(t0, load2d(pA), tt));
        store1(pB+1, t1);
        tt = fnmadd2(t1, load2d(pA+lda+2), fnmadd2(t0, load2d(pA+2), loadu2(pB+2)));
        vec2 t2 = duplo2(tt);
        store1(pB+2, t2);
        vec2 t3 = duphi2(fnmadd2(t2, load2d(pA+2*lda+2), tt));
        store1(pB+3, t3);
#endif
        pA += 4;
        pB += 4;
#if defined(__AVX__)
        vec4 tt0 = duplo2f128(cast4(t0));
        vec4 tt1 = duplo2f128(cast4(t1));
        vec4 tt2 = duplo2f128(cast4(t2));
        vec4 tt3 = duplo2f128(cast4(t3));
        #pragma ivdep
        while ( pB < end-2 )
        {
            vec4 aa = fnmadd4(tt0, load4d(pA), loadu4(pB));
            aa = fnmadd4(tt1, load4d(pA+lda), aa);
            aa = fnmadd4(tt2, load4d(pA+2*lda), aa);
            aa = fnmadd4(tt3, load4d(pA+3*lda), aa); // column K+3
            storeu4(pB, aa);
            pA += 4;
            pB += 4;
        }
#endif
#pragma clang loop unroll_count(2)
        while ( pB < end )
        {
            vec2 x = fnmadd2(t0, load2d(pA), loadu2(pB));
            x = fnmadd2(t1, load2d(pA+lda), x); // column K+1
            x = fnmadd2(t2, load2d(pA+2*lda), x); // column K+2
            x = fnmadd2(t3, load2d(pA+3*lda), x); // column K+3
            storeu2(pB, x);
            pA += 2;
            pB += 2;
        }
        if ( pB <= end )
        {
            vec2 x = fnmadd1(t0, load1d(pA), load1(pB));
            x = fnmadd1(t1, load1d(pA+lda), x);
            x = fnmadd1(t2, load1d(pA+2*lda), x);
            x = fnmadd1(t3, load1d(pA+3*lda), x); // column K+3
            store1(pB, x);
            //pA += 1;
            //pB += 1;
        }
    }
    /*
    for ( ; K < M; ++K )
    {
        const real tmp = B[K];
        for ( int I = K + 1; I < M; ++I )
            B[I] -= tmp * pA[I];
        pA += lda;
    }
     */
    // process remaining columns
    for ( ; K+1 < M; ++K )
    {
        const float * pA = A + K + 1;
        A += lda;
        double * pB = B + K + 1;
        vec2 tt = loaddup2(B+K);
        for ( pB = B + K + 1; pB < end; pB += 2 )
        {
            storeu2(pB, fnmadd2(tt, load2d(pA), loadu2(pB)));
            pA += 2;
        }
        if ( pB <= end )
        {
            store1(pB, fnmadd1(tt, load1d(pA), load1(pB)));
            //pA += 1;
            //pB += 1;
        }
    }
}
#endif


#if USE_SIMD
/// this version works for any M
void alsatian_xtrsmLUN1C_4U_SSE(const int M, const float* A, const int lda, double* B)
{
    assert_true( M <= lda );
    A += M * lda;
    //process columns 4 by 4
    int K = M;
    while ( K > 3 )
    {
        K -= 4;
        A -= 4 * lda;
        double * pB = B + K;
        float const* pA = A + K;
#if 0
        real T3 = pB[3];
        real T2 = pB[2] - T3 * pA[3*lda+2];
        real T1 = pB[1] - T3 * pA[3*lda+1] - T2 * pA[2*lda+1];
        real T0 = pB[0] - T3 * pA[3*lda  ] - T2 * pA[2*lda] - T1 * pA[lda];
        vec2 t0{T0, T0}, t1{T1, T1}, t2{T2, T2}, t3{T3, T3};
#else
        vec2 b2 = loadu2(pB+2);
        vec2 t3 = duphi2(b2);
        vec2 t2 = duplo2(fnmadd1(t3, load1d(pA+3*lda+2), b2));
        vec2 tt = fnmadd2(t2, load2d(pA+2*lda), fnmadd2(t3, load2d(pA+3*lda), loadu2(pB)));
        vec2 t1 = duphi2(tt);
        vec2 t0 = duplo2(fnmadd1(t1, load1d(pA+lda), tt));
#endif
        //pB[3] = pA[3*lda+3] * T3;
        //pB[2] = pA[2*lda+2] * T2;
        //pB[1] = pA[1*lda+1] * T1;
        //pB[0] = pA[      0] * T0;
        store1(pB+3, mul1(t3, load1d(pA+3*lda+3)));
        store1(pB+2, mul1(t2, load1d(pA+2*lda+2)));
        store1(pB+1, mul1(t1, load1d(pA+lda+1)));
        store1(pB, mul1(t0, load1d(pA)));
#if defined(__AVX__)
        vec4 tt0 = duplo2f128(cast4(t0));
        vec4 tt1 = duplo2f128(cast4(t1));
        vec4 tt2 = duplo2f128(cast4(t2));
        vec4 tt3 = duplo2f128(cast4(t3));
        #pragma ivdep
        while ( pB > B+3 )
        {
            pA -= 4;
            pB -= 4;
            vec4 aa = fnmadd4(tt0, load4d(pA), loadu4(pB));
            aa = fnmadd4(tt1, load4d(pA+lda), aa);
            aa = fnmadd4(tt2, load4d(pA+2*lda), aa);
            aa = fnmadd4(tt3, load4d(pA+3*lda), aa);
            storeu4(pB, aa);
        }
#endif
        #pragma ivdep
        while ( pB > B+1 )
        {
            pA -= 2;
            pB -= 2;
            vec2 x = fnmadd2(t0, load2d(pA), loadu2(pB));
            x = fnmadd2(t1, load2d(pA+lda), x); // column K+1
            x = fnmadd2(t2, load2d(pA+lda*2), x); // column K+2
            x = fnmadd2(t3, load2d(pA+lda*3), x); // column K+3
            storeu2(pB, x);
        }
        if ( pB > B )
        {
            --pA;
            --pB;
            vec2 x = fnmadd1(t0, load1d(pA), load1(pB));
            x = fnmadd1(t1, load1d(pA+lda), x); // column K+1
            x = fnmadd1(t2, load1d(pA+lda*2), x); // column K+2
            x = fnmadd1(t3, load1d(pA+lda*3), x); // column K+3
            store1(pB, x);
        }
    }
    /*
    while ( K > 0 ) {
        K -= 1;
        const real tmp = B[K] * A[K];
        B[K] = tmp;
        for ( int I = 0; I < K; ++I )
            B[I] -= tmp * A[I];
    }
    */
    // process remaining columns:
    while ( K > 0 )
    {
        K -= 1;
        A -= lda;
        double * pB = B + K;
        float const* pA = A + K;
        vec2 tt = loaddup2(pB);
        store1(pB, mul1(tt, load1d(pA)));
        #pragma ivdep
        while ( pB > B+1 )
        {
            pA -= 2;
            pB -= 2;
            storeu2(pB, fnmadd2(tt, load2d(pA), loadu2(pB)));
        }
        if ( pB > B )
        {
            pA -= 1;
            pB -= 1;
            store1(pB, fnmadd1(tt, load1d(pA), load1(pB)));
        }
    }
}
#endif

//------------------------------------------------------------------------------
#pragma mark - SIMD ALSATIAN DTRSM for single-precision matrix argument

#if USE_SIMD
/// version for M = multiple of 3
/// Severely suboptimal: float are processed 2 by 2, only using half of the SSE register
void alsatian_xtrsmLLN1U_3D_SSE(const int M, const float* A, const int lda, float* B)
{
    assert_true( M <= lda );
    assert_true( M >= 3 );
    assert_true( M%3 == 0 );
    float *const end = B + M;
    //process columns 3 by 3
    for ( int K = 0; K < M; K += 3 )
    {
        float const* pA = A + K;
        float * pB = B + K;
        const vec2f t0 = loaddupf(pB);  // { T0 T0 }
        vec2f t2 = fnmadd2f(t0, load2f(pA+1), load2f(pB+1));  // { T1 B2-T0*A[2] }
        vec2f t1 = load2f(pA+1+lda);  // using upper value
        t2 = fnmadd2f(unpacklo2f(setzero2f(), t2), t1, t2);  // { T1 T2 }
        store2f(pB+1, t2);
        if ( pB < end )
        {
            t1 = unpacklo2f(t2, t2);  // { T1 T1 }
            t2 = unpackhi2f(t2, t2);  // { T2 T2 }
            pA += 3;
            pB += 3;
            if (( M & 1 ) == ( K & 1 ))
            {
                vec2f x = fnmadd2f(t0, load1f(pA), load1f(pB));
                x = fnmadd2f(t1, load1f(pA+lda), x); // column K+1
                x = fnmadd2f(t2, load1f(pA+lda*2), x); // column K+2
                store1f(pB, x);
                pA += 1;
                pB += 1;
            }
            #pragma ivdep
            while ( pB < end )
            {
                vec2f x = fnmadd2f(t0, load2f(pA), load2f(pB));
                x = fnmadd2f(t1, load2f(pA+lda), x); // column K+1
                x = fnmadd2f(t2, load2f(pA+lda*2), x); // column K+2
                store2f(pB, x);
                pA += 2;
                pB += 2;
            }
        }
        A += 3*lda;
    }
}
#endif

#if USE_SIMD
/// this version works for M multiple of 3, suboptimal since using only 2 floats
/// Severely suboptimal: float are processed 2 by 2, only using half of the SSE register
void alsatian_xtrsmLUN1C_3D_SSE(const int M, const float* A, const int lda, float* B)
{
    assert_true( M <= lda );
    assert_true( M >= 3 );
    assert_true( M%3 == 0 );
    A += M * lda;
    //process columns 3 by 3
    for ( int K = M - 3; K >= 0; K -= 3 )
    {
        A -= 3*lda;
        float * pB = B + K;
        float const* pA = A + K;
        vec2f t2 = loaddupf(pB+2); // { T2, T2 }
        vec2f t0 = fnmadd2f(t2, load2f(pA+2*lda), load2f(pB)); // { -, T1/A }
        vec2f aa = load2f(pA+lda);
        vec2f t1 = duphi2f(t0); // { T1, T1 }
        t0 = duplo2f(fnmadd2f(t1, aa, t0)); // { T0, T0 }
        store1f(pB+2, mul2f(t2, load1f(pA+2*lda+2)));
        store1f(pB+1, mul2f(t1, unpackhi2f(aa, aa)));
        store1f(pB, mul2f(t0, load1f(pA)));
        if ( pB > B )
        {
            if ( K & 1 )
            {
                --pA;
                --pB;
                vec2f x = fnmadd2f(t0, load1f(pA), load1f(pB));
                x = fnmadd2f(t1, load1f(pA+lda), x); // column K+1
                x = fnmadd2f(t2, load1f(pA+lda*2), x); // column K+2
                store1f(pB, x);
            }
            #pragma ivdep
            while ( pB > B )
            {
                pA -= 2;
                pB -= 2;
                vec2f x = fnmadd2f(t0, load2f(pA), load2f(pB));
                x = fnmadd2f(t1, load2f(pA+lda), x); // column K+1
                x = fnmadd2f(t2, load2f(pA+lda*2), x); // column K+2
                store2f(pB, x);
            }
        }
    }
}
#endif


#if USE_SIMD
/// this version works for any M, but there is another version for dimension 3 below
/// Severely suboptimal: float are processed 2 by 2, only using half of the SSE register
void alsatian_strsmLLN1U_SSE(const int M, const float* A, const int lda, float* B)
{
    assert_true( M <= lda );
    assert_true( M >= DIM );
    float *const end = B + M;
    int K = 0;
#if ( DIM & 1 )
    // process one column to leave an even number
    if ( M & 1 )
    {
        vec2f t = loaddupf(B);
        float const* pA = A + 1;
        // there is an even number of scalars, fitting perfectly:
        for ( float * pB = B + 1; pB < end; pB += 2 )
        {
            store2f(pB, fnmadd2f(t, load2f(pA), load2f(pB)));
            pA += 2;
        }
        A += lda;
        ++K;
    }
#endif
    // process columns 2 by 2
    for ( ; K < M; K += 2 )
    {
        float const* pA = A + K;
        A += 2*lda;
        vec2f t0 = load2f(B + K);
        vec2f t1 = load2f(pA); // will use upper value
        t1 = fnmadd2f(unpacklo2f(setzero2f(), t0), t1, t0);
        store2f(B + K, t1);
        t0 = unpacklo2f(t0, t0);
        t1 = unpackhi2f(t1, t1);
        pA += 2;
        for ( float * pB = B + K + 2; pB < end; pB += 2 )
        {
            vec2f x = fnmadd2f(t0, load2f(pA), load2f(pB));
            x = fnmadd2f(t1, load2f(pA+lda), x);
            store2f(pB, x);
            pA += 2;
        }
    }
}
#endif


#if USE_SIMD
/// this version works for any M, but there is another version for dimension 3 below
/// Severely suboptimal: float are processed 2 by 2, only using half of the SSE register
void alsatian_strsmLUN1C_SSE(const int M, const float* A, const int lda, float* B)
{
    assert_true( M <= lda );
    assert_true( M >= DIM );
    A += M * lda;
    int K = M - 1;
#if ( DIM & 1 )
    if ( M & 1 )
    {
        A -= lda;
        float * end = B + K;
        float const* pA = A;
        vec2f t0 = load1f(end);
        t0 = unpacklo2f(t0, t0); // { T0, T0 }
        vec2f aa = mul2f(t0, load1f(A+K));
        // there is an even number of scalars remaining, fitting perfectly:
        for ( float * pB = B ; pB < end; pB += 2 )
        {
            store2f(pB, fnmadd2f(t0, load2f(pA), load2f(pB)));
            pA += 2;
        }
        store1f(end, aa);
        --K;
    }
#endif
    --K;
    //process columns 2 by 2
    for ( ; K >= 0; K -= 2 )
    {
        A -= 2*lda;
        float * end = B + K;
        float const* pA = A;
        //real T1 = B[K+1];
        //real T0 = B[K  ] - T1 * A[lda+K];
        vec2f t0 = load2f(end); // { B[K], T1 }
        vec2f aa = load2f(pA+lda+K);
        vec2f t1 = unpackhi2f(t0, t0); // { T1, T1 }
        t0 = fnmadd2f(t1, aa, t0); // { T0, ? }
        t0 = unpacklo2f(t0, t0); // { T0, T0 }
        // multiply vector component by diagonal term:
        //B[K+1] = T1 * A[lda+K+1];
        //B[K  ] = T0 * A[    K  ];
        aa = mul2f(blend11f(t0, t1), blend11f(load2f(pA+K), aa));
        for ( float * pB = B; pB < end; pB += 2 )
        {
            vec2f x = fnmadd2f(t0, load2f(pA), load2f(pB));
            x = fnmadd2f(t1, load2f(pA+lda), x);
            store2f(pB, x);
            pA += 2;
        }
        store2f(end, aa);
    }
}
#endif

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
        // Solve L**T *X = B, overwriting B with X. ALPHA = 1.0
        blas::xtrsm('L', 'L', 'T', 'N', N, NRHS, 1.0, A, LDA, B, LDB);
    }
}


template < int ORD >
void iso_xpotrsL_lapack(const int N, const real* A, const int LDA, real* B)
{
    /*
     we cannot call lapack::DPOTRS('L', N, 1, A, LDA, B, N, &info);
     because the coordinates of the vector 'Y' are not contiguous but offset by 'ORD'.
     But calling DTBSV gets the required work done.
     */
    real * tmp = new_real(N);
    for ( int d = 0; d < ORD; ++d )
    {
        for ( int u = 0; u < N; ++u )
            tmp[u] = B[d+ORD*u];
        // Solve L*X = B, overwriting B with X. ALPHA = 1.0
        blas::xtrsm('L', 'L', 'N', 'N', N, 1, 1.0, A, LDA, tmp, N);
        // Solve U*X = B, overwriting B with X. ALPHA = 1.0
        blas::xtrsm('L', 'L', 'T', 'N', N, 1, 1.0, A, LDA, tmp, N);
        for ( int u = 0; u < N; ++u )
            B[d+ORD*u] = tmp[u];
    }
    free_real(tmp);
}


template < int ORD >
void iso_xpotrsL(const int N, const real* A, const int LDA, real* B)
{
    //return iso_xpotrsL_lapack<ORD>(N, A, LDA, B);
    // Solve L*X = B, overwriting B with X. ALPHA = 1.0
    iso_xtrsmLLN<ORD,'N'>(N, A, LDA, B);
    // Solve U*X = B, overwriting B with X. ALPHA = 1.0
    iso_xtrsmLLT<ORD,'N'>(N, A, LDA, B);
}



/**
 This calls lapack::xpotf2() and then inverts the diagonal terms
*/
void alsatian_xpotf2L(const int N, real* A, const int LDA, int* INFO)
{
    lapack::xpotf2('L', N, A, LDA, INFO);
    if ( *INFO )
    {
        std::cerr << "lapack::xpotf2 failed (code " << *INFO << ")\n";
    }
    else
    {
        const int S = LDA+1;
        for ( int u = 0; u < N*S; u += S )
            A[u] = real(1) / A[u];
    }
}


template < int ORD >
void alsatian_iso_xpotrsLref(const int N, const real* A, const int LDA, real* B)
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
void alsatian_iso_xpotrsL(const int N, const real* A, const int LDA, real* B)
{
#if XTRSM_USES_AVX
    if ( ORD == 3 )
    {
        alsatian_iso3_xtrsmLLN_AVX<'I'>(N, A, LDA, B);
        alsatian_iso3_xtrsmLLT_AVX<'I'>(N, A, LDA, B);
    } else
#endif
#if XTRSM_USES_SSE3
        if ( ORD == 2 )
    {
        alsatian_iso2_xtrsmLLN_SIMD<'I'>(N, A, LDA, B);
        alsatian_iso2_xtrsmLLT_SIMD<'I'>(N, A, LDA, B);
    } else
#endif
        if ( ORD == 1 )
    {
        alsatian_xtrsmLLN1<'I'>(N, A, LDA, B);
        alsatian_xtrsmLLT1<'I'>(N, A, LDA, B);
    }
    else
    {
        // Solve L*X = B, overwriting B with X. ALPHA = 1.0
        iso_xtrsmLLN<ORD,'I'>(N, A, LDA, B);
        // Solve U*X = B, overwriting B with X. ALPHA = 1.0
        iso_xtrsmLLT<ORD,'I'>(N, A, LDA, B);
    }
}


void alsatian_xpotrsL(const int N, const float* A, const int LDA, real* B)
{
#if 1
    // Solve L*X = B, overwriting B with X
    alsatian_xtrsmLLN1<'I'>(N, A, LDA, B);
    // Solve U*X = B, overwriting B with X
    alsatian_xtrsmLLT1<'I'>(N, A, LDA, B);
#else
    iso_xtrsmLLN<1,'I'>(N, A, LDA, B);
    iso_xtrsmLLT<1,'I'>(N, A, LDA, B);
#endif
}


//------------------------------------------------------------------------------
#pragma mark - LAPACK:LASWP equivalent swap routines

template < int ORD, typename REAL >
inline void xswap(REAL* A, REAL* B)
{
    // this is adjusted for the one-based array index:
    for ( int d = -ORD; d < 0; ++d )
    {
        REAL t = A[d];
        A[d] = B[d];
        B[d] = t;
    }
}


/**
 for ORD=1, this is equivalent to
    DLASWP(N=1, A, LDA=UNUSED, K1, K2, IPIV, INCX)
 for ORD>1, this performs a series of row interchanges on the matrix A,
 always swapping chunks of 'ORD' scalars.

 Note that, following LAPACK's convention, K1, K2 and IPIV contain one-based array indices

 if ( INCX > 0 )
 {
     IX0 = K1
     I1 = K1
     I2 = K2
     INC = 1
 }
 else if ( INCX < 0 )
 {
     IX0 = 1 + ( 1-K2 )*INCX
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
template < int ORD, typename REAL >
void iso_xlaswp(REAL* A, int K1, int K2, const int* IPIV, int INCX)
{
    /*
     as per LAPACK's convention, K1, K2 and IPIV contain one-based array indices
     and by shifting IPIV we fall back on the C-convention of zero-based array
     */
    --IPIV;
    if ( INCX == 1 )
    {
        for ( int I = K1; I <= K2; ++I )
        {
            int P = IPIV[I];
            if ( P != I )
                xswap<ORD>(A+ORD*I, A+ORD*P);
        }
    }
    else if ( INCX > 0 )
    {
        for ( int I = K1; I <= K2; ++I )
        {
            int P = IPIV[INCX*I];
            if ( P != I )
                xswap<ORD>(A+ORD*I, A+ORD*P);
        }
    }
    else if ( INCX < 0 )
    {
        int IX = 1 + ( 1 - K2 )*INCX;
        for ( int I = K2; I >= K1; --I )
        {
            int P = IPIV[IX];
            if ( P != I )
                xswap<ORD>(A+ORD*I, A+ORD*P);
            IX += INCX;
        }
    }
}


/// Should be equivalent to DLASWP(N=1, A, LDA=UNUSED, K1, K2, IPIV, INCX=1)
template < typename REAL >
void xlaswp1(REAL* A, int K1, int K2, const int* IPIV)
{
    /*
     as per LAPACK's convention, K1, K2 and IPIV contain one-based array indices
     and by shifting IPIV we fall back on the C-convention of zero-based array
     */
    --A;
    --IPIV;
    for ( int I = K1; I <= K2; ++I )
    {
        const int P = IPIV[I];
        /*
         Useless swaps could be avoided, but adding branch points to the code
         might decreased performance, so the tradeoff is not clear
         */
        //if ( P != I )
        {
            REAL t = A[I];
            A[I] = A[P];
            A[P] = t;
        }
    }
}


/// Counting the number of permutations performed by xlaswp1()
size_t count_swaps(int K1, int K2, const int* IPIV)
{
    /*
     as per LAPACK's convention, K1, K2 and IPIV contain one-based array indices
     and by shifting IPIV we fall back on the C-convention of zero-based array
     */
    --IPIV;
    size_t cnt = 0;
    for ( int I = K1; I <= K2; ++I )
    {
        const int P = IPIV[I];
        cnt += ( P != I );
    }
    return cnt;
}

void reset_pivot(int K1, int K2, int* IPIV)
{
    /*
     as per LAPACK's convention, K1, K2 and IPIV contain one-based array indices
     and by shifting IPIV we fall back on the C-convention of zero-based array
     */
    --IPIV;
    for ( int I = K1; I <= K2; ++I )
        IPIV[I] = I;
}

//------------------------------------------------------------------------------
#pragma mark - LAPACK-STYLE ROUTINES for General Matrix LU factorization

/**
 Call lapack::xgetf2() to compute the LU factorization, and inverts diagonal terms
*/
void alsatian_xgetf2(const int N, real* A, const int LDA, int* IPIV, int* INFO)
{
    lapack::xgetf2(N, N, A, LDA, IPIV, INFO);
    if ( 0 == *INFO )
    {
        for ( int j = 0; j < N; ++j )
        {
            real tmp = real(1) / A[j];
            for ( int i = 0; i < j; ++i )
                A[i] *= tmp;
            A[j] = tmp;
            A += LDA;
        }
    }
}

/// factorize a single-precision matrix, with inverted diagonal terms
void alsatian_sgetf2(const int N, float* A, const int LDA, int* IPIV, int* INFO)
{
    lapack::sgetf2(N, N, A, LDA, IPIV, INFO);
    if ( 0 == *INFO )
    {
        for ( int j = 0; j < N; ++j )
        {
            float tmp = 1.0f / A[j];
            for ( int i = 0; i < j; ++i )
                A[i] *= tmp;
            A[j] = tmp;
            A += LDA;
        }
    }
}

/// Following LAPACK's interface
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


/// departing from the standard LAPACK interface, specialized for TRANS='N'
template < typename REAL >
void lapack_xgetrsN(int N, int NRHS, const REAL* A, int LDA, const int* IPIV, REAL* B, int LDB, int* INFO)
{
    *INFO = 0;
    for ( int i = 0; i < NRHS; ++i )
    {
        // Apply row interchanges to the right hand side.
        xlaswp1(B, 1, N, IPIV); //iso_xlaswp<1>(B, 1, N, IPIV, 1);
        // Solve L*X = B, overwriting B with X.
        alsatian_xtrsmLLN1<'U'>(N, A, LDA, B);
        // Solve U*X = B, overwriting B with X.
        alsatian_xtrsmLUN1<'N'>(N, A, LDA, B);
        B += LDB;
    }
}

/// simplified calling interface equivalent to xgetrs(TRANS='N', N, NRHS==1, ...)
template < typename REAL >
void lapack_xgetrsN(int N, const REAL* A, int LDA, const int* IPIV, REAL* B)
{
    // Apply row interchanges to the right hand side.
    xlaswp1(B, 1, N, IPIV); //iso_xlaswp<1>(B, 1, N, IPIV, 1);
    // Solve L*X = B, overwriting B with X.
    alsatian_xtrsmLLN1<'U'>(N, A, LDA, B);
    // Solve U*X = B, overwriting B with X.
    alsatian_xtrsmLUN1<'N'>(N, A, LDA, B);
}


/// simplified calling interface equivalent to xgetrs(TRANS='N', N, NRHS==1, ...)
/// this version assumes that diagonal terms have been inverted, and that A is single precision array
/**
 This is used to apply the full block preconditionner, stored in single precision.
 We could use non-temporal loads for the matrix A, since it will not fit in the cache,
 but the vector B should be cached!
 */
void alsatian_xgetrsN(int N, const float* A, int LDA, const int* IPIV, real* B)
{
    // Apply row interchanges to the right hand side.
    xlaswp1(B, 1, N, IPIV); //iso_xlaswp<1>(B, 1, N, IPIV, 1);
    // Solve L*X = B, overwriting B with X.
    alsatian_xtrsmLLN1U(N, A, LDA, B);
    // Solve U*X = B, overwriting B with X.
    alsatian_xtrsmLUN1C(N, A, LDA, B);
}


#if USE_SIMD
/// This is used to apply the full block preconditionner, stored in single precision.
void alsatian_xgetrsN_SSE(int N, const float* A, int LDA, const int* IPIV, real* B)
{
    // Apply row interchanges to the right hand side.
    xlaswp1(B, 1, N, IPIV);
#if REAL_IS_DOUBLE
    // Solve L*X = B, overwriting B with X.
    alsatian_xtrsmLLN1U_4U_SSE(N, A, LDA, B);
    // Solve U*X = B, overwriting B with X.
    alsatian_xtrsmLUN1C_4U_SSE(N, A, LDA, B);
#elif ( DIM == 3 )
    // Solve L*X = B, overwriting B with X.
    alsatian_xtrsmLLN1U_3D(N, A, LDA, B);
    // Solve U*X = B, overwriting B with X.
    alsatian_xtrsmLUN1C_3D(N, A, LDA, B);
#else
    // Solve L*X = B, overwriting B with X.
    alsatian_xtrsmLLN1U(N, A, LDA, B);
    // Solve U*X = B, overwriting B with X.
    alsatian_xtrsmLUN1C(N, A, LDA, B);
#endif
}

/// This compares results obtained by standard and optimized routines
void alsatian_xgetrsN_SSE_CHECK(int N, const float* A, int LDA, const int* IPIV, real* B)
{
    real * T = new_real(N);
    copy_real(N, B, T);

    alsatian_xgetrsN_SSE(N, A, LDA, IPIV, B);

    // try second method:
    alsatian_xgetrsN(N, A, LDA, IPIV, T);
    real err = blas::difference(N, B, T);
    // this test is true even if `err` is not-a-number
    if ( not ( err < 0.01 ) )
    {
        VecPrint::edges("xgetrs", N, B, 3);
        VecPrint::edges("ref---", N, T, 3);
    }// else printf("\n xgetrs %3i okay!", N);
    free_real(T);
}

#endif

template < int ORD >
void iso_getrsN_lapack(const int N, const real* A, const int LDA, const int* IPIV, real* B)
{
    /*
     we cannot call lapack::DGETRS('N', bks, 1, mec->pblock(), bks, mec->pivot(), Y, bks, &info);
     because the coordinates of the vector 'Y' are not contiguous but offset by 'ORD'.
     But calling DTBSV gets the required work done.
     */
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
}

/// version of xgetrs('N', ...) for ORD interleaved vectors in right-hand-side B
template < int ORD >
void iso_xgetrsN(const int N, const real* A, const int LDA, const int* IPIV, real* B)
{
    //return iso_getrsN_lapack<ORD>(N, A, LDA, IPIV, B);
    // Apply row interchanges to the right hand side.
    iso_xlaswp<ORD>(B, 1, N, IPIV, 1);
    // Solve L*X = B, overwriting B with X.
    iso_xtrsmLLN<ORD,'U'>(N, A, LDA, B);
    // Solve U*X = B, overwriting B with X.
    iso_xtrsmLUN<ORD,'N'>(N, A, LDA, B);
}


/// version of xgetrs('N', ...) for ORD interleaved vectors in right-hand-side B
template < int ORD >
void alsatian_iso_xgetrsN(const int N, const real* A, const int LDA, const int* IPIV, real* B)
{
    // Apply row interchanges to the right hand side.
    iso_xlaswp<ORD>(B, 1, N, IPIV, 1);
#if XTRSM_USES_AVX
    if ( ORD == 3 )
    {
        // Solve L*X = B, overwriting B with X.
        alsatian_iso3_xtrsmLLN_AVX<'U'>(N, A, LDA, B);
        // Solve U*X = B, overwriting B with X.
        alsatian_iso3_xtrsmLUN_AVX<'C'>(N, A, LDA, B);
    } else
#endif
#if XTRSM_USES_SSE3
        if ( ORD == 2 )
    {
        // Solve L*X = B, overwriting B with X.
        alsatian_iso2_xtrsmLLN_SIMD<'U'>(N, A, LDA, B);
        // Solve U*X = B, overwriting B with X.
        alsatian_iso2_xtrsmLUN_SIMD<'C'>(N, A, LDA, B);
    } else
#endif
        if ( ORD == 1 )
    {
        // Solve L*X = B, overwriting B with X.
        alsatian_xtrsmLLN1<'U'>(N, (float*)A, LDA, B);
        // Solve U*X = B, overwriting B with X.
        alsatian_xtrsmLUN1<'C'>(N, (float*)A, LDA, B);
    }
    else
    {
        iso_xtrsmLLN<ORD,'U'>(N, (float*)A, LDA, B);
        iso_xtrsmLUN<ORD,'C'>(N, (float*)A, LDA, B);
    }
}

#endif
