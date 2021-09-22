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


/* CLANG offers builtin non temporal load functions */
#ifdef __clang__
template < typename T >
static inline float ntload(T const* x) { return __builtin_nontemporal_load(x); }
#else
template < typename T >
static inline float ntload(T const* x) { return *x; }
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
template < int ORD, char diag, typename REAL >
void iso_xtrsmLLN(const int M, const REAL* A, const int lda, REAL* B)
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
template < int ORD, char diag, typename REAL >
void iso_xtrsmLLT(const int M, const REAL* A, const int lda, REAL* B)
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
template < int ORD, char diag, typename REAL >
void iso_xtrsmLUN(const int M, const REAL* A, const int lda, REAL* B)
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
            } else {
                for ( int d = 0; d < ORD; ++d )
                    tmp[d] = B[ORD*K+d];
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
template < int ORD, char diag, typename REAL >
void iso_xtrsmLUT(const int M, const REAL* A, const int lda, REAL* B)
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
#pragma mark - 3D SIMD ALSATIAN DTRSM

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
void alsatian_xtrsmLLN3(const int M, const double* A, const int lda, double* B)
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
void alsatian_xtrsmLLT3(const int M, const double* A, const int lda, double* B)
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
void alsatian_xtrsmLUN3(const int M, const double* A, const int lda, double* B)
{
    assert_true( M <= lda );
    A += M * lda;
    for ( int K = M-1; K >= 0; --K )
    {
        A -= lda;
        vec4 T = loadu4(B+3*K);
        vec4 ori = T;
        if ( diag == 'N' ) {
            T = div4(ori, broadcast1(A+K));
            ori = blend31(T, ori);  // blend to keep 4th value!
        } else if ( diag == 'I' ) {
            T = mul4(ori, broadcast1(A+K)); // DIV
            ori = blend31(T, ori);
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
void alsatian_xtrsmLLN2(const int M, const double* A, const int lda, double* B)
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
void alsatian_xtrsmLLT2(const int M, const double* A, const int lda, double* B)
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
void alsatian_xtrsmLUN2(const int M, const double* A, const int lda, double* B)
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
        }
        for ( int I = 0; I < K; ++I )
            store2(B+2*I, fnmadd2(tmp, loaddup2(A+I), load2(B+2*I)));
    }
}

#endif

//------------------------------------------------------------------------------
#pragma mark - 1D ALSATIAN DTRSM


/// specialized DTRSM('L','L','N', diag, M, N=1, ALPHA=1.0, A, LDA, B, LDB=UNUSED);
template < char diag, typename REAL >
void alsatian_xtrsmLLN1(const int M, const REAL* A, const int lda, REAL* B)
{
    assert_true( M <= lda );
    for ( int K = 0; K < M; ++K )
    {
        REAL tmp = B[K];
        if ( diag == 'N' ) {
            tmp /= A[K];
            B[K] = tmp;
        } else if ( diag == 'I' ) {
            tmp *= A[K]; // DIV
            B[K] = tmp;
        }
        # pragma ivdep
        for ( int I = K + 1; I < M; ++I )
            B[I] -= tmp * A[I];
        A += lda;
    }
}


/// specialized DTRSM('L','L','T', diag, M, N=1, ALPHA=1.0, A, LDA, B, LDB=UNUSED);
template < char diag, typename REAL >
void alsatian_xtrsmLLT1(const int M, const REAL* A, const int lda, REAL* B)
{
    assert_true( M <= lda );
    A += M * lda;
    for ( int I = M-1; I >= 0; --I )
    {
        A -= lda;
        REAL tmp = B[I];
        # pragma ivdep
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
template < char diag, typename REAL >
void alsatian_xtrsmLUN1(const int M, const REAL* A, const int lda, REAL* B)
{
    assert_true( M <= lda );
    A += M * lda;
    for ( int K = M-1; K >= 0; --K )
    {
        A -= lda;
        REAL tmp = B[K];
        if ( diag == 'N' ) {
            tmp /= A[K];
            B[K] = tmp;
        } else if ( diag == 'I' ) {
            tmp *= A[K]; // DIV
            B[K] = tmp;
        }
        # pragma ivdep
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
        # pragma ivdep
        for ( int I = K + 1; I < M; ++I )
            B[I] -= tmp * A[I];
        A += lda;
    }
}

/**
 Modified DTRSM('L','U','N','I', M, N=1, ALPHA=1.0, A, LDA, B, LDB=UNUSED),
 with multiplication by the diagonal terms, when BLAS normally uses divisions
 */
void alsatian_xtrsmLUN1I(const int M, const float* A, const int lda, real* B)
{
    assert_true( M <= lda );
    A += M * lda;
    for ( int K = M-1; K >= 0; --K )
    {
        A -= lda;
        const real tmp = B[K] * A[K];
        B[K] = tmp;
        # pragma ivdep
        for ( int I = 0; I < K; ++I )
            B[I] -= tmp * A[I];
    }
}

//------------------------------------------------------------------------------
#pragma mark - Unrolled ALSATIAN DTRSM for DIM == 3

#if ( DIM == 3 )
/// version for M = multiple of 3
void alsatian_xtrsmLLN1U_3D(const int M, const float* A, const int lda, real* B)
{
    assert_true( M <= lda );
    assert_true( M >= 2*DIM ); // 3D:6
    assert_true( M%3 == 0 );
    //process columns 3 by 3
    for ( int K = 0; K < M; K += 3 )
    {
        real T0 = B[K];
        real T1 = B[K+1] - T0 * A[K+1];
        real T2 = B[K+2] - T0 * A[K+2] - T1 * A[lda+K+2];
        B[K+1] = T1;
        B[K+2] = T2;
        # pragma ivdep
        for ( int I = K + 3; I < M; ++I )
            B[I] -= T0 * A[I] + T1 * A[I+lda] + T2 * A[I+lda*2];
        A += 3*lda;
    }
}

/// this version works for M multiple of 3
void alsatian_xtrsmLUN1I_3D(const int M, const float* A, const int lda, real* B)
{
    assert_true( M <= lda );
    assert_true( M >= 2*DIM ); // 2D:4 3D:6
    assert_true( M%3 == 0 );
    A += M * lda;
    //process columns 3 by 3
    for ( int K = M - 3; K >= 0; K -= 3 )
    {
        A -= 3*lda;
        real T2 = A[lda*2+K+2] * B[K+2];
        real T1 = A[lda  +K+1] * ( B[K+1] - T2 * A[lda*2+K+1] );
        real T0 = A[      K  ] * ( B[K  ] - T2 * A[lda*2+K  ] - T1 * A[lda+K] );
        B[K  ] = T0;
        B[K+1] = T1;
        B[K+2] = T2;
        # pragma ivdep
        for ( int I = 0; I < K; ++I )
            B[I] -= T0 * A[I] + T1 * A[I+lda] + T2 * A[I+lda*2];
    }
}
#endif

//------------------------------------------------------------------------------
#pragma mark - SIMD ALSATIAN DTRSM for mixed double/single-precision matrix argument

#if ( DIM < 3 ) && defined(__SSE3__)
/// this version works for any M, but there is another version for DIM==3 below
void alsatian_xtrsmLLN1U_SSE(const int M, const float* A, const int lda, double* B)
{
    assert_true( M <= lda );
    assert_true( M >= 2*DIM ); // 2D:4
    int K = 0;
#if ( DIM & 1 )
    // process one column to leave an even number
    if ( M & 1 )
    {
        vec2 t = loaddup2(B);
        double * pB = B + 1;
        double * end = B + M;
        float const* pA = A + 1;
        // there is an even number of scalars, fitting perfectly:
        # pragma ivdep
        while ( pB < end )
        {
            storeu2(pB, fnmadd2(t, load2d(pA), loadu2(pB)));
            pA += 2;
            pB += 2;
        }
        A += lda;
        ++K;
    }
#endif
    // process columns 2 by 2
    for ( ; K < M; K += 2 )
    {
        double * pB = B + K;
        double * end = B + M;
        float const* pA = A + K;
        /*
        const real tmp = B[K];
        for ( int I = K + 1; I < M; ++I )
            B[I] -= tmp * A[I];
         */
        vec2 b = loadu2(pB);
        vec2 c = load2d(pA); // will use upper value
        c = fnmadd2(unpacklo2(setzero2(), b), c, b);
        storeu2(pB, c);
        b = unpacklo2(b, b);
        c = unpackhi2(c, c);
        pB += 2;
        pA += 2;
        # pragma ivdep
        while ( pB < end )
        {
            vec2 x = fnmadd2(b, load2d(pA), loadu2(pB));
            x = fnmadd2(c, load2d(pA+lda), x);
            storeu2(pB, x);
            pA += 2;
            pB += 2;
        }
        A += 2*lda;
    }
}
#endif

#if ( DIM == 3 ) && defined(__SSE3__)
/// version for M = multiple of 3
void alsatian_xtrsmLLN1U_SSE(const int M, const float* pA, const int lda, double* B)
{
    assert_true( M <= lda );
    assert_true( M >= 2*DIM ); // 3D:6
    assert_true( M%3 == 0 );
    //process columns 3 by 3
    for ( int K = 0; K < M; K += 3 )
    {
        double * pB = B + K;
        double * end = B + M;
        pA += K;
        /*
        const real tmp = B[K];
        for ( int I = K + 1; I < M; ++I )
            B[I] -= tmp * A[I];
         */
        /*
         real T0 = B[K];
         real T1 = B[K+1] - T0 * A[K+1];
         real T2 = B[K+2] - T0 * A[K+2] - T1 * A[lda+K+2];
        */
        const vec2 t0 = loaddup2(pB);  // { T0 T0 }
        vec2 t2 = fnmadd2(t0, load2d(pA+1), loadu2(pB+1));  // { T1 B2-T0*A[2] }
        vec2 t1 = load2d(pA+1+lda);  // using upper value
        t2 = fnmadd2(unpacklo2(setzero2(), t2), t1, t2);  // { T1 T2 }
        storeu2(pB+1, t2);
        t1 = unpacklo2(t2, t2);  // { T1 T1 }
        t2 = unpackhi2(t2, t2);  // { T2 T2 }
        pA += 3;
        pB += 3;
        if ( pB < end )
        {
            if ( ( end - pB ) & 1 )
            {
                vec2 x = fnmadd1(t0, load1d(pA), load1(pB));
                x = fnmadd1(t1, load1d(pA+lda), x);
                x = fnmadd1(t2, load1d(pA+lda*2), x);
                store1(pB, x);
                pA += 1;
                pB += 1;
            }
#if defined(__AVX__)
            if ( ( end - pB ) & 2 )
            {
                vec2 x = sub2(loadu2(pB), mul2(t0, load2d(pA)));
                vec2 y = add2(mul2(t1, load2d(pA+lda)), mul2(t2, load2d(pA+lda*2)));
                storeu2(pB, sub2(x, y));
                pA += 2;
                pB += 2;
            }
            vec4 tt0 = duplo2f128(cast4(t0));
            vec4 tt1 = duplo2f128(cast4(t1));
            vec4 tt2 = duplo2f128(cast4(t2));
            while ( pB < end-3 )
            {
                vec4 aa = fnmadd4(tt0, load4d(pA), loadu4(pB));
                aa = fnmadd4(tt1, load4d(pA+lda), aa);
                aa = fnmadd4(tt2, load4d(pA+lda*2), aa);
                storeu4(pB, aa);
                pA += 4;
                pB += 4;
            }
#else
            # pragma ivdep
            while ( pB < end )
            {
                vec2 x = fnmadd2(t0, load2d(pA), loadu2(pB));
                x = fnmadd2(t1, load2d(pA+lda), x);
                x = fnmadd2(t2, load2d(pA+lda*2), x);
                storeu2(pB, x);
                pA += 2;
                pB += 2;
            }
#endif
        }
        pA += lda*2;
    }
}
#endif


#if ( DIM < 3 ) && defined(__SSE3__)
/// this version works for any M, but there is another version for DIM==3 below
void alsatian_xtrsmLUN1I_SSE(const int M, const float* A, const int lda, double* B)
{
    assert_true( M <= lda );
    assert_true( M >= 2*DIM ); // 2D:4 3D:6
    A += M * lda;
    int K = M - 1;
#if ( DIM & 1 )
    if ( M & 1 )
    {
        A -= lda;
        double * pB = B;
        double * end = B + K;
        float const* pA = A;
        /*
         const real tmp = B[K] * A[K];
         B[K] = tmp;
         for ( int I = 0; I < K; ++I )
             B[I] -= tmp * A[I];
         */
        vec2 t = mul1(load1(end), load1d(pA+K));
        t = unpacklo2(t, t); // { T0, T0 }
        // there is an even number of scalars remaining, fitting perfectly:
        # pragma ivdep
        while ( pB < end )
        {
            storeu2(pB, fnmadd2(t, load2d(pA), loadu2(pB)));
            pA += 2;
            pB += 2;
        }
        store1(end, t);
        --K;
    }
#endif
    --K;
    //process columns 2 by 2
    for ( ; K >= 0; K -= 2 )
    {
        A -= 2*lda;
        /*
         const real tmp = B[K] * A[K];
         B[K] = tmp;
         for ( int I = 0; I < K; ++I )
             B[I] -= tmp * A[I];
         */
#if 0
        real T1 = A[lda+K+1] * B[K+1];
        real T0 = A[    K  ] * ( B[K  ] - T1 * A[lda+K] );
        B[K  ] = T0;
        B[K+1] = T1;
        for ( int I = 0; I < K; ++I )
            B[I] -= T0 * pA[I] + T1 * pA[I+lda];
        continue;
#endif
        double * pB = B;
        double * end = B + K;
        float const* pA = A;
        vec2 t0 = loadu2(end);
        vec2 a = load2d(pA+lda+K);
        vec2 t1 = mul2(a, t0); // { --, T1 }
        t1 = unpackhi2(t1, t1); // { T1, T1 }
        t0 = mul2(fnmadd2(t1, a, t0), load2d(pA+K)); // { T0, ? }
        storeu2(end, blend11(t0, t1));
        t0 = unpacklo2(t0, t0); // { T0, T0 }
        # pragma ivdep
        while ( pB < end )
        {
            vec2 x = fnmadd2(t0, load2d(pA), loadu2(pB));
            x = fnmadd2(t1, load2d(pA+lda), x);
            storeu2(pB, x);
            pA += 2;
            pB += 2;
        }
    }
}
#endif

#if ( DIM == 3 ) && defined(__SSE3__)
/// this version works for M multiple of 3
void alsatian_xtrsmLUN1I_SSE(const int M, const float* A, const int lda, double* B)
{
    assert_true( M <= lda );
    assert_true( M >= 2*DIM ); // 2D:4 3D:6
    assert_true( M%3 == 0 );
    A += M * lda;
    //process columns 3 by 3
    for ( int K = M - 3; K >= 0; K -= 3 )
    {
        A -= 3*lda;
        /*
         const real tmp = B[K] * A[K];
         B[K] = tmp;
         for ( int I = 0; I < K; ++I )
         B[I] -= tmp * A[I];
         */
#if 0
        real T2 = A[lda*2+K+2] * B[K+2];
        real T1 = A[lda  +K+1] * ( B[K+1] - T2 * A[lda*2+K+1] );
        real T0 = A[      K  ] * ( B[K  ] - T2 * A[lda*2+K  ] - T1 * A[lda+K] );
#endif
#if 0
        B[K  ] = T0;
        B[K+1] = T1;
        B[K+2] = T2;
        for ( int I = 0; I < K; ++I )
            B[I] -= T0 * A[I] + T1 * A[I+lda] + T2 * A[I+lda*2];
        continue;
#endif
        double * pB = B + K;
        float const* pA = A + K;
        vec2 t0, t1, t2;
        {
            vec4f f = loadu4f(pA+2*lda); // last value not used
            vec2 a0 = _mm_cvtps_pd(f);
            vec2 a1 = _mm_cvtps_pd(broadcastZf(f));
            t2 = mul2(loaddup2(pB+2), a1); // { T2, T2 }
            a1 = load2d(pA+lda);
            t0 = fnmadd2(t2, a0, loadu2(pB)); // { preT0, T1/A }
            a0 = load1d(pA);
            t1 = duphi2(mul2(a1, t0)); // { T1, T1 }
            t0 = fnmadd1(t1, a1, t0); // { T0/A, - }
            t0 = duplo2(mul1(a0, t0)); // { T0, T0 }
        }
        storeu2(pB+1, blend11(t1, t2));
        store1(pB, t0);
        if ( pB > B )
        {
            // there could be a odd number of scalar remaining:
            if ( ( pB - B ) & 1 )
            {
                --pA;
                --pB;
                vec2 x = fnmadd1(t0, load1d(pA), load1(pB));
                x = fnmadd1(t1, load1d(pA+lda), x);
                x = fnmadd1(t2, load1d(pA+lda*2), x);
                store1(pB, x);
            }
#if defined(__AVX__)
            if ( ( pB - B ) & 2 )
            {
                pA -= 2;
                pB -= 2;
                vec2 x = fnmadd2(t0, load2d(pA), loadu2(pB));
                x = fnmadd2(t1, load2d(pA+lda), x);
                x = fnmadd2(t2, load2d(pA+lda*2), x);
                storeu2(pB, x);
            }
            vec4 tt0 = duplo2f128(cast4(t0));
            vec4 tt1 = duplo2f128(cast4(t1));
            vec4 tt2 = duplo2f128(cast4(t2));
            while ( pB > B+3 )
            {
                pA -= 4;
                pB -= 4;
                vec4 aa = fnmadd4(tt0, load4d(pA), loadu4(pB));
                aa = fnmadd4(tt1, load4d(pA+lda), aa);
                aa = fnmadd4(tt2, load4d(pA+lda*2), aa);
                storeu4(pB, aa);
            }
#else
            while ( pB > B )
            {
                pA -= 2;
                pB -= 2;
                vec2 x = fnmadd2(t0, load2d(pA), loadu2(pB));
                x = fnmadd2(t1, load2d(pA+lda), x);
                x = fnmadd2(t2, load2d(pA+lda*2), x);
                storeu2(pB, x);
            }
#endif
        }
    }
}
#endif

//------------------------------------------------------------------------------
#pragma mark - SIMD ALSATIAN DTRSM for single-precision matrix argument

#if ( DIM < 3 ) && defined(__SSE3__)
/// this version works for any M, but there is another version for DIM==3 below
void alsatian_xtrsmLLN1U_SSE(const int M, const float* A, const int lda, float* B)
{
    assert_true( M <= lda );
    assert_true( M >= 2*DIM ); // 2D:4
    int K = 0;
#if ( DIM & 1 )
    // process one column to leave an even number
    if ( M & 1 )
    {
        vec4f t = broadcast1f(B);
        float * pB = B + 1;
        float * end = B + M;
        float const* pA = A + 1;
        // there is an even number of scalars, fitting perfectly:
        # pragma ivdep
        while ( pB < end )
        {
            store2f(pB, fnmadd4f(t, load2f(pA), load2f(pB)));
            pA += 2;
            pB += 2;
        }
        A += lda;
        ++K;
    }
#endif
    // process columns 2 by 2
    for ( ; K < M; K += 2 )
    {
        float * pB = B + K;
        float * end = B + M;
        float const* pA = A + K;
        vec4f t0 = load2f(pB);
        vec4f t1 = load2f(pA); // will use upper value
        t1 = fnmadd4f(unpacklo4f(setzero4f(), t0), t1, t0);
        store2f(pB, t1);
        t0 = duplo4f(t0);
        t1 = duphi4f(t1);
        pB += 2;
        pA += 2;
        # pragma ivdep
        while ( pB < end )
        {
            vec4f x = fnmadd4f(t0, load2f(pA), load2f(pB));
            x = fnmadd4f(t1, load2f(pA+lda), x);
            store2f(pB, x);
            pA += 2;
            pB += 2;
        }
        A += 2*lda;
    }
}
#endif

#if ( DIM == 3 ) && defined(__SSE3__)
/// version for M = multiple of 3
void alsatian_xtrsmLLN1U_SSE(const int M, const float* pA, const int lda, float* B)
{
    assert_true( M <= lda );
    assert_true( M >= 2*DIM ); // 3D:6
    assert_true( M%3 == 0 );
    //process columns 3 by 3
    for ( int K = 0; K < M; K += 3 )
    {
        float * pB = B + K;
        float * end = B + M;
        pA += K;
        const vec4f t0 = broadcast1f(pB);  // { T0 T0 }
        vec4f t2 = fnmadd4f(t0, load2f(pA+1), load2f(pB+1));  // { T1 B2-T0*A[2] }
        vec4f t1 = load2f(pA+1+lda);  // using upper value
        t2 = fnmadd4f(unpacklo4f(setzero4f(), t2), t1, t2);  // { T1 T2 }
        store2f(pB+1, t2);
        t1 = duplo4f(t2);  // { T1 T1 }
        t2 = duphi4f(t2);  // { T2 T2 }
        pA += 3;
        pB += 3;
        if ( pB < end )
        {
            if ( ( end - pB ) & 1 )
            {
                vec4f x = fnmadd4f(t0, load1f(pA), load1f(pB));
                x = fnmadd4f(t1, load1f(pA+lda), x);
                x = fnmadd4f(t2, load1f(pA+lda*2), x);
                store1f(pB, x);
                pA += 1;
                pB += 1;
            }
            # pragma ivdep
            while ( pB < end )
            {
                vec4f x = fnmadd4f(t0, load2f(pA), load2f(pB));
                x = fnmadd4f(t1, load2f(pA+lda), x);
                x = fnmadd4f(t2, load2f(pA+lda*2), x);
                store2f(pB, x);
                pA += 2;
                pB += 2;
            }
        }
        pA += lda*2;
    }
}
#endif


#if ( DIM < 3 ) && defined(__SSE3__)
/// this version works for any M, but there is another version for DIM==3 below
void alsatian_xtrsmLUN1I_SSE(const int M, const float* A, const int lda, float* B)
{
    assert_true( M <= lda );
    assert_true( M >= 2*DIM ); // 2D:4 3D:6
    A += M * lda;
    int K = M - 1;
#if ( DIM & 1 )
    if ( M & 1 )
    {
        A -= lda;
        float * pB = B;
        float * end = B + K;
        float const* pA = A;
        vec4f t0 = mul4f(load1f(end), load1f(pA+K));
        t0 = unpacklo4f(t0, t0); // { T0, T0 }
        // there is an even number of scalars remaining, fitting perfectly:
        # pragma ivdep
        while ( pB < end )
        {
            store2f(pB, fnmadd4f(t0, load2f(pA), load2f(pB)));
            pA += 2;
            pB += 2;
        }
        store1f(end, t0);
        --K;
    }
#endif
    --K;
    //process columns 2 by 2
    for ( ; K >= 0; K -= 2 )
    {
        A -= 2*lda;
        float * pB = B;
        float * end = B + K;
        float const* pA = A;
        vec4f t0 = load2f(end);
        vec4f z = load2f(pA+lda+K);
        vec4f t1 = mul4f(z, t0); // { --, T1 }
        t1 = duphi4f(t1); // { T1, T1 }
        t0 = mul4f(fnmadd4f(t1, z, t0), load2f(pA+K)); // { T0, ? }
        store2f(end, blend13f(t0, t1));
        t0 = duplo4f(t0); // { T0, T0 }
        # pragma ivdep
        while ( pB < end )
        {
            vec4f x = fnmadd4f(t0, load2f(pA), load2f(pB));
            x = fnmadd4f(t1, load2f(pA+lda), x);
            store2f(pB, x);
            pA += 2;
            pB += 2;
        }
    }
}
#endif

#if ( DIM == 3 ) && defined(__SSE3__)
/// this version works for M multiple of 3
void alsatian_xtrsmLUN1I_SSE(const int M, const float* A, const int lda, float* B)
{
    assert_true( M <= lda );
    assert_true( M >= 2*DIM ); // 2D:4 3D:6
    assert_true( M%3 == 0 );
    A += M * lda;
    //process columns 3 by 3
    for ( int K = M - 3; K >= 0; K -= 3 )
    {
        A -= 3*lda;
        float * pB = B + K;
        float const* pA = A + K;
        vec4f t2 = mul4f(broadcast1f(pB+2), duplo4f(load1f(pA+2*lda+2))); // { T2, T2 }
        vec4f t0 = fnmadd4f(t2, load2f(pA+2*lda), load2f(pB)); // { -, T1/A }
        vec4f b = load2f(pA+lda);
        vec4f t1 = duphi4f(mul4f(t0, b)); // { T1, T1 }
        t0 = duplo4f(mul4f(fnmadd4f(t1, b, t0), load2f(pA))); // { T0, T0 }
        store2f(pB, blend13f(t0, t1));
        store1f(pB+2, t2);
        if ( pB > B )
        {
            if ( ( pB - B ) & 1 )
            {
                --pA;
                --pB;
                vec4f x = fnmadd4f(t0, load1f(pA), load1f(pB));
                x = fnmadd4f(t1, load1f(pA+lda), x);
                x = fnmadd4f(t2, load1f(pA+lda*2), x);
                store1f(pB, x);
            }
            # pragma ivdep
            while ( pB > B )
            {
                pA -= 2;
                pB -= 2;
                vec4f x = fnmadd4f(t0, load2f(pA), load2f(pB));
                x = fnmadd4f(t1, load2f(pA+lda), x);
                x = fnmadd4f(t2, load2f(pA+lda*2), x);
                store2f(pB, x);
            }
        }
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
    iso_xtrsmLLN<ORD,'N'>(N, A, LDA, B);
    // Solve U*X = B, overwriting B with X. ALPHA = 1.0
    iso_xtrsmLLT<ORD,'N'>(N, A, LDA, B);
#endif
}



/**
 This calls lapack::xpotf2() and then inverts the diagonal terms
*/
void alsatian_xpotf2L(const int N, real* A, const int LDA, int* INFO)
{
    lapack::xpotf2('L', N, A, LDA, INFO);
    if ( 0 == *INFO )
    {
        const int S = LDA+1;
        for ( int u = 0; u < N*S; u += S )
            A[u] = real(1) / A[u];
    }
}


template < int ORD >
void alsatian_xpotrsLref(const int N, const real* A, const int LDA, real* B)
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
    // Solve L*X = B, overwriting B with X. ALPHA = 1.0
    iso_xtrsmLLN<ORD,'I'>(N, A, LDA, B);
    // Solve U*X = B, overwriting B with X. ALPHA = 1.0
    iso_xtrsmLLT<ORD,'I'>(N, A, LDA, B);
#endif
}


void alsatian_xpotrsL(const int N, const real* A, const int LDA, real* B)
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
    if ( INCX == 1 )
    {
        for ( int I = K1; I <= K2; ++I )
        {
            int P = IPIV[I-1];  // P is a one-based array index!
            if ( P != I )
                xswap<ORD>(A+ORD*I, A+ORD*P);
        }
    }
    else if ( INCX > 0 )
    {
        for ( int I = K1; I <= K2; ++I )
        {
            int P = IPIV[INCX*I-1];  // P is a one-based array index!
            if ( P != I )
                xswap<ORD>(A+ORD*I, A+ORD*P);
        }
    }
    else if ( INCX < 0 )
    {
        int IX = 1 + ( 1 - K2 )*INCX;
        for ( int I = K2; I >= K1; --I )
        {
            int P = IPIV[IX-1];  // P is a one-based array index!
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
     but we can shift the pointers to use one-based array indices
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


//------------------------------------------------------------------------------
#pragma mark - LAPACK-STYLE ROUTINES for General LU factorization

/**
 Call lapack::xgetf2() to compute the LU factorization, and inverts diagonal terms
*/
void alsatian_xgetf2(const int N, real* A, const int LDA, int* IPIV, int* INFO)
{
    lapack::xgetf2(N, N, A, LDA, IPIV, INFO);
    if ( 0 == *INFO )
    {
        const int S = LDA+1;
        for ( int u = 0; u < N*S; u += S )
            A[u] = real(1) / A[u];
    }
}

/// factorize a single-precision matrix
void alsatian_sgetf2(const int N, float* A, const int LDA, int* IPIV, int* INFO)
{
    lapack::sgetf2(N, N, A, LDA, IPIV, INFO);
    if ( 0 == *INFO )
    {
        const int S = LDA+1;
        for ( int u = 0; u < N*S; u += S )
            A[u] = float(1) / A[u];
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


/// following the standard LAPACK interface
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

/// version of xgetrs('N', ...) for NRHS==1
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


/// version of xgetrs('N', ...) for NRHS==1 and inverted diagonal terms
/**
 This is used to apply the full block preconditionner and we could use
 Non-temporal loads for the matrix A, since it will not fit in the cache,
 but the vector B should be cached!
 */
void alsatian_xgetrsN(int N, const real* A, int LDA, const int* IPIV, real* B)
{
    // Apply row interchanges to the right hand side.
    xlaswp1(B, 1, N, IPIV); //iso_xlaswp<1>(B, 1, N, IPIV, 1);
    // Solve L*X = B, overwriting B with X.
    alsatian_xtrsmLLN1U(N, (float*)A, LDA, B);
    // Solve U*X = B, overwriting B with X.
    alsatian_xtrsmLUN1I(N, (float*)A, LDA, B);
}


#if defined(__SSE3__)
void alsatian_xgetrsN_SSE(int N, const real* A, int LDA, const int* IPIV, real* B)
{
#if 0
    // checking result against other routine:
    real * T = new_real(N); copy_real(N, B, T);
#endif
    // Apply row interchanges to the right hand side.
    xlaswp1(B, 1, N, IPIV);
    // Solve L*X = B, overwriting B with X.
    alsatian_xtrsmLLN1U_SSE(N, (float*)A, LDA, B);
    // Solve U*X = B, overwriting B with X.
    alsatian_xtrsmLUN1I_SSE(N, (float*)A, LDA, B);
#if 0
    // checking result against other routine:
    alsatian_xgetrsN(N, A, LDA, IPIV, T);
    real err = blas::difference(N, B, T);
    if ( std::isnan(err) || err > 0.01 )
    {
        int S = std::min(N/(2*DIM), 2)*DIM, E = N - 1 - S;
        printf("\n xgetrs %3i ",N); VecPrint::print(S, B, 3); VecPrint::print(S, B+E, 3);
        printf("\n    ref --- "); VecPrint::print(S, T, 3); VecPrint::print(S, T+E, 3);
    }// else printf("\n xgetrs %3i okay!", N);
    free_real(T);
#endif
}
#endif

/// version of xgetrs('N', ...) for ORD interleaved vectors in right-hand-side B
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


/// version of xgetrs('N', ...) for ORD interleaved vectors in right-hand-side B
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
    iso_xtrsmLLN<ORD,'U'>(N, A, LDA, B);
    iso_xtrsmLLT<ORD,'N'>(N, A, LDA, B);
#endif
}

#endif
