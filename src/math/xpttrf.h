// Cytosim was created by Francois Nedelec.
// Copyright 2019-2020 Sainsbury Laboratory, Cambridge University

#ifndef XPTTRF_H
#define XPTTRF_H

/**
 This is a C-translation of the LAPACK reference implementation of dpttrf()
 to factorize a symmetric definite-positive tridiagonal matrix.
 The method is known as Thomas' algorithm.
*/
void lapack_xpttrf(int size, real* D, real* E, int* INFO)
{
    *INFO = 0;
    for ( int i = 0; i < size-1; ++i )
    {
        if ( D[i] < 0 )
        {
            *INFO = i;
            return;
        }
        real e = E[i];
        E[i] = e / D[i];
        D[i+1] = D[i+1] - e * E[i];
    }
}

/**
 This is a C-translation of the LAPACK reference implementation of dptts2()
 
 *     Solve A * X = B using the factorization A = L*D*L**T,
 *     overwriting each right hand side vector with its solution.
 
     DO I = 2, N
         B( I ) = B( I ) - B( I-1 ) * E( I-1 )
     CONTINUE
 
     B( N ) = B( N ) / D( N )
 
     DO I = N - 1, 1, -1
         B( I ) = B( I ) / D( I ) - B( I+1 ) * E( I )
     CONTINUE
 */
void lapack_xptts2(int size, int NRHS, const real* D, const real* E, real* B, int LDB)
{
    assert_true( NRHS == 1 ); // in this case, LDB is not used

    for ( int i = 1; i < size; ++i )
        B[i] = B[i] - B[i-1] * E[i-1];
    
    B[size-1] = B[size-1] / D[size-1];
    
    for ( int i = size-2; i >= 0; --i )
        B[i] = B[i] / D[i] - B[i+1] * E[i];
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 Italian factorization that uses different operations
 From
     Numerical Mathematics
     Springer (2000) ISBN 0-387-98959-5
     A. Quarteroni, R. Sacco and F. Saleri
     Page 93
 
 PSEUDO code with indices starting at 1:
     a_i = diagonal terms i in [1...N]
     b_i = [0;b]; lower diagonal
     c_i = [c;0]; upper diagonal
     gamma(1) = 1.0 / a(1);
     for i = 2:N
         gamma(i) = 1.0 / ( a(i) - b(i) * gamma(i-1) * c(i-1) );
 */
void italian_xpttrf(int size, real* D, real const* E, int* INFO)
{
    *INFO = 0;
    D[0] = 1.0 / D[0];
    
    for ( int n = 1; n < size; ++n )
        D[n] = 1.0 / ( D[n] - ( E[n-1] * E[n-1] ) * D[n-1] );
}


/**
 Italian version without divisions
 
 PSEUDO code with indices starting at 1:
     y(1) = gamma(1) * f(1);
     for i = 2:N
         y(i) = gamma(i) * ( f(i) - b(i) * y(i-1) );
     x(N) = y(N);
     for i = N-1:-1:1
         x(i) = y(i) - gamma(i) * c(i) * x(i+1);
 */
void italian_xptts2(int size, int nrhs, real const* D, real const* E, real* B, int LDB)
{
    assert_true(size > 0);
    assert_true(nrhs == 1); // in this case, LDB is not used
 
    B[0] = D[0] * B[0];
    if ( size > 1 )
    {
        // upward recursion on B[]
        for ( int n = 1; n < size; ++n )
            B[n] = D[n] * ( B[n] - E[n-1] * B[n-1] );
        
        // downward recursion on B[]
        for ( int n = size-2; n > 0; --n )
            B[n] = B[n] - ( D[n] * E[n] ) * B[n+1];
        B[0] = B[0] - ( D[0] * E[0] ) * B[1];
    }
}


/**
 This implements Thomas's algorithm to solve a linear system
 with a tridiagonal matrix {L, D, U} and right-hand side vector 'B'
 
 L is lower diagonal, with valid index in [0, size-2]
 D is diagonal, with valid index in [0, size-1]
 U is upper diagonal, with valid index in [0, size-2]
 B is input/output vector, with valid index in [0, size-1]
 
 This works even if 'L == U'
 
 This code was not tested for L != U!!!
 Modified from:
 https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm
 */
void italian_thomas(size_t size, real* L, real* D, real* U, real* B)
{
#if 0
    D[0] = 1.0 / D[0];
    B[0] = D[0] * B[0];
    if ( size > 1 )
    {
        // upward recursion on B[]
        for ( size_t n = 1; n < size; ++n )
        {
            D[n] = 1.0 / ( D[n] - ( L[n-1] * U[n-1] ) * D[n-1] );
            B[n] = D[n] * ( B[n] - U[n-1] * B[n-1] );
        }
        // downward recursion on B[]
        for ( size_t i = size-2; i > 0; --i )
            B[i] = B[i] - ( D[i] * U[i] ) * B[i+1];
        B[0] = B[0] - ( D[0] * U[0] ) * B[1];
    }
#else
    real l = L[0];
    real e = U[0] * L[0];
    real d = 1 / D[0];
    U[0] = U[0] / D[0];
    B[0] = B[0] / D[0];
    if ( size > 1 )
    {
        // upward recursion on B[]
        for ( size_t n = 1; n < size; ++n )
        {
            //D[n] = 1.0 / ( D[n] - ( L[n-1] * U[n-1] ) * D[n-1] );
            d = 1.0 / ( D[n] - e * d );
            //B[n] = D[n] * ( B[n] - U[n-1] * B[n-1] );
            B[n] = d * ( B[n] - l * B[n-1] );
            l = L[n];
            e = l * U[n];
            U[n] = U[n] * d;
        }
        // downward recursion on B[]
        for ( size_t n = size-2; n > 0; --n )
            B[n] = B[n] - U[n] * B[n+1];
        B[0] = B[0] - U[0] * B[1];
    }
#endif
}


/*
 Adapted from
 https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
 */
void tridiagonal_solve(size_t N, real* A, real* B, real* C, real* X)
{
#if 1
    // revised version to supress divisions in the downward recursion
    X[0] = X[0] / B[0];
    B[0] = C[0] / B[0];
    for ( size_t i = 1; i < N; ++i )
    {
        real W = 1.0 / ( B[i] - B[i-1] * A[i-1] );
        X[i] = W * ( X[i] - A[i-1] * X[i-1] );
        B[i] = W * C[i];
    }
    for ( size_t i = N-2; i > 0; --i )
        X[i] = X[i] - B[i] * X[i+1];
    X[0] = X[0] - B[0] * X[1];
#else
    // wikipedia's version translated to C
    for ( size_t i = 1; i < N; ++i )
    {
        real W = A[i-1] / B[i-1];
        B[i] = B[i] - W * C[i-1];
        X[i] = X[i] - W * X[i-1];
    }
    X[N-1] = X[N-1] / B[N-1];
    for ( size_t i = N-2; i > 0; --i )
        X[i] = (X[i] - C[i] * X[i+1]) / B[i];
    X[0] = (X[0] - C[0] * X[1]) / B[0];
#endif
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 Based on the 'Italian' version, with an additional substitution:
 
     E'[n] <-  D[n] * E[n]

 Improved on 02.04.2020: reduced memory use
 */
void alsatian_xpttrf(size_t size, real* D, real* E, int* INFO)
{
    *INFO = 0;
    if ( size < 1 ) return;

    real x = 0;
    real e = 0;

    for ( size_t n = 0; n < size-1; ++n )
    {
        //D[n] = 1.0 / ( D[n] - E[n-1] * E[n-1] * D[n-1] );
        //x = 1.0 / ( D[n] - ( E[n-1] * E[n-1] ) * x );
        if ( D[n] < 0 )
        {
            *INFO = n;
            return;
        }
        x = 1.0 / ( D[n] - e * x );
        e = E[n];
        D[n] = x;
        x = x * e;
        E[n] = x;
    }
    D[size-1] = 1.0 / ( D[size-1] - e * x );
}


/**
 Based on the 'Italian' version, using precalculated constant terms
 
 Improved on 02.04.2020: removed one multiplication, now only using FMAs
 */
void alsatian_xptts2(size_t size, size_t nrhs, real const* D, real const* DE, real* B, size_t LDB)
{
    assert_true(size > 0);
    assert_true(nrhs == 1); // in this case, LDB is not used

    // upward recursion on B[]
    real x = B[0];
    for ( size_t n = 0; n < size-1; ++n )
    {
        //B[n] = D[n] * ( B[n] - B[n-1] * E[n-1] );
        B[n] = D[n] * x;
        x = B[n+1] - x * DE[n];  // = B[n+1] - B[n] * E[n]
    }
    x = D[size-1] * x;
    B[size-1] = x;
    
    if ( size > 1 )
    {
        // downward recursion on B[]
        for ( size_t n = size-2; n > 0; --n )
        {
            // B[n] = B[n] - ( D[n] * E[n] ) * B[n+1];
            x = B[n] - DE[n] * x;
            B[n] = x;
        }
        B[0] = B[0] - DE[0] * x;
    }
}


/**
This implements Thomas's algorithm to solve a linear system
with a tridiagonal symmetric matrix {E, D, E} and right-hand side 'B'
*/
void alsatian_thomas(size_t size, real* D, real* E, real* B)
{
    real x = 0;
    real e = 0;
    real y = B[0];

    // upward recursion
    for ( size_t n = 0; n < size-1; ++n )
    {
        x = 1.0 / ( D[n] - e * x );
        e = E[n];
        D[n] = x;
        x = x * e;
        E[n] = x;
        B[n] = D[n] * y;
        y = B[n+1] - y * E[n];
    }
    x = 1.0 / ( D[size-1] - e * x );
    D[size-1] = x;
    E[size-1] = x * E[size-1];
    y = D[size-1] * y;
    B[size-1] = y;

    if ( size > 1 )
    {
        // downward recursion on B[]
        for ( size_t n = size-2; n > 0; --n )
        {
            y = B[n] - E[n] * y;
            B[n] = y;
        }
        B[0] = B[0] - E[0] * y;
    }
}

#endif
