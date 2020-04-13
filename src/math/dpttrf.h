// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef DPTTRF_H
#define DPTTRF_H

/**
 This is a C-translation of the LAPACK reference implementation of dpttrf()
 (the method is known as Thomas' algorithm to factorize a tridiagonal matrix)
*/
void lapack_xpttrf(int N, real* D, real* E, int* INFO)
{
    for ( int i = 0; i < N-1; ++i )
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
void lapack_xptts2(int N, int NRHS, const real* D, const real* E, real* B, int LDB)
{
    assert_true( NRHS == 1 ); // in this case, LDB is not used

    for ( int i = 1; i < N; ++i )
        B[i] = B[i] - B[i-1] * E[i-1];
    
    B[N-1] = B[N-1] / D[N-1];
    
    for ( int i = N-2; i >= 0; --i )
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
     b=[0;b];
     c=[c;0];
     gamma(1) = 1.0 / a(1);
     for i = 2:N
         gamma(i) = 1.0 / ( a(i) - b(i) * gamma(i-1) * c(i-1) );
 */
void italian_xpttrf(int size, real* D, real const* E, int* INFO)
{
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
    
    // upward recursion on B[]
    for ( int n = 1; n < size; ++n )
        B[n] = D[n] * ( B[n] - E[n-1] * B[n-1] );
    
    // downward recursion on B[]
    if ( size > 1 )
    {
        for ( int n = size-2; n > 0; --n )
            B[n] = B[n] - ( D[n] * E[n] ) * B[n+1];
        B[0] = B[0] - ( D[0] * E[0] ) * B[1];
    }
}


/**
 This implements Thomas's algorithm to solve a linear system
 with a tridiagonal matrix {L, D, U} and right-hand side 'B'
 
 L is lower diagonal, with valid index in [0, size-2]
 D is diagonal, with valid index in [0, size-1]
 U is upper diagonal, with valid index in [0, size-2]
 B is input/output, with valid index in [0, size-1]
 
 This works even if 'L == U'
 
 Modified from:
 https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm
 */
void italian_thomas(size_t size, real const*L, real const* D, real* U, real* B)
{
    real e = L[0];
    U[0] = U[0] / D[0];
    B[0] = B[0] / D[0];
    
    // upward recursion: U is changed, D and L are not changed
    for ( size_t i = 1; i < size; ++i )
    {
        //const real m = 1.0 / (D[i] - L[i-1] * U[i-1]);
        const real m = 1.0 / ( D[i] - e * U[i-1] );
        //B[i] = (B[i] - L[i-1] * B[i-1]) * m;
        B[i] = ( B[i] - e * B[i-1] ) * m;
        e = L[i];
        U[i] = U[i] * m;
    }
    
    // downward recursion: D[] and L[] are not used
    if ( size > 1 )
    {
        for ( size_t i = size-2; i > 0; --i )
            B[i] = B[i] - U[i] * B[i+1];
        B[0] = B[0] - U[0] * B[1];
    }
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 Based on the 'Italian' version, with an additional substitution:
 
     E'[n] <-  D[n] * E[n]

 */
void alsatian_xpttrf(size_t size, real* D, real* E, int* INFO)
{
    real x = 1.0 / D[0];
    real e = E[0];
    D[0] = x;
    x = x * e;
    E[0] = x;

    for ( size_t n = 1; n < size; ++n )
    {
        //D[n] = 1.0 / ( D[n] - E[n-1] * E[n-1] * D[n-1] );
        //x = 1.0 / ( D[n] - ( E[n-1] * E[n-1] ) * x );
        x = 1.0 / ( D[n] - e * x );
        e = E[n];
        D[n] = x;
        x = x * e;
        E[n] = x;
    }
}


/**
 Based on the 'Italian' version, using precalculated constant terms
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
    y = y * D[size-1];
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
