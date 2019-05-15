// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef DPTTRF_H
#define DPTTRF_H

/**
 This is a C-translation of the LAPACK reference implementation of dpttrf()
 (that is Thomas' algorithm to factorize a tridiagonal matrix)
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
void italian_xpttrf(int size, real* D, real* E, int* INFO)
{
    D[0] = 1.0 / D[0];
    
    for ( int n = 1; n < size; ++n )
        D[n] = 1.0 / ( D[n] - E[n-1] * E[n-1] * D[n-1] );
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
    assert_true(nrhs == 1); // in this case, LDB is not used
 
    B[0] = D[0] * B[0];
    
    // upward recursion on B[]
    for ( int n = 1; n < size; ++n )
        B[n] = D[n] * ( B[n] - B[n-1] * E[n-1] );
    
    // downward recursion on B[]
    for ( int n = size-2; n >= 0; --n )
        B[n] = B[n] - D[n] * E[n] * B[n+1];
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 Based on the 'Italian' version, precalculating constant products:
 
     DEL[n] <-  D[n] * E[n-1]
     E[n]   <-  D[n] * E[n]
 */
void alsatian_xpttrf(int size, real* D, real* E, int* INFO)
{
    real * DEL = D + size;
    D[0] = 1.0 / D[0];
    
    real x = D[0];
    for ( int n = 1; n < size; ++n )
    {
        //D[n] = 1.0 / ( D[n] - E[n-1] * E[n-1] * D[n-1] );
        x = 1.0 / ( D[n] - E[n-1] * E[n-1] * x );
        D[n] = x;
        // precalculate product that is constant:
        DEL[n] = x * E[n-1];
    }
    
    // precalculate two products that are constant:
    for ( int n = 0; n < size-1; ++n )
        E[n] = D[n] * E[n];
}


/**
 Based on the 'Italian' version, using precalculated constant terms
 */
void alsatian_xptts2(int size, int nrhs, real const* D, real const* E, real* B, int LDB)
{
    assert_true(nrhs == 1); // in this case, LDB is not used
    real const* DEL = D + size;

    B[0] = D[0] * B[0];
    
    // upward recursion on B[]
    for ( int n = 1; n < size; ++n )
    {
        // B[n] = D[n] * ( B[n] - B[n-1] * E[n-1] );
        B[n] = D[n] * B[n] - DEL[n] * B[n-1];
    }
    
    // downward recursion on B[]
    for ( int n = size-2; n >= 0; --n )
    {
        // B[n] = B[n] - D[n] * E[n] * B[n+1];
        B[n] = B[n] - E[n] * B[n+1];
    }
}

#endif
