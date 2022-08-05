// Cytosim was created by Francois Nedelec.
// Copyright 2019-2020 Sainsbury Laboratory, Cambridge University

#ifndef XPTTRF_H
#define XPTTRF_H

/// returns 1/x
static inline real inverse(real x) { return real(1) / x; }


/**
 This is a C-translation of the LAPACK reference implementation of dpttrf()
 to factorize a symmetric definite-positive tridiagonal matrix.
 The method is known as Thomas' algorithm.
*/
void lapack_xpttrf_ori(int size, real* D, real* E, int* INFO)
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

/// modified version only loads D[] at one place
void lapack_xpttrf(int size, real* D, real* E, int* INFO)
{
    *INFO = 0;
    real d = D[0];
    for ( int i = 0; i < size-1; ++i )
    {
        if ( D[i] < 0 )
        {
            *INFO = i;
            return;
        }
        real e = E[i] / d;
        d = D[i+1] - e * E[i];
        D[i+1] = d;
        E[i] = e;
    }
}

/**
 This is a C-translation of the LAPACK reference implementation of dpttrf()
 to factorize a symmetric definite-positive tridiagonal matrix.
 The method is known as Thomas' algorithm. There is no test for positivity.
*/
void xpttrf(int size, real* D, real* E)
{
    real d = D[0];
    for ( int i = 0; i < size-1; ++i )
    {
        real e = E[i] / d;
        d = D[i+1] - e * E[i];
        D[i+1] = d;
        E[i] = e;
    }
}

/**
 This is a C-translation of the LAPACK reference implementation of dptts2()
 
 *     Solve A * X = B using the factorization A = L*D*L**T,
 *     overwriting each right hand side vector with its solution.
 
     DO I = 2, N
         B(I) = B(I) - B(I-1) * E(I-1)
     CONTINUE
 
     B(N) = B(N) / D(N)
 
     DO I = N - 1, 1, -1
         B(I) = B(I) / D(I) - B(I+1) * E(I)
     CONTINUE
 */
void lapack_xptts2(int size, const real* D, const real* E, real* B)
{
    assert_true(size > 0);

    for ( int i = 1; i < size; ++i )
        B[i] = B[i] - B[i-1] * E[i-1];
    
    B[size-1] = B[size-1] / D[size-1];
    
    for ( int i = size-2; i >= 0; --i )
        B[i] = B[i] / D[i] - B[i+1] * E[i];
}


inline void lapack_xptts2(int size, int nrhs, real const* D, real const* E, real* B, int LDB)
{
    assert_true( nrhs == 1 ); // in this case, the last argument (LDB) is ignored
    return lapack_xptts2(size, D, E, B);
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 Translated from pseudo code to C
 https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
     A is the sub diagonal of size N-1
     B is the diagonal of size N
     C is the upper diagonal of size N-1
 B is modified in this version!
 ATTENTION: Wikipedia version has A starting at index 2.
 This code has all arrays starting at index 0.
 */
void wikipedia_solve(int N, real const* A, real* B, real const* C, real* X)
{
    for ( int i = 1; i < N; ++i )
    {
        real W = A[i-1] / B[i-1];
        B[i] = B[i] - W * C[i-1];
        X[i] = X[i] - W * X[i-1];
    }
    X[N-1] = X[N-1] / B[N-1];
    for ( int i = N-2; i >= 0; --i )
        X[i] = ( X[i] - C[i] * X[i+1] ) / B[i];
}

/**
 From Numerical Recipee, 3rd Ed, page 56
 2.4 Tridiagonal and Band-Diagonal Systems of Equations

 Solves for a vector u[0..n-1] the tridiagonal linear set given by equation (2.4.1).
 a[0..n-2], b[0..n-1], c[0..n-2], and r[0..n-1] are input vectors and are not modified.
 */
void nr_tridag(const int size, const real a[], real b[], const real c[], const real r[], real u[])
{
    real * gam = new real[size];
    real bet = b[0];
    if (bet == 0.0) throw("Error 1 in tridag");
    u[0] = r[0] / bet;
    for ( int j = 1; j < size; ++j )
    {
        gam[j] = c[j-1] / bet;
        bet = b[j] - a[j-1] * gam[j];
        if (bet == 0.0) throw("Error 2 in tridag");
        u[j] = ( r[j] - a[j-1] * u[j-1] ) / bet;
    }
    for ( int j = size-2; j >= 0; --j )
        u[j] -= gam[j+1] * u[j+1];
    delete[] gam;
}

/**
 From Numerical Recipee, 3rd Ed, page 56
 2.4 Tridiagonal and Band-Diagonal Systems of Equations

 Modified FJN 4.8.2022 to follow LAPACK's interface.
 Solves for a vector u[0..n-1] the tridiagonal linear set given by equation (2.4.1).
 a[1..n-1] and c[0..n-2] are input vectors and are not modified.
 b[0..n-1], is input and modified.
 x[0..n-1] is both input (right-hand side) and output (result)
 */

void nr_tridag(const int size, const real a[], real b[], const real c[], real x[])
{
    real bet = b[0];
    //if (bet == 0.0) throw("Error 1 in tridag");
    x[0] = x[0] / bet;
    for ( int j = 1; j < size; ++j )
    {
        real gam = c[j-1] / bet;
        bet = b[j] - a[j-1] * gam;
        b[j] = gam;
        //if (bet == 0.0) throw("Error 2 in tridag");
        x[j] = ( x[j] - a[j-1] * x[j-1] ) / bet;
    }
    for ( int j = size-2; j >= 0; --j )
        x[j] -= b[j+1] * x[j+1];
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
    D[0] = inverse(D[0]);
    
    for ( int n = 1; n < size; ++n )
        D[n] = inverse( D[n] - ( E[n-1] * E[n-1] ) * D[n-1] );
}


/**
 Italian solve without divisions
 
 PSEUDO code with indices starting at 1:
     y(1) = gamma(1) * f(1);
     for i = 2:N
         y(i) = gamma(i) * ( f(i) - b(i) * y(i-1) );
     x(N) = y(N);
     for i = N-1:-1:1
         x(i) = y(i) - gamma(i) * c(i) * x(i+1);
 */
void italian_xptts2(int size, real const* D, real const* E, real* B)
{
    assert_true(size > 0);
 
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


inline void italian_xptts2(int size, int nrhs, real const* D, real const* E, real* B, int LDB)
{
    assert_true(nrhs == 1); // in this case, the last argument (LDB) is ignored
    return italian_xptts2(size, D, E, B);
}

/**
 This implements Thomas's algorithm to solve a linear system
 with a tridiagonal matrix {L, D, U} and right-hand side vector 'B'
 
     L is lower diagonal, with valid index in [0, size-2]
     D is diagonal, with valid index in [0, size-1]
     U is upper diagonal, with valid index in [0, size-2]
     B is input/output vector, with valid index in [0, size-1]
 
 This works even if 'L == U'
 
 WARNING: This code was not tested for L != U
 Modified from:
 https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm
 */
void italian_solve(int size, real* L, real* D, real* U, real* X)
{
    D[0] = inverse(D[0]);
    X[0] = D[0] * X[0];
    if ( size > 1 )
    {
        // upward recursion
        for ( int n = 1; n < size; ++n )
        {
            D[n] = inverse( D[n] - ( L[n-1] * U[n-1] ) * D[n-1] );
            X[n] = D[n] * ( X[n] - U[n-1] * X[n-1] );
        }
        // downward recursion
        for ( int i = size-2; i > 0; --i )
            X[i] = X[i] - ( D[i] * U[i] ) * X[i+1];
        X[0] = X[0] - ( D[0] * U[0] ) * X[1];
    }
}


void italian2_solve(int size, real* L, real* D, real* U, real* X)
{
    real d = inverse(D[0]);
    real e = U[0] * L[0];
    real l = L[0];
    D[0] = U[0] * d;
    X[0] = X[0] * d;
    if ( size > 1 )
    {
        // upward recursion
        for ( int n = 1; n < size; ++n )
        {
            //D[n] = inverse( D[n] - ( L[n-1] * U[n-1] ) * D[n-1] );
            d = inverse( D[n] - e * d );
            //X[n] = D[n] * ( X[n] - U[n-1] * X[n-1] );
            X[n] = ( X[n] - l * X[n-1] ) * d;
            D[n] = U[n] * d;
            l = L[n];
            e = l * U[n];
        }
        // downward recursion
        for ( int n = size-2; n > 0; --n )
            X[n] = X[n] - D[n] * X[n+1];
        X[0] = X[0] - D[0] * X[1];
    }
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 Based on the 'Italian' version, with an additional substitution:
 
     E'[n] <-  D[n] * E[n]

 D[] is addressed within [0, size-1]
 E[] is addressed within [0, size-2]
 Improved on 02.04.2020: reduced memory use
 */
void alsatian_xpttrf(const int size, real* D, real* E)
{
    if ( size < 1 || D[0] < 0 )
        return;
    real w = inverse(D[0]);
    D[0] = w;
    if ( size > 1 )
    {
        real e = E[0] * E[0];
        E[0] = w * E[0];
        const int end = size-1;
        for ( int n = 1; n < end; ++n )
        {
            //if ( D[n] < 0 ) return n;
            w = inverse( D[n] - w * e );
            e = E[n] * E[n];
            D[n] = w;
            E[n] = w * E[n];
        }
        //if ( D[end] < 0 ) return end;
        D[end] = inverse( D[end] - w * e );
    }
}


void alsatian_xpttrf(int size, real* D, real* E, int* INFO)
{
    if ( size < 1 || D[0] < 0 )
    {
        *INFO = 1;
        return;
    }
    real w = inverse(D[0]);
    D[0] = w;
    if ( size > 1 )
    {
        real e = E[0] * E[0];
        E[0] = w * E[0];
        const int end = size-1;
        for ( int n = 1; n < end; ++n )
        {
            if ( D[n] < 0 )
            {
                *INFO = n+1;
                return;
            }
            w = inverse( D[n] - w * e );
            e = E[n] * E[n];
            D[n] = w;
            E[n] = w * E[n];
        }
        if ( D[end] < 0 )
        {
            *INFO = end+1;
            return;
        }
        D[end] = inverse( D[end] - w * e );
    }
    *INFO = 0;
}

/**
 Based on the 'Italian' version, using precalculated constant terms
 
 Improved on 02.04.2020: removed one multiplication, now only using FMAs
 */
void alsatian2_xptts2(int size, real const* D, real const* DE, real* B)
{
    assert_true(size > 0);
    real x = B[0];
    B[0] = x * D[0];
    // upward recursion on B[]
    #pragma nounroll
    for ( int n = 1; n < size; ++n )
    {
        x = B[n] - x * DE[n-1];
        B[n] = x * D[n];
    }
    // downward recursion on B[]
    x = B[size-1];
    #pragma nounroll
    for ( int n = size-2; n >= 0; --n )
    {
        x = B[n] - x * DE[n];
        B[n] = x;
    }
}

/**
 Based on the 'Italian' version, using precalculated constant terms
 
 Improved on 02.04.2020: removed one multiplication, now only using FMAs
 */
void alsatian_xptts2(int size, real const* D, real const* DE, real* B)
{
    assert_true(size > 0);
    real x = B[0];
    B[0] = x * D[0];
    // upward recursion on B[]
    #pragma nounroll
    for ( int n = 1; n < size; ++n )
    {
        x = B[n] - x * DE[n-1];
        B[n] = x * D[n];
    }
    // downward recursion on B[]
    x = B[size-1];
    #pragma nounroll
    for ( int n = size-2; n >= 0; --n )
    {
        x = B[n] - x * DE[n];
        B[n] = x;
    }
}


/// offering LAPACK's standard arguments
inline void alsatian_xptts2(int size, int nrhs, real const* D, real const* DE, real* B, int LDB)
{
    assert_true(nrhs == 1); // in this case, the last argument is ignored
    return alsatian_xptts2(size, D, DE, B);
}

/**
This implements Thomas's algorithm to solve a linear system
with a tridiagonal symmetric matrix {E, D, E} and right-hand side 'B'
*/
void alsatian_solve(int size, real* D, real* E, real* B)
{
    assert_true(size > 0);
    real w = inverse(D[0]);
    real y = B[0];
    real e = E[0] * E[0];
    real x = w * E[0];
    D[0] = w;
    E[0] = x;
    B[0] = y * w;

    // upward recursion
    for ( int n = 1; n < size; ++n )
    {
        w = inverse( D[n] - w * e );
        e = E[n] * E[n];
        y = B[n] - y * x;
        D[n] = w;
        x = w * E[n];
        E[n] = x;
        B[n] = y * w;
    }

    y = B[size-1];
    // downward recursion on B[]
    for ( int n = size-2; n >= 0; --n )
    {
        y = B[n] - E[n] * y;
        B[n] = y;
    }
}

#endif
