// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "quartic_solver.h"
#include <cstdio>


#define QUARTIC_VERBOSE 0

    
/// apply one step of Halley's method (convergence is cubic in general)
void QuarticSolver::refineQuadratic(const real A, const real B, const real C, real& x)
{
    real   F = quadratic(A, B, C, x);
    real  dF = 2*A*x + B;
    
    if ( dF != 0 )
        x -= ( F * dF ) / ( dF * dF - A * F );
}


/// apply one step of Halley's method (convergence is cubic in general)
void QuarticSolver::refineCubic(const real A, const real B, const real C, const real D, real& x)
{
    real   F = cubic(A, B, C, D, x);
    real  dF = quadratic(3*A, 2*B, C, x);
    
    if ( dF != 0 )
    {
        real ddF2 = 3*A*x + B;
        x -= ( F * dF ) / ( dF * dF - F * ddF2 );
    }
}


/// apply one step of Halley's method (convergence is cubic in general)
void QuarticSolver::refineQuartic(const real A, const real B, const real C, const real D, const real E, real& x)
{
    real   F = quartic(A, B, C, D, E, x);
    real  dF = cubic(4*A, 3*B, 2*C, D, x);
    
    if ( dF != 0 )
    {
        real ddF2 = quadratic(6*A, 3*B, C, x);
        x -= ( F * dF ) / ( dF * dF - F * ddF2 );
    }
}

//----------------------------------------------------------------------------
#pragma mark Quadratic

/**
 return number of roots found (0, 1, 2).
 if ( a > 0 ) then r1 > r2
 
 In the particular case of a double root, this returns 2
 
 http://en.wikipedia.org/wiki/Quadratic_function
 */
int QuarticSolver::solveQuadratic(const real A, const real B, const real C,
                                  real& r1, real& r2 )
{
    if ( A == 0.0 )
    {
        if ( B == 0.0 )
            return -1;
        r1 = -C/B;
        return 1;
    }
    
    real A2 = A * 2;
    real discriminant = B * B - 4 * A * C;
    
    if ( discriminant >= 0.0 )
    {
        discriminant = std::sqrt(discriminant);
        r1 =  ( discriminant - B ) / A2;
        r2 = -( discriminant + B ) / A2;
        return 2;
    }
    
    return 0;
}
    
    
int QuarticSolver::solveQuadratic(const real A, const real B, const real C,
                                  cplx& r1, cplx& r2 )
{
    if ( A == 0.0 )
    {
        if ( B == 0.0 )
            return -1;
        r1 = cplx(-C/B, 0);
        return 1;
    }
    
    real A2 = A * 2;
    real discriminant = B * B - 4 * A * C;
    
    if ( discriminant >= 0.0 )
    {
        discriminant = std::sqrt(discriminant);
        r1 = cplx( ( discriminant - B ) / A2, 0);
        r2 = cplx(-( discriminant + B ) / A2, 0);
    }
    else
    {
        discriminant = std::sqrt(-discriminant);
        r1 = cplx( -B / A2,  discriminant / A2 );
        r2 = cplx( -B / A2, -discriminant / A2 );
    }
    
    return 2;
}

//----------------------------------------------------------------------------
#pragma mark Cubic

/**
 return number of real roots found (0, 1, 3).
 The roots are returned unsorted.
 
 http://en.wikipedia.org/wiki/Cubic_function
 */
int QuarticSolver::solveCubicUnsorted(real A, real B, real C, real D,
                                      real& r1, real& r2, real& r3)
{
    if ( A == 0.0 )
        return solveQuadratic(B, C, D, r1, r2);
    
    if ( D == 0.0 )
    {
#if ( QUARTIC_VERBOSE > 0 )
        fprintf(stderr, " cubic has null root\n");
#endif
        r1 = 0;
        return 1 + solveQuadratic(A, B, C, r2, r3);
    }
    
    if ( A != 1.0 )
    {
        B /= A;
        C /= A;
        D /= A;
    }
    
#if ( QUARTIC_VERBOSE > 1 )
    fprintf(stderr, "  cubic x^3  %+f x^2  %+f x  %+f\n", B, C, D);
#endif
    
    // calculate determinant of cubic
    real B3 = B / 3.0;
    real Q = B3*B3 - C / 3.0, QQQ = Q*Q*Q;
    real R = B3 * ( B3*B3 - 0.5 * C ) + 0.5 * D, RR = R*R;
    
    if ( R == 0 && Q == 0 )
    {
        // 3 identical real roots
        r1 = -B3;
        r2 = -B3;
        r3 = -B3;
        return 3;
    }
    else if ( RR < QQQ )
    {
        /* This sqrt and division is safe, since RR >= 0, so QQQ > RR,    */
        /* so QQQ > 0.  The acos is also safe, since RR/QQQ < 1, and      */
        /* thus R/std::sqrt(QQQ) < 1.                                          */
        real theta = std::acos(R/std::sqrt(QQQ));
        /* This sqrt is safe, since QQQ >= 0, and thus Q >= 0             */
        
        real x = -2*std::sqrt(Q);
        real t = theta/3.0;
        r1 = x * std::cos(t)            - B3;
        r2 = x * std::cos(t+2*M_PI/3.0) - B3;
        r3 = x * std::cos(t-2*M_PI/3.0) - B3;

        return 3;
    }
    else if ( RR == QQQ )
    {
        // 3 real roots
        real sq = std::sqrt(Q);
        
        if ( R > 0 )
        {
            r1 = sq - B3;
            r2 = sq - B3;
            r3 = -2 * sq - B3;
        }
        else
        {
            r1 = 2 * sq - B3;
            r2 = -sq - B3;
            r3 = -sq - B3;
        }
        return 3;
    }
    else
    {
        // 1 real root
        real sgn = -sign_real(R);
        real x = sgn * std::pow( abs_real(R)+std::sqrt(RR-QQQ), 1.0/3.0 );
        real y = Q / x;

        r1 = x + y - B3;
        return 1;
    }
}


    
int QuarticSolver::solveCubicUnsorted(real A, real B, real C, real D,
                                      cplx& r1, cplx& r2, cplx& r3)
{
    if ( A == 0.0 )
        return solveQuadratic(B, C, D, r1, r2);
    
    if ( D == 0.0 )
    {
        r1 = cplx( 0, 0 );
        return 1 + solveQuadratic(A, B, C, r2, r3);
    }
    
    if ( A != 1.0 )
    {
        B /= A;
        C /= A;
        D /= A;
    }
    
    // calculate determinant of cubic
    real B3 = B / 3.0;
    real Q = B3*B3 - C / 3.0, QQQ = Q*Q*Q;
    real R = B3 * ( B3*B3 - 0.5 * C ) + 0.5 * D, RR = R*R;
    
    if ( R == 0 && Q == 0 )
    {
        // 3 identical real roots
        r1 = cplx( -B3, 0 );
        r2 = cplx( -B3, 0 );
        r3 = cplx( -B3, 0 );
        
        return 3;
    }
    else if ( RR < QQQ )
    {
        /* This sqrt and division is safe, since RR >= 0, so QQQ > RR,    */
        /* so QQQ > 0.  The acos is also safe, since RR/QQQ < 1, and      */
        /* thus R/std::sqrt(QQQ) < 1.                                          */
        real theta = std::acos(R/std::sqrt(QQQ));
        /* This sqrt is safe, since QQQ >= 0, and thus Q >= 0             */
        
        real x = -2*std::sqrt(Q);
        real t = theta/3.0;
        
        r1 = cplx( x * std::cos(t)            - B3, 0 );
        r2 = cplx( x * std::cos(t+2*M_PI/3.0) - B3, 0 );
        r3 = cplx( x * std::cos(t-2*M_PI/3.0) - B3, 0 );
        
        return 3;
    }
    else if ( RR == QQQ )
    {
        // 3 real roots
        real sgn = sign_real(R);
        real x = sgn * std::sqrt(Q);
        
        r1 = cplx( x - B3, 0 );
        r2 = cplx( x - B3, 0 );
        r3 = cplx( -2 * x - B3, 0 );
        
        return 3;
    }
    else
    {
        // 1 real root
        real sgn = -sign_real(R);
        real x = sgn * std::pow( abs_real(R)+std::sqrt(RR-QQQ), 1.0/3.0 );
        real y = Q / x;
        
        r1 = cplx( x + y - B3, 0 );
        r2 = cplx( -0.5 * ( x + y ) - B3, -0.5 * std::sqrt(3) * abs_real( x - y ) );
        r3 = cplx( -0.5 * ( x + y ) - B3,  0.5 * std::sqrt(3) * abs_real( x - y ) );
        
        return 3;
    }
}


//----------------------------------------------------------------------------
#pragma mark Quartic


/**
 return number of roots found (0, 1, 2, 3, 4).
 The roots are unsorted

 
 METHOD:
 
 x^4 + Bx^3 + Cx^2 + Dx + E = 0
 
 B = b/a
 C = c/a
 D = d/a
 E = e/a
 
 (depressed quartic, like in Ferrari's method)
 alpha = I = -3(B^2) / 8 + C
 beta  = J = (B^3)/8 - BC/2 + D
 gamma = K = -3(B^4)/256 + C(B^2)/16 - BD/4 + E
 
 Solve equation for z (one solution is enough):
 z^3 + 2Iz^2 + (I^2 - 4K)z - J^2 = 0
 
 p = std::sqrt(z)
 r = -p
 q = (I + z - J/p)/2
 s = (I + z + J/p)/2
 
 Solve (there are four u's, two for each equation):
 u^2 + pu + q = 0
 u^2 + ru + s = 0
 
 Solution (x[]):
 for each root u[i], x[i] = u[i] - B/4
 
 http://en.wikipedia.org/wiki/Quartic_function
*/
int QuarticSolver::solveQuarticUnsorted(real A, real B, real C, real D, real E,
                                        real& r1, real& r2, real& r3, real& r4)
{
    if ( A == 0.0 )
        return solveCubicUnsorted(B, C, D, E, r1, r2, r3);
    
    if ( E == 0.0 )
    {
#if ( QUARTIC_VERBOSE > 0 )
        fprintf(stderr, " quartic has null root\n");
#endif
        // zero is root
        r1 = 0;
        return 1 + solveCubicUnsorted(A, B, C, D, r2, r3, r4);
    }
    
    if ( A != 1.0 )
    {
        B /= A;
        C /= A;
        D /= A;
        E /= A;
    }

#if ( QUARTIC_VERBOSE > 1 )
    fprintf(stderr, "quartic x^4  %+f x^3  %+f x^2  %+f x  %+f\n", B, C, D, E);
#endif
    
    real BB = B*B;
    
    // define new coefficients that lead to the equation X^4 + I*X^2 + J*X + K = 0,
    
    real I = C - BB * 0.375;
    real J = D + ( BB * 0.25 - C ) * B * 0.5;
    real K = E - (( (3.0/16.00) * BB - C ) * B * 0.25 + D ) * B * 0.25;
 
    real z1 = 0, z2 = 0, z3 = 0;
    
    if ( K == 0 )
    {
        // zero is a trivial root, and we can factorize: X * [ X^3 + I*X + J ] = 0,
        r1 = -0.25*B;
        int n = 1 + solveCubicUnsorted(1.0, 0.0, I, J, z1, z2, z3);
        r2 = z1 + r1;
        r3 = z2 + r1;
        r4 = z3 + r1;
#if ( QUARTIC_VERBOSE > 0 )
        fprintf(stderr, " K == 0\n");
#endif
        return n;
    }
    
    if ( J == 0 )
    {
        // the quartic is quadratic in Z=X^2:  Z^2 + I*Z + K = 0
        int m = solveQuadratic(1.0, I, K, z1, z2);
        int n = 0;
        
        if ( m > 0 && z1 >= 0 )
        {
            r1 =  std::sqrt(z1) - 0.25*B;
            r2 = -std::sqrt(z1) - 0.25*B;
            n  = 2;
        }
        
        if ( m > 1 && z2 >= 0 )
        {
            if ( n == 2 )
            {
                r3 =  std::sqrt(z2) - 0.25*B;
                r4 = -std::sqrt(z2) - 0.25*B;
                n = 4;
            }
            else
            {
                r1 =  std::sqrt(z2) - 0.25*B;
                r2 = -std::sqrt(z2) - 0.25*B;
                n = 2;
            }
        }
#if ( QUARTIC_VERBOSE > 0 )
        fprintf(stderr, " J == 0\n");
#endif
        return n;
    }
    
    real C1 = I + I, C2 = I*I - 4*K, C3 = -J*J;
    
    int nz = solveCubic(1.0, C1, C2, C3, z1, z2, z3);
    if ( nz )
    {
        // z1 should be positive, because cubic(0) = -J*J < 0.
        // but this may not be true with computer precision:
        real ss = ( z1 > 0 ) ? z1 : 0;
        refineCubic(1.0, C1, C2, C3, ss);
        refineCubic(1.0, C1, C2, C3, ss);
        refineCubic(1.0, C1, C2, C3, ss);

#if ( QUARTIC_VERBOSE > 0 )
        fprintf(stderr, " real %+f\n", ss);
#endif

        real p = std::sqrt(ss);
        real w = J/p;
        // alternative formula, which has lower precision:
        //real w = std::sqrt( abs_real( ( I + z1 )*( I + z1 ) - 4*K ) ) * ( J < 0 ? -1 : 1 );
        
        real q = ( I + ss - w ) * 0.5;
        real s = ( I + ss + w ) * 0.5;
        real x, y;
        
        int n = solveQuadratic(1.0, p, q, x, y);
        if ( n == 2 )
        {
            r1 = x - 0.25*B;
            r2 = y - 0.25*B;
        }
        
        int m = solveQuadratic(1.0,-p, s, x, y);
        if ( m == 2 )
        {
            if ( n == 2 )
            {
                r3 = x - 0.25*B;
                r4 = y - 0.25*B;
                n = 4;
            }
            else
            {
                r1 = x - 0.25*B;
                r2 = y - 0.25*B;
                n = 2;
            }
        }
        return n;
    }
    return 0;
}

    
    
int QuarticSolver::solveQuarticUnsorted(real A, real B, real C, real D, real E,
                                        cplx& r1, cplx& r2, cplx& r3, cplx& r4)
{
    if ( A == 0.0 )
        return solveCubicUnsorted(B, C, D, E, r1, r2, r3);
    
    if ( E == 0.0 )
    {
#if ( QUARTIC_VERBOSE > 0 )
        fprintf(stderr, " quartic has null root\n");
#endif
        // zero is root
        r1 = cplx(0, 0);
        return 1 + solveCubicUnsorted(A, B, C, D, r2, r3, r4);
    }
    
    if ( A != 1.0 )
    {
        B /= A;
        C /= A;
        D /= A;
        E /= A;
    }
    
#if ( QUARTIC_VERBOSE > 1 )
    fprintf(stderr, "quartic x^4  %+f x^3  %+f x^2  %+f x  %+f\n", B, C, D, E);
#endif
    
    real BB = B*B;
    
    // define new coefficients that lead to the equation X^4 + I*X^2 + J*X + K = 0,
    
    real I = C - BB * 0.375;
    real J = D + ( BB * 0.25 - C ) * B * 0.5;
    real K = E - (( (3.0/16.00) * BB - C ) * B * 0.25 + D ) * B * 0.25;

    cplx z1, z2, z3;
    
    if ( K == 0 )
    {
        // zero is a trivial root, and we can factorize: X * [ X^3 + I*X + J ] = 0,
        r1 = cplx(-0.25*B, 0);
        int n = 1 + solveCubicUnsorted(1.0, 0.0, I, J, z1, z2, z3);
        r2 = z1 + r1;
        r3 = z2 + r1;
        r4 = z3 + r1;
#if ( QUARTIC_VERBOSE > 0 )
        fprintf(stderr, " K == 0\n");
#endif
        return n;
    }
    
    if ( J == 0 )
    {
        // the quartic is quadratic in Z=X^2:  Z^2 + I*Z + K = 0
        int n = solveQuadratic(1.0, I, K, z1, z2);
        
        if ( n > 0 )
        {
            z3 = z1.root();
            r1 = cplx( z3.r - 0.25*B, z3.i );
            r2 = cplx(-z3.r - 0.25*B,-z3.i );
        }
        if ( n > 1 )
        {
            z3 = z2.root();
            r3 = cplx( z3.r - 0.25*B, z3.i );
            r4 = cplx(-z3.r - 0.25*B,-z3.i );
        }
#if ( QUARTIC_VERBOSE > 0 )
        fprintf(stderr, " J == 0\n");
#endif
        return 2*n;
    }
    
    real C1 = I + I, C2 = I*I - 4*K, C3 = -J*J;

    int nz = solveCubic(1.0, C1, C2, C3, z1, z2, z3);
    if ( nz )
    {
        // chose root with smallest imaginary part:
        if ( abs_real(z1.i) > abs_real(z2.i) ) z1 = z2;
        if ( abs_real(z1.i) > abs_real(z3.i) ) z1 = z3;
        real ss = z1.r;

#if ( QUARTIC_VERBOSE > 0 )
        fprintf(stderr, " real %+f\n", ss);
#endif

        real p = std::sqrt(ss);
        real w = J/p;
        // alternative formula, which has lower precision:
        //real w = std::sqrt( abs_real( ( I + z1 )*( I + z1 ) - 4*K ) ) * ( J < 0 ? -1 : 1 );
        
        real q = ( I + ss - w ) * 0.5;
        real s = ( I + ss + w ) * 0.5;
        
        int n = solveQuadratic(1.0, p, q, z1, z2);
        cplx off(-0.25*B, 0);
        
        if ( n == 2 )
        {
            r1 = z1 + off;
            r2 = z2 + off;
        }
        
        int m = solveQuadratic(1.0,-p, s, z1, z2);
        if ( m == 2 )
        {
            if ( n == 2 )
            {
                r3 = z1 + off;
                r4 = z2 + off;
                n = 4;
            }
            else
            {
                r1 = z1 + off;
                r2 = z2 + off;
                n = 2;
            }
        }
        return n;
    }
    return 0;
}


