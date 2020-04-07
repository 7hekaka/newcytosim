// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <cstdio>
#include <sys/types.h>

#include "quartic_solver.h"
#include "random.h"

using namespace QuarticSolver;

/// expands ( u - s )
void expand(real& A, real& B, real s)
{
    A = 1.0;
    B = -s;
}

/// expands ( u - s )( u - s1 )
void expand(real& A, real& B, real& C, real s, real s1)
{
    expand(A, B, s1);
    C  = -s*B;
    B += -s*A;
}

/// expands ( u - s )( u - s1 )( u - s2 )
void expand(real& A, real& B, real& C, real& D, real s, real s1, real s2)
{
    expand(A, B, C, s1, s2);
    D  = -s*C;
    C += -s*B;
    B += -s*A;
}

/// expands ( u - s )( u - s1 )( u - s2 )( u - s3 )
void expand(real& A, real& B, real& C, real& D, real& E, real s, real s1, real s2, real s3)
{
    expand(A, B, C, D, s1, s2, s3);
    E  = -s*D;
    D += -s*C;
    C += -s*B;
    B += -s*A;
}

/// expands ( u*u + s )( u - s1 )( u - s2 )
void expand(real& A, real& B, real& C, real& D, real& E, real s, real s1, real s2)
{
    expand(A, B, C, s1, s2);
    E  = s*C;
    D  = s*B;
    C += s*A;
}


const real epsilon = 0.001;

void checkCubic(real a, real b, real c, real d, cplx x)
{
    cplx r = cubic(a, b, c, d, x);
    if ( abs_real(r.r) > epsilon || abs_real(r.i) > epsilon )
        fprintf(stderr, "      %+f  %+f -> %+f %+f\n", x.r, x.i, r.r, r.i);
}

void checkQuartic(real a, real b, real c, real d, real e, cplx x)
{
    cplx r = quartic(a, b, c, d, e, x);
    if ( abs_real(r.r) > epsilon || abs_real(r.i) > epsilon )
        fprintf(stderr, "      %+f  %+f -> %+f %+f\n", x.r, x.i, r.r, r.i);
}


void testQuartic(size_t cnt, const int DEG)
{
    real res = 0;
    unsigned miss = 0;
    
    const real p = 0.3;
    real A, B, C, D, E = 0;
    real s1, s2, s3, s4;
    real x1, x2, x3, x4;
    cplx z1, z2, z3, z4;
    
    for ( size_t u = 0; u < cnt; ++u )
    {
        x1 = RNG.sreal();
        x2 = RNG.test(p) ? x1 : RNG.sreal();
        x3 = RNG.test(p) ? x2 : RNG.sreal();
        x4 = RNG.test(p) ? x3 : RNG.sreal();
        
        int n = 0, m = 0;
        s1 = 0; s2 = 0; s3 = 0; s4 = 0;
        
#if ( 0 )
        A = 1.0;
        B = RNG.sreal();
        C = RNG.sreal();
        D = RNG.sreal();
        E = RNG.sreal();
#else
        if ( DEG == 3 )
        {
            expand(A, B, C, D, x1, x2, x3);
            //printf("cubic   %+f x^3  %+f x^2  %+f x  %+f\n", A, B, C, D);
        }
        else
        {
            expand(A, B, C, D, E, x1, x2, x3, x4);
            //printf("quartic %+f x^4  %+f x^3  %+f x^2  %+f x  %+f\n", A, B, C, D, E);
        }
#endif
        
        if ( DEG == 3 )
            n = solveCubic(A, B, C, D, s1, s2, s3);
        else
            n = solveQuartic(A, B, C, D, E, s1, s2, s3, s4);
        
        
        if ( n != DEG )
        {
            ++miss;
            printf("missed:\n");
            if ( DEG == 3 )
            {
                printf("x   %i :  %+f   %+f   %+f\n", DEG, x1, x2, x3);
                printf("s   %i :  %+f   %+f   %+f\n", n,   s1, s2, s3);
            } else {
                printf("x   %i :  %+f   %+f   %+f   %+f\n", DEG, x1, x2, x3, x4);
                printf("s   %i :  %+f   %+f   %+f   %+f\n", n,   s1, s2, s3, s4);
            }
        }

        if ( 0 )
        {
            //check the order of the roots:
            if ( n > 1 && s1 < s2 ) fprintf(stderr, " disorder s1 s2\n");
            if ( n > 2 && s2 < s3 ) fprintf(stderr, " disorder s2 s3\n");
            if ( n > 3 && s3 < s4 ) fprintf(stderr, " disorder s3 s4\n");
        }
        
        if ( DEG == 3 )
        {
            m = solveCubic(A, B, C, D, z1, z2, z3);
            //fprintf(stderr, "  equation has %i solutions:\n", n);
            if ( m > 0 ) checkCubic(A, B, C, D, z1);
            if ( m > 1 ) checkCubic(A, B, C, D, z2);
            if ( m > 2 ) checkCubic(A, B, C, D, z3);
        }
        if ( DEG == 4 )
        {
            m = solveQuartic(A, B, C, D, E, z1, z2, z3, z4);
            //fprintf(stderr, "  equation has %i solutions:\n", n);
            if ( m > 0 ) checkQuartic(A, B, C, D, E, z1);
            if ( m > 1 ) checkQuartic(A, B, C, D, E, z2);
            if ( m > 2 ) checkQuartic(A, B, C, D, E, z3);
            if ( m > 3 ) checkQuartic(A, B, C, D, E, z4);
        }
  
        real e1 = abs_real( x1 - z1.r );
        real e2 = abs_real( x2 - z2.r );
        real e3 = abs_real( x3 - z3.r );
        real e4 = abs_real( x4 - z4.r );
        
        bool stop = ( DEG > 3 && e4 > 0.01 );
        if ( e1 > 0.01 ) stop = true;
        if ( e2 > 0.01 ) stop = true;
        if ( e3 > 0.01 ) stop = true;
        if ( 0 ) //stop
        {
            printf("error:\n");
            if ( DEG == 3 )
            {
                printf("x   %i :  %+f   %+f   %+f\n", DEG, x1, x2, x3);
                printf("re  %i :  %+f   %+f   %+f\n", m,   z1.r, z2.r, z3.r);
                printf("im  %i :  %+f   %+f   %+f\n", m,   z1.i, z2.i, z3.i);
            } else {
                printf("x   %i :  %+f   %+f   %+f   %+f\n", DEG, x1, x2, x3, x4);
                printf("re  %i :  %+f   %+f   %+f   %+f\n", m,   z1.r, z2.r, z3.r, z4.r);
                printf("im  %i :  %+f   %+f   %+f   %+f\n", m,   z1.i, z2.i, z3.i, z4.i);
            }
        }
        
        
        if ( n > 0 )
        {
            real r;
            if ( DEG == 3 )
                r = abs_real(cubic(A, B, C, D, s1));
            else
                r = abs_real(quartic(A, B, C, D, E, s1));
            
            if ( r > res )
                res = r;
            
            if ( r > 0.5  ||  s1 != s1 )
            {
                printf("ERROR: q(%e) = %e\n", s1, r);
                printf("A = %.36f\n", A);
                printf("B = %.36f\n", B);
                printf("C = %.36f\n", C);
                printf("D = %.36f\n", D);
                printf("E = %.36f\n", E);
                return;
            }
        }
    }
    printf("max quartic residual = %e  misses %i\n", res, miss);
}


int main(int argc, char* argv[])
{
    RNG.seed();
    testQuartic(1<<14, 4);
    printf("test complete\n");
    return 0;
}
