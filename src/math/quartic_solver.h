// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef QUARTIC_SOLVER_H
#define QUARTIC_SOLVER_H

#ifndef REAL_H
    #include "real.h"
#endif

#include <cmath>

/// Data structures and function to solve polynomial equations of order 4 or less
namespace QuarticSolver
{
    /// complex type built from two real
    class cplx
    {
    public:
        real r, i;
        cplx() { r = 0; i = 0; }
        cplx(real rs, real is) { r = rs; i = is; }
        
        bool operator <(cplx c) const { return r < c.r; }
        bool operator >(cplx c) const { return r > c.r; }
        cplx operator -() const { return cplx(-r, -i); }
        real norm() const { return sqrt(r*r+i*i); }
        cplx root() const
        {
            real n = norm();
            if  ( i >= 0 )
                return cplx(sqrt(0.5*( n + r )),  sqrt(0.5*( n - r )));
            else
                return cplx(sqrt(0.5*( n + r )), -sqrt(0.5*( n - r )));
        }
    };
        
    inline cplx operator +(cplx a, cplx b) { return cplx(a.r+b.r, a.i+b.i); }
    inline cplx operator -(cplx a, cplx b) { return cplx(a.r-b.r, a.i-b.i); }
    inline cplx operator *(cplx a, cplx b) { return cplx(a.r*b.r-a.i*b.i, a.r*b.i+a.i*b.r); }
        
    inline cplx operator +(cplx a, real b) { return cplx(a.r+b, a.i); }
    inline cplx operator -(cplx a, real b) { return cplx(a.r-b, a.i); }
    inline cplx operator *(cplx a, real b) { return cplx(a.r*b, a.i*b); }
    inline cplx operator *(real a, cplx b) { return cplx(a*b.r, a*b.i); }

        
    //----------------------------------------------------------------------------
#pragma mark Utilities
        
        
    template < typename TYPE >
    TYPE quadratic(const real a, const real b, const real c, const TYPE x)
    {
        return ( a * x + b ) * x + c;
    }
    
    template < typename TYPE >
    TYPE cubic(const real a, const real b, const real c, const real d, const TYPE x)
    {
        return (( a * x + b ) * x + c ) * x + d;
    }
    
    template < typename TYPE >
    TYPE quartic(const real a, const real b, const real c, const real d, const real e, const TYPE x)
    {
        return ((( a * x + b ) * x + c ) * x + d ) * x + e;
    }
    
    /**
     Sort in decreasing order
     */
    template < typename TYPE >
    void sortInPlace(int n, TYPE& x1, TYPE& x2, TYPE& x3)
    {
        TYPE x;
        if ( n > 2 && x1 < x3 )
        {
            x  = x3;
            x3 = x1;
            x1 = x;
        }
        if ( n > 1 && x1 < x2 )
        {
            x  = x2;
            x2 = x1;
            x1 = x;
        }
        if ( n > 2 && x2 < x3 )
        {
            x  = x3;
            x3 = x2;
            x2 = x;
        }
    }
    
    /**
     Sort in decreasing order
     */
    template < typename TYPE >
    void sortInPlace(int n, TYPE& x1, TYPE& x2, TYPE& x3, TYPE& x4)
    {
        TYPE x;
        if ( n > 3 && x2 < x4 )
        {
            x  = x4;
            x4 = x2;
            x2 = x;
        }
        if ( n > 2 && x1 < x3 )
        {
            x  = x3;
            x3 = x1;
            x1 = x;
        }
        if ( n > 3 && x3 < x4 )
        {
            x  = x4;
            x4 = x3;
            x3 = x;
        }
        if ( n > 1 && x1 < x2 )
        {
            x  = x2;
            x2 = x1;
            x1 = x;
        }
        if ( n > 2 && x2 < x3 )
        {
            x  = x3;
            x3 = x2;
            x2 = x;
        }
    }

    //----------------------------------------------------------------------------
#pragma mark Solvers

    /// Solve Quadratic, return number of real roots found (0, 1, 2).
    int solveQuadratic(const real A, const real B, const real C, real& r1, real& r2 );
    
    /// Solve Quadratic, return number of roots found.
    int solveQuadratic(const real A, const real B, const real C, cplx& r1, cplx& r2 );
    

    /// Solve Cubic, return number of real roots found (0, 1, 3).
    int solveCubicUnsorted(real A, real B, real C, real D, real& r1, real& r2, real& r3);
    
    /// Solve Cubic in complex space, return number of roots found (0, 1, 3).
    int solveCubicUnsorted(real A, real B, real C, real D, cplx& r1, cplx& r2, cplx& r3);
        
    /**
     return number of real roots found (0, 1, 3).
     The roots are returned in decreasing order:
     r1 >= r2 >= r3
     */
    template < typename TYPE >
    int solveCubic(const real A, const real B, const real C, const real D,
                   TYPE& r1, TYPE& r2, TYPE& r3)
    {
        int n = solveCubicUnsorted(A, B, C, D, r1, r2, r3);
        
        sortInPlace(n, r1, r2, r3);
        
        return n;
    }
    

    /// Solve Quartic, return number of real roots found (0, 1, 3).
    int solveQuarticUnsorted(real A, real B, real C, real D, real E, real& r1, real& r2, real& r3, real& r4);
    
    /// Solve Quartic in complex space, return number of roots found (0, 1, 3).
    int solveQuarticUnsorted(real A, real B, real C, real D, real E, cplx& r1, cplx& r2, cplx& r3, cplx& r4);
    
    /**
     return number of roots found (0, 1, 2, 3, 4).
     The roots are returned in decreasing order:
     r1 >= r2 >= r3 >= r4
     */
    template < typename TYPE >
    int solveQuartic(const real A, const real B, const real C, const real D, const real E,
                     TYPE& r1, TYPE& r2, TYPE& r3, TYPE& r4)
    {
        int n = solveQuarticUnsorted(A, B, C, D, E, r1, r2, r3, r4);
        
        sortInPlace(n, r1, r2, r3, r4);
        
        return n;
    }
        
            
    //----------------------------------------------------------------------------
    
    /// apply one step of Halley's method (convergence is cubic in general)
    void refineQuadratic(const real A, const real B, const real C, real& x);
    
    /// apply one step of Halley's method (convergence is cubic in general)
    void refineCubic(const real A, const real B, const real C, const real D, real& x);
    
    /// apply one step of Halley's method (convergence is cubic in general)
    void refineQuartic(const real A, const real B, const real C, const real D, const real E, real& x);

}

#endif
