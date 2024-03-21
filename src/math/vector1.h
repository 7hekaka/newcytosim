// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef VECTOR1_H
#define VECTOR1_H


#include "real.h"
#include "assert_macro.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cmath>

/// Vector1 is a vector with 1 `real` component.
class Vector1 final
{
    
public:
    
    /// dimensionality is 1
    static size_t dimensionality() { return 1; }
    
    /// coordinate is public
    real XX;
    
    
    /// by default, coordinates are not initialized
    Vector1() {}
    
    /// construct from 3 values
    constexpr Vector1(real x, real, real) : XX(x) {}
    
    /// construct from 1 value
    constexpr explicit Vector1(real x) : XX(x) {}

    /// construct from address
    Vector1(const real v[]) : XX(v[0]) {}    
    
    
    /// address of coordinate array
    real * data()                { return &XX; }
    
    /// constant address of coordinate array
    real const* data()     const { return &XX; }
    
    /// implicit conversion to a modifiable real pointer
    operator real*()             { return &XX; }
    
    /// implicit conversion to a constant real pointer
    operator const real*() const { return &XX; }
#if ( 0 )
    /// value of a coordinate
    real const& operator[](size_t i) const
    {
        assert_true(i==0);
        return XX;
    }
    
    /// modifiable access to individual coordinates
    real & operator[](size_t i)
    {
        assert_true(i==0);
        return XX;
    }
#endif
    
    /// return x-component
    real x() const { return XX; }
    /// return y-component
    real y() const { return 0; }
    /// return z-component
    real z() const { return 0; }

    /// copy at most one coordinate from array of size d
    void load(const real v[], int d)
    {
        XX = ( d > 0 ) ? v[0] : 0;
    }
    
    /// load from memory: X = b[0]
    void load(const float b[])
    {
        XX = b[0];
    }
    
    /// load from memory: X = b[0]
    void load(const double b[])
    {
        XX = b[0];
    }
    
    /// load difference: X = b[1] - b[0]
    void load_diff(const float b[])
    {
        XX = b[1] - b[0];
    }
    
    /// load difference: X = b[1] - b[0]
    void load_diff(const double b[])
    {
        XX = b[1] - b[0];
    }
    
    /// load difference: X = a[0] - b[0]
    void load_diff(const float a[], const float b[])
    {
        XX = a[0] - b[0];
    }

    /// load difference: X = a[0] - b[0]
    void load_diff(const double a[], const double b[])
    {
        XX = a[0] - b[0];
    }

    /// copy coordinates to given array
    void store(float b[]) const
    {
        b[0] = (float)XX;
    }
    
    /// copy coordinates to given array
    void store(double b[]) const
    {
        b[0] = (double)XX;
    }
    
    /// add content to given address
    void add_to(real b[]) const
    {
        b[0] += XX;
    }
    
    /// add content scaled by `alpha` to given address
    void add_to(real alpha, real b[]) const
    {
        b[0] += alpha * XX;
    }
    
    /// add content `n` times to array `b` of size `ldd*n`
    void add_to(real b[], int n, int ldd) const
    {
        for ( int i = 0; i < n; ++i )
            b[ldd*i] += XX;
    }
    
    /// subtract content to given address
    void sub_to(real b[]) const
    {
        b[0] -= XX;
    }
    
    /// subtract content scaled by `alpha` to given address
    void sub_to(real alpha, real b[]) const
    {
        b[0] -= alpha * XX;
    }
    
    /// set coordinates to zero
    void reset()
    {
        XX = 0;
    }
    
    /// change coordinates
    void set(const real x)
    {
        XX = x;
    }
    
    /// change coordinates (last 2 arguments are discarded)
    void set(const real x, const real, const real)
    {
        XX = x;
    }
    
    /// change signs of all coordinates
    void negate()
    {
        XX = -XX;
    }
    
    //------------------------------------------------------------------
    
    /// the square of the standard norm
    real normSqr() const
    {
        return XX*XX;
    }
    
    /// the square of the standard norm, minus TT*TT
    real normSqrSubtracted(const real& TT) const
    {
        return XX*XX - TT*TT;
    }

    /// the square of the norm
    friend real normSqr(Vector1 const& V)
    {
        return V.normSqr();
    }

    
    /// the standard norm = std::sqrt(x^2)
    real norm() const
    {
        return abs_real(XX);
    }

    /// the standard norm = std::sqrt(x^2)
    friend real norm(Vector1 const& V)
    {
        return V.norm();
    }
    
    /// the inversed magnitude = 1.0 / abs(x)
    real inv_norm() const
    {
        return 1 / abs_real(XX);
    }
    
    /// the 2D norm = std::sqrt(x^2+y^2)
    real normXYSqr() const
    {
        return XX*XX;
    }

    /// the 2D norm = std::sqrt(x^2+y^2)
    real normXY() const
    {
        return abs_real(XX);
    }
    
    /// the 2D norm = 0 since Y = Z = 0
    real normYZ() const
    {
        return 0;
    }
    
    /// the 2D norm = 0 since Y = Z = 0
    real normYZSqr() const
    {
        return 0;
    }

    /// square of the distance between two points, equivalent to (a-b).normSqr()
    friend real distanceSqr(Vector1 const& a, Vector1 const& b)
    {
        real x = a.XX - b.XX;
        return x*x;
    }
    
    /// distance between two points, equivalent to (a-b).norm()
    friend real distance(Vector1 const& a, Vector1 const& b)
    {
        return abs_real(a.XX-b.XX);
    }
 
    /// absolute values: (|x|)
    Vector1 abs() const
    {
        return Vector1(abs_real(XX));
    }

    /// the infinite norm = |x|
    real norm_inf() const
    {
        return abs_real(XX);
    }
    
    /// true if no component is NaN
    bool valid() const
    {
        return ( XX == XX );
    }
    
    /// true if component is not zero
    bool is_not_zero() const
    {
        return XX != 0.0;
    }
    
    /// scale to unit norm
    void normalize()
    {
        XX = sign_real(XX);
    }

    /// scale to obtain norm `n`
    void normalize(const real n)
    {
        XX = std::copysign(n, XX);
    }
    
    /// returns the colinear vector of norm `n` (default 1.0)
    Vector1 normalized(const real n = 1.0) const
    {
        return Vector1(std::copysign(n, XX));
    }
    
    /// returns vector parallel to argument of unit norm
    friend Vector1 normalize(Vector1 const& V)
    {
        return Vector1(sign_real(V.XX));
    }
    
    /// returns a perpendicular vector, of comparable but unspecified norm
    Vector1 orthogonal() const
    {
        ABORT_NOW("Vector::orthogonal() is meaningless in 1D");
        return Vector1(0.0);
    }
    
    /// returns a perpendicular vector, of norm `n`
    Vector1 orthogonal(const real) const
    {
        ABORT_NOW("Vector::orthogonal() is meaningless in 1D");
        return Vector1(0.0);
    }
    
    /// returns a vector perpendicular to *this, close to `d` and of norm = `n`
    Vector1 orthogonal(Vector1 const&, const real n) const
    {
        ABORT_NOW("Vector::orthogonal() is meaningless in 1D");
        return Vector1(n);
    }
    
    /// convert from cartesian to spherical coordinates ( r, theta, phi )
    Vector1 spherical() const { return Vector1(XX); }
    
    /// convert from spherical to cartesian coordinates ( x, y, z )
    Vector1 cartesian() const { return Vector1(XX); }
    
    //------------------------------------------------------------------
    
    /// Calculate intermediate position = A + C * ( B - A )
    void interpolate(const float a[], const float C, const float b[])
    {
        XX = a[0] + C * ( b[0] - a[0] );
    }
    
    /// Calculate intermediate position = A + C * ( B - A )
    void interpolate(const double a[], const double C, const double b[])
    {
        XX = a[0] + C * ( b[0] - a[0] );
    }
    
    /// linear interpolation, returning a + C * b
    friend Vector1 interpolated(const Vector1& a, real C, const Vector1& b)
    {
        return Vector1(a.XX+C*b.XX);
    }
    
    /// Calculate intermediate position = A + C * ( B - A )
    static Vector1 interpolated(const float a[], const float C, const float b[])
    {
        return Vector1(a[0]+C*(b[0]-a[0]));
    }

    /// Calculate intermediate position = A + C * ( B - A )
    static Vector1 interpolated(const double a[], const double C, const double b[])
    {
        return Vector1(a[0]+C*(b[0]-a[0]));
    }

    //------------------------------------------------------------------
    
    /// addition of two vectors
    friend Vector1 operator +(Vector1 const& a, Vector1 const& b)
    {
        return Vector1(a.XX+b.XX);
    }
    
    /// subtraction of two vectors
    friend Vector1 operator -(Vector1 const& a, Vector1 const& b)
    {
        return Vector1(a.XX-b.XX);
    }
    
    /// unary + operator does nothing
    friend Vector1 operator +(Vector1 const& b)
    {
        return b;
    }
    
    /// opposition of a vector
    friend Vector1 operator -(Vector1 const& b)
    {
        return Vector1(-b.XX);
    }
    
    /// returns the element-by-element product
    Vector1 e_mul(const Vector1& b) const
    {
        return Vector1(XX*b.XX);
    }
    
    /// returns the element-by-element division
    Vector1 e_div(const Vector1& b) const
    {
        return Vector1(XX/b.XX);
    }
    
    /// returns a vector with each element squared
    Vector1 e_squared() const
    {
        return Vector1(XX*XX);
    }
    
    /// returns sum of all coordinates
    real e_sum() const
    {
        return XX;
    }
    
    /// returns X
    real e_min() const
    {
        return XX;
    }
    
    /// returns X
    real e_max() const
    {
        return XX;
    }
    
    /// returns the element-by-element minimum
    Vector1 e_min(Vector1 const& v) const
    {
        return Vector1(std::min(XX, v.XX));
    }
    
    /// returns the element-by-element maximum
    Vector1 e_max(Vector1 const& v) const
    {
        return Vector1(std::max(XX, v.XX));
    }
    
    
    /**
     In dimension 1, the vector product is not really useful,
     but it is defined for completeness with the other class Vector2, Vector3.
     */
    
    /// the cross product of two vectors is a Z-Vector
    friend real cross(Vector1 const&, Vector1 const&)
    {
        return 0;
    }
    
    /// cross product of a vector with a Z-Vector
    friend Vector1 cross(Vector1 const&, const real)
    {
        return Vector1(0.0);
    }
    
    /// cross product of a Z-vector with a Vector
    friend Vector1 cross(const real, Vector1 const&)
    {
        return Vector1(0.0);
    }
    
    /// scalar product of two vectors
    friend real dot(Vector1 const& a, Vector1 const& b)
    {
        return a.XX * b.XX;
    }
    
    /// multiplication by scalar
    friend Vector1 operator *(Vector1 const& a, const real s)
    {
        return Vector1(s*a.XX);
    }
    
    /// mutiplication by scalar
    friend Vector1 operator *(const real s, Vector1 const& a)
    {
        return Vector1(s*a.XX);
    }
    
    /// division by scalar
    friend Vector1 operator /(Vector1 const& a, const real s)
    {
        return Vector1(a.XX/s);
    }
    
    /// addition of another vector
    void operator +=(Vector1 const& b)
    {
        XX += b.XX;
    }
    
    /// subtraction of another vector
    void operator -=(Vector1 const& b)
    {
        XX -= b.XX;
    }
    
    /// multiplication by a scalar
    void operator *=(const real s)
    {
        XX *= s;
    }
    
    /// division by a scalar
    void operator /=(const real s)
    {
        XX /= s;
    }
    
    //------------------------------------------------------------------
    
    /// equality test
    friend bool operator ==(Vector1 const& a, Vector1 const& b)
    {
        return ( a.XX==b.XX );
    }
    
    /// non-equality test
    friend bool operator !=(Vector1 const& a, Vector1 const& b)
    {
        return ( a.XX!=b.XX );
    }
    
    //------------------------------------------------------------------
    
    /// output
    void print(std::ostream& os) const
    {
        os << XX;
    }
    
    /// output using width 'w' and precision 'p'
    void print(std::ostream& os, int w, int p) const
    {
        os.precision(p);
        os << std::setw(w) << XX;
    }

    /// conversion to a string
    std::string to_string() const
    {
        std::ostringstream oss;
        print(oss);
        return oss.str();
    }
    
    /// conversion to a string with given precision
    std::string to_string(int w, int p) const
    {
        std::ostringstream oss;
        print(oss, w, p);
        return oss.str();
    }

    /// print to a file
    void print(FILE * out = stdout) const
    {
        fprintf(out, "  %+9.3f", XX);
    }
    
    /// print to a file, surrounded by parenthesis
    void pprint(FILE * out = stdout) const
    {
        fprintf(out, "( %+9.3f )", XX);
    }
    
    /// print, followed by a new line
    void println(FILE * out = stdout) const
    {
        fprintf(out, "  %+9.3f\n", XX);
    }
    
    //------------------------------------------------------------------
    
    /// add a random component in [-s, s] to each coordinate
    void addRand(real s);
    
    
    /// return null vector
    Vector1 randOrthoU(real) const;
    
    /// return null vector
    Vector1 randOrthoB(real) const;

    /// Vector with random independent coordinates in [0,+1]
    static Vector1 randP();
    
    /// Vector with random independent coordinates in [0,+n]
    static Vector1 randP(real n);
    
    /// Vector with random independent coordinates in [-1,+1]
    static Vector1 randS();
    
    /// Vector with random independent coordinates in [-1/2,+1/2]
    static Vector1 randH();
    
    /// Vector with random independent coordinates in [-n,+n]
    static Vector1 randS(real n);
    
    
    /// random Vector of norm = 1; sampling is uniform
    static Vector1 randU();
    
    /// return a random vector of norm = n; sampling is uniform
    static Vector1 randU(real n);
    
    
    /// return a random vector of norm <= 1; sampling is uniform
    static Vector1 randB();
    
    /// return a random vector of norm <= n; sampling is uniform
    static Vector1 randB(real n);
    
    
    /// return a random vector with Normally distributed coordinates ~ N(0,n)
    static Vector1 randG(real n);
    
};


//-------------------------- associated global functions -----------------------

/// stream input operator
std::istream& operator >> (std::istream&, Vector1&);

/// output operator
inline std::ostream& operator << (std::ostream& os, Vector1 const& arg)
{
    std::ios::fmtflags fgs = os.flags();
    arg.print(os);
    os.setf(fgs);
    return os;
}

#endif

