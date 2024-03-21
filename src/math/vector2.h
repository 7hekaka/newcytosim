// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef VECTOR2_H
#define VECTOR2_H


#include "real.h"
#include "assert_macro.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cmath>

#ifdef __SSE3__
#  define VECTOR2_USES_SSE REAL_IS_DOUBLE
#  include "simd.h"
#else
#  define VECTOR2_USES_SSE 0
#endif


/// Vector2 is a vector with 2 `real` components.
/**
 Note: We assume that the coordinates XX and YY are adjacent in memory,
 allowing easy conversion operators to and from C-array.
 Although this is not guaranteed by the C-standard, this is usually the case.
 */
class Vector2 final
{
    
public:
    
    /// dimensionality is 2
    static size_t dimensionality() { return 2; }
    
    /// coordinates are public
#if VECTOR2_USES_SSE
    union {
        struct {
            real XX;
            real YY;
        };
        vec2 vec;
    };
#else
    real XX;
    real YY;
#endif

    /// by default, coordinates are not initialized
    Vector2() {}
    
    /// construct from 3 values
    explicit constexpr Vector2(real x, real y, real) : XX(x), YY(y) {}
    
    /// construct from 2 values
    explicit constexpr Vector2(real x, real y) : XX(x), YY(y) {}
    
    /// construct from address
    Vector2(const real v[]) : XX(v[0]), YY(v[1]) {}
    
#if VECTOR2_USES_SSE
    /// construct from SIMD vector
    Vector2(vec2 const& v) { vec = v; }
#elif defined(__SSE3__)
    /// construct from SIMD vector
    Vector2(vec2 const& v) { XX = v[0]; YY = v[1]; }
#endif
    
    
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
        assert_true(i<2);
        return (&XX)[i];
    }
    
    /// modifiable access to individual coordinates
    real& operator[](size_t i)
    {
        assert_true(i<2);
        return (&XX)[i];
    }
#endif

    /// return x-component
    real x() const { return XX; }
    /// return y-component
    real y() const { return YY; }
    /// return z-component
    real z() const { return 0; }

    /// copy coordinates from array of size d
    void load(const real v[], const int& d)
    {
        XX = ( d > 0 ) ? v[0] : 0;
        YY = ( d > 1 ) ? v[1] : 0;
    }
    
    /// load from memory: X = b[0]; Y = b[1]
    void load(const float b[])
    {
        XX = b[0];
        YY = b[1];
    }
    
    /// load from memory: X = b[0]; Y = b[1]
    void load(const double b[])
    {
#if VECTOR2_USES_SSE
        vec = loadu2(b);
#else
        XX = b[0];
        YY = b[1];
#endif
    }
    
    /// load difference: X = b[2] - b[0]; Y = b[3] - b[1]
    void load_diff(const float b[])
    {
        XX = b[2] - b[0];
        YY = b[3] - b[1];
    }
    
    /// load difference: X = b[2] - b[0]; Y = b[3] - b[1]
    void load_diff(const double b[])
    {
#if VECTOR2_USES_SSE
        vec = sub2(loadu2(b+2), loadu2(b));
#else
        XX = b[2] - b[0];
        YY = b[3] - b[1];
#endif
    }
    
    /// load difference: X = a[0] - b[0]; Y = a[1] - b[1]
    void load_diff(const float a[], const float b[])
    {
        XX = a[0] - b[0];
        YY = a[1] - b[1];
    }
    
    /// load difference: X = a[0] - b[0]; Y = a[1] - b[1]
    void load_diff(const double a[], const double b[])
    {
#if VECTOR2_USES_SSE
        vec = sub2(loadu2(a), loadu2(b));
#else
        XX = a[0] - b[0];
        YY = a[1] - b[1];
#endif
    }
    
    /// copy coordinates to given array
    void store(float b[]) const
    {
        b[0] = (float)XX;
        b[1] = (float)YY;
    }
    
    /// copy coordinates to given array
    void store(double b[]) const
    {
#if VECTOR2_USES_SSE
        store2(b, vec);
#else
        b[0] = (double)XX;
        b[1] = (double)YY;
#endif
    }
    
    /// add content to given address
    void add_to(real b[]) const
    {
#if VECTOR2_USES_SSE
        store2(b, add2(vec, load2(b)));
#else
        b[0] += XX;
        b[1] += YY;
#endif
    }
    
    /// add content scaled by `alpha` to given address
    void add_to(real alpha, real b[]) const
    {
#if VECTOR2_USES_SSE
        store2(b, fmadd2(set2(alpha), vec, load2(b)));
#else
        b[0] += alpha * XX;
        b[1] += alpha * YY;
#endif
    }
    
    /// add content `n` times to array `b` of size `ldd*n`
    void add_to(real b[], int n, int ldd) const
    {
        for ( int i = 0; i < n; ++i )
        {
            b[ldd*i  ] += XX;
            b[ldd*i+1] += YY;
        }
    }
    
    /// subtract to given address
    void sub_to(real b[]) const
    {
#if VECTOR2_USES_SSE
        store2(b, sub2(load2(b), vec));
#else
        b[0] -= XX;
        b[1] -= YY;
#endif
    }
    
    /// subtract content scaled by `alpha` to given address
    void sub_to(real alpha, real b[]) const
    {
#if VECTOR2_USES_SSE
        store2(b, fnmadd2(set2(alpha), vec, load2(b)));
#else
        b[0] -= alpha * XX;
        b[1] -= alpha * YY;
#endif
    }
    
    /// set coordinates to zero
    void reset()
    {
#if VECTOR2_USES_SSE
        vec = setzero2();
#else
        XX = 0;
        YY = 0;
#endif
    }
    
    /// change coordinates
    void set(const real x, const real y)
    {
#if VECTOR2_USES_SSE
        vec = setr2(x, y);
#else
        XX = x;
        YY = y;
#endif
    }
    
    /// change coordinates (last argument is discarded)
    void set(const real x, const real y, const real)
    {
#if VECTOR2_USES_SSE
        vec = setr2(x, y);
#else
        XX = x;
        YY = y;
#endif
    }
    
    /// change signs of all coordinates
    void negate()
    {
#if VECTOR2_USES_SSE
        vec = flipsign2(vec);
#else
        XX = -XX;
        YY = -YY;
#endif
    }
    
    //------------------------------------------------------------------
    
    /// the square of the standard norm
    real normSqr() const
    {
        return XX*XX + YY*YY;
    }
    
    /// the square of the standard norm, minus TT*TT
    real normSqrSubtracted(const real& TT) const
    {
        return ( XX*XX + YY*YY ) - TT*TT;
    }

    /// the square of the norm
    friend real normSqr(Vector2 const& V)
    {
        return V.normSqr();
    }

    
    /// the standard norm = std::sqrt(x^2+y^2)
    real norm() const
    {
        return std::sqrt(XX*XX + YY*YY);
    }
    
    /// the standard norm = std::sqrt(x^2+y^2)
    friend real norm(Vector2 const& V)
    {
        return V.norm();
    }
    
    /// the inversed magnitude = 1.0 / std::sqrt(x^2+y^2)
    real inv_norm() const
    {
        return 1 / std::sqrt(XX*XX + YY*YY);
    }
    
    /// the 2D norm = std::sqrt(x^2+y^2)
    real normXYSqr() const
    {
        return XX*XX + YY*YY;
    }

    /// the 2D norm = std::sqrt(x^2+y^2)
    real normXY() const
    {
        return std::sqrt(XX*XX + YY*YY);
    }
    
    /// the 2D norm = |y| since Z = 0
    real normYZ() const
    {
        return abs_real(YY);
    }
    
    /// the 2D norm = y^2 since Z = 0
    real normYZSqr() const
    {
        return YY*YY;
    }
    
    /// square of the distance between two points, equivalent to (a-b).normSqr()
    friend real distanceSqr(Vector2 const& a, Vector2 const& b)
    {
#if VECTOR2_USES_SSE
        return normsqr2(sub2(a.vec, b.vec))[0];
#else
        real x = a.XX - b.XX;
        real y = a.YY - b.YY;
        return x*x + y*y;
#endif
    }

    /// distance between two points, equivalent to (a-b).norm()
    friend real distance(Vector2 const& a, Vector2 const& b)
    {
        return std::sqrt(distanceSqr(a, b));
    }
    
    /// absolute values: (|x|, |y|)
    Vector2 abs() const
    {
        return Vector2(abs_real(XX), abs_real(YY));
    }

    /// the infinite norm = max(|x|, |y|)
    real norm_inf() const
    {
        return std::max(abs_real(XX), abs_real(YY));
    }
    
    /// true if no component is NaN
    bool valid() const
    {
        return ( XX == XX ) & ( YY == YY );
    }
    
    /// true if some component is not zero
    bool is_not_zero() const
    {
        return ( XX != 0.0 ) | ( YY != 0.0 );
    }
    
    /// scale to unit norm
    void normalize()
    {
#if VECTOR2_USES_SSE
        vec = normalize2(vec);
#else
        real s = norm();
        XX /= s;
        YY /= s;
#endif
    }

    /// scale to obtain norm `n`
    void normalize(const real n)
    {
#if VECTOR2_USES_SSE
        vec = normalize2(vec, n);
#else
        real s = n / norm();
        XX *= s;
        YY *= s;
#endif
    }
    
    /// returns the colinear vector of norm `n` (default 1.0)
    Vector2 normalized(const real n = 1.0) const
    {
#if VECTOR2_USES_SSE
        return Vector2(normalize2(vec, n));
#else
        real s = n / norm();
        return Vector2(s*XX, s*YY);
#endif
    }
    
    /// returns vector parallel to argument of unit norm
    friend Vector2 normalize(Vector2 const& V)
    {
#if VECTOR2_USES_SSE
        return Vector2(normalize2(V.vec));
#else
        const real s = V.norm();
        return Vector2(V.XX/s, V.YY/s);
#endif
    }

    //------------------------------------------------------------------
    
    /// returns a perpendicular vector, of same norm
    Vector2 orthogonal() const
    {
        return Vector2(-YY, XX);
    }
    
    /// returns a perpendicular vector, of norm `n`
    Vector2 orthogonal(const real n) const
    {
        real s = n / std::sqrt( XX*XX + YY*YY );
        return Vector2(-s*YY, s*XX);
    }
    
    /// returns a vector perpendicular to *this, close to `d` and of norm = `n`
    Vector2 orthogonal(Vector2 const& d, const real n) const
    {
        real s = dot(*this, d) / normSqr();
        return ( d - s * (*this) ).normalized(n);
    }
    
    /// convert from cartesian to spherical coordinates ( r, theta, phi )
    Vector2 spherical() const
    {
        return Vector2(std::sqrt(XX*XX+YY*YY), std::atan2(YY, XX));
    }
    
    /// convert from spherical to cartesian coordinates ( x, y, z )
    Vector2 cartesian() const
    {
        return Vector2(XX*std::cos(YY), XX*std::sin(YY));
    }
    
    //------------------------------------------------------------------
    
    /// Calculate intermediate position = A + C * ( B - A )
    void interpolate(const float a[], const float C, const float b[])
    {
        XX = a[0] + C * ( b[0] - a[0] );
        YY = a[1] + C * ( b[1] - a[1] );
    }
    
    /// Calculate intermediate position = A + C * ( B - A )
    void interpolate(const double a[], const double C, const double b[])
    {
#if VECTOR2_USES_SSE
        vec2 A = loadu2(a), B = loadu2(b);
#  ifdef __FMA__
        vec = fmadd2(set2(C), sub2(B, A), A);
#  else
        vec = add2(mul2(set2(C), sub2(B, A)), A);
#  endif
#else
        XX = a[0] + C * ( b[0] - a[0] );
        YY = a[1] + C * ( b[1] - a[1] );
#endif
    }
    
    /// linear interpolation, returning a + x * b
    friend Vector2 interpolated(const Vector2& a, real C, const Vector2& b)
    {
        return Vector2(a.XX+C*b.XX, a.YY+C*b.YY);
    }

    /// Calculate intermediate position = A + C * ( B - A )
    static Vector2 interpolated(const float a[], const float C, const float b[])
    {
        return Vector2(a[0]+C*(b[0]-a[0]), a[1]+C*(b[1]-a[1]));
    }

    /// Calculate intermediate position = A + C * ( B - A )
    static Vector2 interpolated(const double a[], const double C, const double b[])
    {
#if VECTOR2_USES_SSE
        vec2 A = loadu2(a), B = loadu2(b);
#  ifdef __FMA__
        return Vector2(fmadd2(set2(C), sub2(B, A), A));
#  else
        return Vector2(add2(mul2(set2(C), sub2(B, A)), A));
#  endif
#else
        return Vector2(a[0]+C*(b[0]-a[0]), a[1]+C*(b[1]-a[1]));
#endif
    }

    //------------------------------------------------------------------

    /// addition of two vectors
    friend Vector2 operator +(Vector2 const& a, Vector2 const& b)
    {
        return Vector2(a.XX+b.XX, a.YY+b.YY);
    }
    
    /// subtraction of two vectors
    friend Vector2 operator -(Vector2 const& a, Vector2 const& b)
    {
        return Vector2(a.XX-b.XX, a.YY-b.YY);
    }
    
    /// unary + operator does nothing
    friend Vector2 operator +(Vector2 const& b)
    {
        return b;
    }
    
    /// opposition of a vector
    friend Vector2 operator -(Vector2 const& b)
    {
        return Vector2(-b.XX, -b.YY);
    }
    
    /// returns the element-by-element product
    Vector2 e_mul(Vector2 const& b) const
    {
        return Vector2(XX*b.XX, YY*b.YY);
    }
    
    /// returns the element-by-element division
    Vector2 e_div(Vector2 const& b) const
    {
        return Vector2(XX/b.XX, YY/b.YY);
    }
    
    /// returns a vector with each element squared
    Vector2 e_squared() const
    {
#if VECTOR2_USES_SSE
        return Vector2(mul2(vec, vec));
#else
        return Vector2(XX*XX, YY*YY);
#endif
    }
    
    /// returns sum of all coordinates
    real e_sum() const
    {
        return XX + YY;
    }
    
    /// returns min(x, y)
    real e_min() const
    {
        return std::min(XX, YY);
    }
    
    /// returns max(x, y)
    real e_max() const
    {
        return std::max(XX, YY);
    }
    
    /// returns the element-by-element minimum
    Vector2 e_min(Vector2 const& v) const
    {
        return Vector2(std::min(XX, v.XX), std::min(YY, v.YY));
    }
    
    /// returns the element-by-element maximum
    Vector2 e_max(Vector2 const& v) const
    {
        return Vector2(std::max(XX, v.XX), std::max(YY, v.YY));
    }
    
    /**
     In dimension 2, we define a cross-product operator which returns a real,
     which in this case represents a Vector aligned with the Z axis.
     We also define the cross-product with a scalar, also corresponding to a
     Vector aligned with Z. This is a fair contraction of the 3D vector product.
     */
    
    /// the cross product of two vectors is a Z-Vector
    friend real cross(Vector2 const& a, Vector2 const& b)
    {
        return a.XX * b.YY - a.YY * b.XX;
    }
    
    /// cross product of a vector with a Z-Vector
    friend Vector2 cross(Vector2 const& a, const real b)
    {
        return Vector2(a.YY*b, -a.XX*b);
    }
    
    /// cross product of a Z-vector with a Vector
    friend Vector2 cross(const real a, Vector2 const& b)
    {
        return Vector2(-a*b.YY, a*b.XX);
    }
    
    /// scalar product of two vectors
    friend real dot(Vector2 const& a, Vector2 const& b)
    {
        return a.XX * b.XX + a.YY * b.YY;
    }
    
    /// multiplication by scalar
    friend Vector2 operator *(Vector2 const& a, const real s)
    {
        return Vector2(s*a.XX, s*a.YY);
    }
    
    /// mutiplication by scalar
    friend Vector2 operator *(const real s, Vector2 const& a)
    {
        return Vector2(s*a.XX, s*a.YY);
    }
    
    /// division by scalar
    friend Vector2 operator /(Vector2 const& a, const real s)
    {
        return Vector2(a.XX/s, a.YY/s);
    }
    
    /// addition of another vector
    void operator +=(Vector2 const& b)
    {
        XX += b.XX;
        YY += b.YY;
    }
    
    /// subtraction of another vector
    void operator -=(Vector2 const& b)
    {
        XX -= b.XX;
        YY -= b.YY;
    }
    
    /// multiplication by a scalar
    void operator *=(const real s)
    {
        XX *= s;
        YY *= s;
    }
    
    /// division by a scalar
    void operator /=(const real s)
    {
        XX /= s;
        YY /= s;
    }
    
    //------------------------------------------------------------------
    
    /// equality test
    friend bool operator ==(Vector2 const& a, Vector2 const& b)
    {
        return ( a.XX==b.XX  &&  a.YY==b.YY );
    }
    
    /// non-equality test
    friend bool operator !=(Vector2 const& a, Vector2 const& b)
    {
        return ( a.XX!=b.XX  ||  a.YY!=b.YY );
    }
    
    //------------------------------------------------------------------
    
    /// output
    void print(std::ostream& os) const
    {
        const int w = (int)os.width();
        os << std::showpos << XX << " ";
        os << std::setw(w) << YY << std::noshowpos;
    }
    
    /// output using width 'w' and precision 'p'
    void print(std::ostream& os, int w, int p) const
    {
        os.precision(p);
        os.setf(std::ios::showpos);
        os << std::setw(w) << XX << " ";
        os << std::setw(w) << YY;
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
        fprintf(out, "  %+9.3f %+9.3f", XX, YY);
    }
    
    /// print to a file, surrounded by parenthesis
    void pprint(FILE * out = stdout) const
    {
        fprintf(out, "( %+9.3f %+9.3f )", XX, YY);
    }
    
    /// print, followed by a new line
    void println(FILE * out = stdout) const
    {
        fprintf(out, "  %+9.3f %+9.3f\n", XX, YY);
    }
    
    //------------------------------------------------------------------
    
    /// add a random component in [-s, s] to each coordinate
    void addRand(real s);
    
    
    /// a vector of norm n, orthogonal to *this, assuming `norm(*this)==1`
    Vector2 randOrthoU(real n) const;
    
    /// a vector of norm <= n, orthogonal to *this, assuming `norm(*this)==1`
    Vector2 randOrthoB(real n) const;
    
    
    /// Vector with random independent coordinates in [0,+1]
    static Vector2 randP();
    
    /// Vector with random independent coordinates in [0,+n]
    static Vector2 randP(real n);
    
    /// Vector with random independent coordinates in [-1,+1]
    static Vector2 randS();
    
    /// Vector with random independent coordinates in [-1/2,+1/2]
    static Vector2 randH();
    
    /// Vector with random independent coordinates in [-n,+n]
    static Vector2 randS(real n);
    
    
    /// random Vector of norm = 1; sampling is uniform
    static Vector2 randU();
    
    /// return a random vector of norm = n; sampling is uniform
    static Vector2 randU(real n);
    
    
    /// return a random vector of norm <= 1; sampling is uniform
    static Vector2 randB();
    
    /// return a random vector of norm <= n; sampling is uniform
    static Vector2 randB(real n);
    
    
    /// return a random vector with Normally distributed coordinates ~ N(0,n)
    static Vector2 randG(real n);
    
};


//-------------------------- associated global functions -----------------------

/// stream input operator
std::istream& operator >> (std::istream&, Vector2&);

/// output operator
inline std::ostream& operator << (std::ostream& os, Vector2 const& arg)
{
    std::ios::fmtflags fgs = os.flags();
    arg.print(os);
    os.setf(fgs);
    return os;
}


#endif

