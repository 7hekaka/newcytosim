// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef VECTOR4_H
#define VECTOR4_H


#include "real.h"
#include "assert_macro.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cmath>

class Vector1;
class Vector2;
class Vector3;

#ifdef __AVX__
#  define VECTOR4_USES_AVX REAL_IS_DOUBLE
#  include "simd.h"
#else
#  define VECTOR4_USES_AVX 0
#endif

/// Vector4 is a vector with 4 `real` components.
/**
 Note: We assume that the coordinates XX, YY, ZZ and TT are adjacent in memory,
 allowing easy conversion operators to and from C-array.
 Although this is not guaranteed by the C-standard, this is usually the case.
 */
class Vector4 final
{
    
public:
    
    /// dimensionality is 4
    static size_t dimensionality() { return 4; }
    
    /// coordinates are public
#if VECTOR4_USES_AVX
    union {
        struct {
            real XX;
            real YY;
            real ZZ;
            real TT;
        };
        vec4 vec;
    };
#else
    real XX;
    real YY;
    real ZZ;
    real TT;
#endif

    /// by default, coordinates are not initialized
    Vector4() {}
    
    /// construct from 3 values
    explicit constexpr Vector4(real x, real y, real z) : XX(x), YY(y), ZZ(z), TT(0.0) {}

    /// construct from 4 values
    explicit constexpr Vector4(real x, real y, real z, real t) : XX(x), YY(y), ZZ(z), TT(t) {}
    
    /// construct from address
    Vector4(const real v[]) : XX(v[0]), YY(v[1]), ZZ(v[2]), TT(v[3]) {}
    
#ifdef __AVX__
    /// construct from SIMD vector
    Vector4(vec4 const& v) { XX = v[0]; YY = v[1]; ZZ = v[2]; TT = v[3]; }
#endif
    
    /// copy 1 coordinate from Vector1
    explicit Vector4(const Vector1&);
    
    /// copy 2 coordinates from Vector2
    explicit Vector4(const Vector2&);
    
    /// copy 3 coordinates from Vector3
    explicit Vector4(const Vector3&);
    
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
        assert_true(i<4);
        return (&XX)[i];
    }
    
    /// modifiable access to individual coordinates
    real& operator[](size_t i)
    {
        assert_true(i<4);
        return (&XX)[i];
    }
#endif

    /// return x-component
    real x() const { return XX; }
    /// return y-component
    real y() const { return YY; }
    /// return z-component
    real z() const { return ZZ; }
    /// return t-component
    real t() const { return TT; }

    /// copy coordinates from array of size d
    void load(const real v[], const int& d)
    {
        XX = ( d > 0 ) ? v[0] : 0;
        YY = ( d > 1 ) ? v[1] : 0;
        ZZ = ( d > 2 ) ? v[2] : 0;
        TT = ( d > 3 ) ? v[3] : 0;
    }
    
    /// load from memory: X = b[0]; Y = b[1]; Z = b[2]; T = b[3]
    void load(const float b[])
    {
        XX = b[0];
        YY = b[1];
        ZZ = b[2];
        TT = b[3];
    }
    
    /// load from memory: X = b[0]; Y = b[1]; Z = b[2]; T = b[3]
    void load(const double b[])
    {
#if VECTOR4_USES_AVX
        vec = loadu4(b);
#else
        XX = b[0];
        YY = b[1];
        ZZ = b[2];
        TT = b[3];
#endif
    }
    
    /// load difference: X = b[4] - b[0]; Y = b[5] - b[1]; Z = b[6] - b[2]; T = b[7] - b[3]
    void load_diff(const float b[])
    {
        XX = b[4] - b[0];
        YY = b[5] - b[1];
        ZZ = b[6] - b[2];
        TT = b[7] - b[3];
    }
    
    /// load difference: X = b[4] - b[0]; Y = b[5] - b[1]; Z = b[6] - b[2]; T = b[7] - b[3]
    void load_diff(const double b[])
    {
#if VECTOR4_USES_AVX
        vec = sub4(loadu4(b+3), loadu4(b));
#else
        XX = b[4] - b[0];
        YY = b[5] - b[1];
        ZZ = b[6] - b[2];
        TT = b[7] - b[3];
#endif
    }
        
    /// load difference: X = a[0] - b[0]; Y = a[1] - b[1]; Z = a[2] - b[2]; T = a[3] - b[3]
    void load_diff(const float a[], const float b[])
    {
        XX = a[0] - b[0];
        YY = a[1] - b[1];
        ZZ = a[2] - b[2];
        TT = a[3] - b[3];
    }
    
    /// load difference: X = a[0] - b[0]; Y = a[1] - b[1]; Z = a[2] - b[2]; T = a[3] - b[3]
    void load_diff(const double a[], const double b[])
    {
#if VECTOR4_USES_AVX
        vec = sub4(loadu4(a), loadu4(b));
#else
        XX = a[0] - b[0];
        YY = a[1] - b[1];
        ZZ = a[2] - b[2];
        TT = a[3] - b[3];
#endif
    }
    
    /// copy coordinates to given array
    void store(float b[]) const
    {
        b[0] = (float)XX;
        b[1] = (float)YY;
        b[2] = (float)ZZ;
        b[3] = (float)TT;
    }
    
    /// copy coordinates to given array
    void store(double b[]) const
    {
#if VECTOR4_USES_AVX
        storeu4(b, vec);
#else
        b[0] = (double)XX;
        b[1] = (double)YY;
        b[2] = (double)ZZ;
        b[3] = (double)TT;
#endif
    }
    
    /// add content to given address
    void add_to(real b[]) const
    {
        b[0] += XX;
        b[1] += YY;
        b[2] += ZZ;
        b[3] += TT;
    }
    
    /// add content scaled by `alpha` to given address
    void add_to(real alpha, real b[]) const
    {
        b[0] += alpha * XX;
        b[1] += alpha * YY;
        b[2] += alpha * ZZ;
        b[3] += alpha * TT;
    }
    
    /// add content `n` times to array `b` of size `ldd*n`
    void add_to(real b[], int n, int ldd) const
    {
        for ( int i = 0; i < n; ++i )
        {
            b[ldd*i  ] += XX;
            b[ldd*i+1] += YY;
            b[ldd*i+2] += ZZ;
            b[ldd*i+3] += TT;
        }
    }
    
    /// subtract to given address
    void sub_to(real b[]) const
    {
        b[0] -= XX;
        b[1] -= YY;
        b[2] -= ZZ;
        b[3] -= TT;
    }
    
    /// subtract content scaled by `alpha` to given address
    void sub_to(real alpha, real b[]) const
    {
        b[0] -= alpha * XX;
        b[1] -= alpha * YY;
        b[2] -= alpha * ZZ;
        b[3] -= alpha * TT;
    }
    
    /// set coordinates to zero
    void reset()
    {
#if VECTOR4_USES_AVX
        vec = setzero4();
#else
        XX = 0;
        YY = 0;
        ZZ = 0;
        TT = 0;
#endif
    }
    
    /// change coordinates
    void set(const real x, const real y, const real z)
    {
        XX = x;
        YY = y;
        ZZ = z;
        TT = 0;
    }

    /// change coordinates
    void set(const real x, const real y, const real z, const real t)
    {
        XX = x;
        YY = y;
        ZZ = z;
        TT = t;
    }
    
    /// change signs of all coordinates
    void negate()
    {
#if VECTOR4_USES_AVX
        vec = flipsign4(vec);
#else
        XX = -XX;
        YY = -YY;
        ZZ = -ZZ;
        TT = -TT;
#endif
    }
    
    //------------------------------------------------------------------
    
    /// the square of the standard norm
    real normSqr() const
    {
        return XX*XX + YY*YY + ZZ*ZZ + TT*TT;
    }
    
    /// the square of the norm
    friend real normSqr(Vector4 const& V)
    {
        return V.normSqr();
    }

    /// the standard norm = std::sqrt(x^2+y^2+z^2+t^2)
    real norm() const
    {
        return std::sqrt(XX*XX + YY*YY + ZZ*ZZ + TT*TT);
    }
    
    /// the standard norm = std::sqrt(x^2+y^2+z^2+t^2)
    friend real norm(Vector4 const& V)
    {
        return V.norm();
    }
    
    /// the inversed magnitude = 1.0 / std::sqrt(x^2+y^2+z^2+t^2)
    real inv_norm() const
    {
        return 1 / std::sqrt(XX*XX + YY*YY + ZZ*ZZ + TT*TT);
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
    
    /// the 2D norm = std::sqrt(x^2+z^2)
    real normXZ() const
    {
        return std::sqrt(XX*XX + ZZ*ZZ);
    }
    
    /// the 2D norm = std::sqrt(y^2+z^2)
    real normYZ() const
    {
        return std::sqrt(YY*YY + ZZ*ZZ);
    }
    
    /// the 2D norm = std::sqrt(y^2+z^2)
    real normYZSqr() const
    {
        return YY*YY + ZZ*ZZ;
    }
    
    /// square of the distance between two points, equivalent to (a-b).normSqr()
    friend real distanceSqr(Vector4 const& a, Vector4 const& b)
    {
#if VECTOR4_USES_AVX
        return normsqr4(sub4(a.vec, b.vec))[0];
#else
        real x = a.XX - b.XX;
        real y = a.YY - b.YY;
        real z = a.ZZ - b.ZZ;
        real t = a.TT - b.TT;
        return x*x + y*y + z*z + t*t;
#endif
    }

    /// distance between two points, equivalent to (a-b).norm()
    friend real distance(Vector4 const& a, Vector4 const& b)
    {
        return std::sqrt(distanceSqr(a, b));
    }
    
    /// absolute values: (|x|, |y|, |z|, |t|)
    Vector4 abs() const
    {
        return Vector4(abs_real(XX), abs_real(YY), abs_real(ZZ), abs_real(TT));
    }

    /// the infinite norm = max(|x|, |y|, |z|)
    real norm_inf() const
    {
        return std::max(std::max(abs_real(XX), abs_real(YY)), std::max(abs_real(ZZ), abs_real(TT)));
    }
    
    /// true if no component is NaN
    bool valid() const
    {
        return ( XX == XX ) & ( YY == YY ) & ( ZZ == ZZ ) & ( TT == TT );
    }
    
    /// true if some component is not zero
    bool is_not_zero() const
    {
        return ( XX != 0.0 ) | ( YY != 0.0 ) | ( ZZ != 0.0 ) | ( TT != 0.0 );
    }
    
    /// scale to unit norm
    void normalize()
    {
#if VECTOR4_USES_AVX
        vec = normalize4(vec);
#else
        real s = norm();
        XX /= s;
        YY /= s;
        ZZ /= s;
        TT /= s;
#endif
    }
    
    /// scale to obtain norm `n`
    void normalize(const real n)
    {
#if VECTOR4_USES_AVX
        vec = normalize4(vec, n);
#else
        real s = n / norm();
        XX *= s;
        YY *= s;
        ZZ *= s;
        TT *= s;
#endif
    }

    /// returns the colinear vector of norm `n` (default 1.0)
    Vector4 normalized(const real n = 1.0) const
    {
#if VECTOR4_USES_AVX
        return Vector4(normalize4(vec, n));
#else
        real s = n / norm();
        return Vector4(s*XX, s*YY, s*ZZ, s*TT);
#endif
    }
    
    /// returns vector parallel to argument of unit norm
    friend Vector4 normalize(Vector4 const& V)
    {
#if VECTOR4_USES_AVX
        return Vector4(normalize4(V.vec));
#else
        const real s = V.norm();
        return Vector4(V.XX/s, V.YY/s, V.ZZ/s, V.TT/s);
#endif
    }

    //------------------------------------------------------------------
    
    /// Calculate intermediate position = A + C * ( B - A )
    void interpolate(const float a[], const float C, const float b[])
    {
        XX = a[0] + C * ( b[0] - a[0] );
        YY = a[1] + C * ( b[1] - a[1] );
        ZZ = a[2] + C * ( b[2] - a[2] );
        TT = a[3] + C * ( b[3] - a[3] );
    }

    /// Calculate intermediate position = A + C * ( B - A )
    void interpolate(const double a[], const double C, const double b[])
    {
#if VECTOR4_USES_AVX
        vec4 A = loadu4(a), B = loadu4(b);
        vec = fmadd4(set4(C), sub4(B, A), A);
#else
        XX = a[0] + C * ( b[0] - a[0] );
        YY = a[1] + C * ( b[1] - a[1] );
        ZZ = a[2] + C * ( b[2] - a[2] );
        TT = a[3] + C * ( b[3] - a[3] );
#endif
    }
    
    /// linear interpolation, returning a + x * b
    friend Vector4 interpolated(const Vector4& a, real C, const Vector4& b)
    {
        return Vector4(a.XX+C*b.XX, a.YY+C*b.YY, a.ZZ+C*b.ZZ, a.TT+C*b.TT);
    }
    
    /// Calculate intermediate position = A + C * ( B - A )
    static Vector4 interpolated(const float a[], const float C, const float b[])
    {
        return Vector4(a[0]+C*(b[0]-a[0]), a[1]+C*(b[1]-a[1]), a[2]+C*(b[2]-a[2]), a[3]+C*(b[3]-a[3]));
    }

    /// Calculate intermediate position = A + C * ( B - A )
    static Vector4 interpolated(const double a[], const double C, const double b[])
    {
#if VECTOR4_USES_AVX
        vec4 A = loadu4(a), B = loadu4(b);
        return Vector4(fmadd4(set4(C), sub4(B, A), A));
#else
        return Vector4(a[0]+C*(b[0]-a[0]), a[1]+C*(b[1]-a[1]), a[2]+C*(b[2]-a[2]), a[3]+C*(b[3]-a[3]));
#endif
    }

    //------------------------------------------------------------------
    
    /// addition of two vectors
    friend Vector4 operator +(Vector4 const& a, Vector4 const& b)
    {
        return Vector4(a.XX+b.XX, a.YY+b.YY, a.ZZ+b.ZZ, a.TT+b.TT);
    }
    
    /// subtraction of two vectors
    friend Vector4 operator -(Vector4 const& a, Vector4 const& b)
    {
        return Vector4(a.XX-b.XX, a.YY-b.YY, a.ZZ-b.ZZ, a.TT-b.TT);
    }
    
    /// unary + operator does nothing
    friend Vector4 operator +(Vector4 const& b)
    {
        return b;
    }
    
    /// opposition of a vector
    friend Vector4 operator -(Vector4 const& b)
    {
        return Vector4(-b.XX, -b.YY, -b.ZZ, -b.TT);
    }
    
    /// returns the element-by-element product
    Vector4 e_mul(Vector4 const& b) const
    {
        return Vector4(XX*b.XX, YY*b.YY, ZZ*b.ZZ, TT*b.TT);
    }

    /// returns the element-by-element division
    Vector4 e_div(Vector4 const& b) const
    {
        return Vector4(XX/b.XX, YY/b.YY, ZZ/b.ZZ, TT/b.TT);
    }
    
    /// returns a vector with each element squared
    Vector4 e_squared() const
    {
#if VECTOR4_USES_AVX
        return Vector4(mul4(vec, vec));
#else
        return Vector4(XX*XX, YY*YY, ZZ*ZZ, TT*TT);
#endif
    }
    
    /// returns sum of all coordinates
    real e_sum() const
    {
        return XX + YY + ZZ + TT;
    }
    
    /// returns min(x, y, z, t)
    real e_min() const
    {
        return std::min(std::min(XX, YY), std::min(ZZ, TT));
    }
    
    /// returns max(x, y, z, t)
    real e_max() const
    {
        return std::max(std::max(XX, YY), std::max(ZZ, TT));
    }
    
    /// returns the element-by-element minimum
    Vector4 e_min(Vector4 const& v) const
    {
        return Vector4(std::min(XX, v.XX), std::min(YY, v.YY), std::min(ZZ, v.ZZ), std::min(TT, v.TT));
    }
    
    /// returns the element-by-element maximum
    Vector4 e_max(Vector4 const& v) const
    {
        return Vector4(std::max(XX, v.XX), std::max(YY, v.YY), std::max(ZZ, v.ZZ), std::max(TT, v.TT));
    }
    
    /// scalar product of two vectors
    friend real dot(Vector4 const& a, Vector4 const& b)
    {
        return a.XX * b.XX + a.YY * b.YY + a.ZZ * b.ZZ + a.TT * b.TT;
    }
    
    /// 3D cross product of two vector (TT components are ignored)
    friend Vector4 cross(Vector4 const& a, Vector4 const& b)
    {
        return Vector4(a.YY * b.ZZ - a.ZZ * b.YY,
                       a.ZZ * b.XX - a.XX * b.ZZ,
                       a.XX * b.YY - a.YY * b.XX);
    }

    /// multiplication by scalar
    friend Vector4 operator *(Vector4 const& a, const real s)
    {
        return Vector4(s*a.XX, s*a.YY, s*a.ZZ, s*a.TT);
    }
    
    /// mutiplication by scalar
    friend Vector4 operator *(const real s, Vector4 const& a)
    {
        return Vector4(s*a.XX, s*a.YY, s*a.ZZ, s*a.TT);
    }
    
    /// division by scalar
    friend Vector4 operator /(Vector4 const& a, const real s)
    {
        return Vector4(a.XX/s, a.YY/s, a.ZZ/s, a.TT/s);
    }
    
    /// addition of another vector
    void operator +=(Vector4 const& b)
    {
        XX += b.XX;
        YY += b.YY;
        ZZ += b.ZZ;
        TT += b.TT;
    }
    
    /// subtraction of another vector
    void operator -=(Vector4 const& b)
    {
        XX -= b.XX;
        YY -= b.YY;
        ZZ -= b.ZZ;
        TT -= b.TT;
    }
    
    /// multiplication by a scalar
    void operator *=(const real s)
    {
        XX *= s;
        YY *= s;
        ZZ *= s;
        TT *= s;
    }
    
    /// division by a scalar
    void operator /=(const real s)
    {
        XX /= s;
        YY /= s;
        ZZ /= s;
        TT /= s;
    }
    
    //------------------------------------------------------------------
    
    /// equality test
    friend bool operator ==(Vector4 const& a, Vector4 const& b)
    {
        return ( a.XX==b.XX  &&  a.YY==b.YY  &&  a.ZZ==b.ZZ  &&  a.TT==b.TT );
    }
    
    /// non-equality test
    friend bool operator !=(Vector4 const& a, Vector4 const& b)
    {
        return ( a.XX!=b.XX  ||  a.YY!=b.YY  ||  a.ZZ!=b.ZZ  ||  a.TT!=b.TT );
    }
    
    //------------------------------------------------------------------
    
    /// output
    void print(std::ostream& os) const
    {
        const int w = (int)os.width();
        os << std::showpos << XX << " ";
        os << std::setw(w) << YY << " ";
        os << std::setw(w) << ZZ << " ";
        os << std::setw(w) << TT << std::noshowpos;
    }

    /// output using width 'w' and precision 'p'
    void print(std::ostream& os, int w, int p) const
    {
        os.precision(p);
        os.setf(std::ios::showpos);
        os << std::setw(w) << XX << " ";
        os << std::setw(w) << YY << " ";
        os << std::setw(w) << ZZ << " ";
        os << std::setw(w) << TT;
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
        fprintf(out, "  %+9.3f %+9.3f %+9.3f %+9.3f", XX, YY, ZZ, TT);
    }
    
    /// print to a file, surrounded by parenthesis
    void pprint(FILE * out = stdout) const
    {
        fprintf(out, "( %+9.3f %+9.3f %+9.3f %+9.3f )", XX, YY, ZZ, TT);
    }
    
    /// print, followed by a new line
    void println(FILE * out = stdout) const
    {
        fprintf(out, "  %+9.3f %+9.3f %+9.3f %+9.3f\n", XX, YY, ZZ, TT);
    }
    
    //------------------------------------------------------------------
    
    /// add a random component in [-s, s] to each coordinate
    void addRand(real s);
    
    
    /// Vector with random independent coordinates in [0,+1]
    static Vector4 randP();
    
    /// Vector with random independent coordinates in [0,+n]
    static Vector4 randP(real n);
    
    /// Vector with random independent coordinates in [-1,+1]
    static Vector4 randS();
    
    /// Vector with random independent coordinates in [-1/2,+1/2]
    static Vector4 randH();
    
    /// Vector with random independent coordinates in [-n,+n]
    static Vector4 randS(real n);
    
    
    /// random Vector of norm = 1; sampling is uniform
    static Vector4 randU();
    
    /// return a random vector of norm = n; sampling is uniform
    static Vector4 randU(real n);
    
    
    /// return a random vector of norm <= 1; sampling is uniform
    static Vector4 randB();
    
    /// return a random vector of norm <= n; sampling is uniform
    static Vector4 randB(real n);
    
    
    /// return a random vector with Normally distributed coordinates ~ N(0,n)
    static Vector4 randG(real n);
    
};

//-------------------------- associated global functions -----------------------

/// stream input operator
std::istream& operator >> (std::istream&, Vector4&);

/// output operator
inline std::ostream& operator << (std::ostream& os, Vector4 const& arg)
{
    std::ios::fmtflags fgs = os.flags();
    arg.print(os);
    os.setf(fgs);
    return os;
}

#endif
