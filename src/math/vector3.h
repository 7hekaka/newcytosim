// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef VECTOR3_H
#define VECTOR3_H


#include "real.h"
#include "assert_macro.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cmath>

class Vector1;
class Vector2;

#ifdef __AVX__
#  define VECTOR3_USES_AVX REAL_IS_DOUBLE
#  include "simd.h"
#else
#  define VECTOR3_USES_AVX 0
#endif

/// Vector3 is a vector with 3 `real` components.
/**
 Note: We assume that the coordinates XX, YY and ZZ are adjacent in memory,
 allowing easy conversion operators to and from C-array.
 Although this is not guaranteed by the C-standard, this is usually the case.
 */
class Vector3 final
{
    
public:
    
    /// dimensionality is 3
    static size_t dimensionality() { return 3; }
    
    /// coordinates are public
#if VECTOR3_USES_AVX
    union {
        struct {
            real XX;
            real YY;
            real ZZ;
            real TT;  // should be zero!
        };
        vec4 vec;
    };
#else
    real XX;
    real YY;
    real ZZ;
#endif

    /// by default, coordinates are not initialized
    Vector3() { }

#if VECTOR3_USES_AVX
    /// construct from 3 values
    explicit constexpr Vector3(real x, real y, real z) : XX(x), YY(y), ZZ(z), TT(0.0) {}
    
    /// construct from address
    Vector3(const real v[]) : vec(load3Z(v)) {}
#else
    /// construct from 3 values
    explicit constexpr Vector3(real x, real y, real z) : XX(x), YY(y), ZZ(z) {}

    /// construct from address
    Vector3(const real v[]) : XX(v[0]), YY(v[1]), ZZ(v[2]) {}
#endif

#if VECTOR3_USES_AVX
    /// construct from SIMD vector
    Vector3(vec4 const& v) : vec(v) { } //{ assert_true(v[3]==0); }
    /// conversion to SIMD vector
    operator vec4 () const { assert_true(vec[3]==0); return vec; }
#elif defined(__AVX__) && REAL_IS_DOUBLE
    /// construct from SIMD vector
    Vector3(vec4 const& v) : XX(v[0]), YY(v[1]), ZZ(v[2]) { }
    /// conversion to SIMD vector
    operator vec4 () const { return load3(&XX); }
#endif
    
    /// copy 1 coordinate from Vector1
    explicit Vector3(const Vector1&);
    
    /// copy 2 coordinates from Vector2
    explicit Vector3(const Vector2&);
    
    
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
        assert_true(i<3);
        return (&XX)[i];
    }
    
    /// modifiable access to individual coordinates
    real& operator[](size_t i)
    {
        assert_true(i<3);
        return (&XX)[i];
    }
#endif

    /// return x-component
    real x() const { return XX; }
    /// return y-component
    real y() const { return YY; }
    /// return z-component
    real z() const { return ZZ; }

    /// copy coordinates from array of size d
    void load(const real v[], const int& d)
    {
        XX = ( d > 0 ) ? v[0] : 0;
        YY = ( d > 1 ) ? v[1] : 0;
        ZZ = ( d > 2 ) ? v[2] : 0;
#if VECTOR3_USES_AVX
        vec[3] = 0;
#endif
    }
    
    /// load from memory: X = b[0]; Y = b[1]; Z = b[2]
    void load(const float b[])
    {
        XX = b[0];
        YY = b[1];
        ZZ = b[2];
#if VECTOR3_USES_AVX
        vec[3] = 0;
#endif
    }
    
    /// load from memory: X = b[0]; Y = b[1]; Z = b[2]
    void load(const double b[])
    {
#if VECTOR3_USES_AVX
        vec = load3(b);
#else
        XX = b[0];
        YY = b[1];
        ZZ = b[2];
#endif
    }
    
    /// load difference: X = b[3] - b[0]; Y = b[4] - b[1]; Z = b[5] - b[2]
    void load_diff(const float b[])
    {
        XX = b[3] - b[0];
        YY = b[4] - b[1];
        ZZ = b[5] - b[2];
#if VECTOR3_USES_AVX
        vec[3] = 0;
#endif
    }
    
    /// load difference: X = b[3] - b[0]; Y = b[4] - b[1]; Z = b[5] - b[2]
    void load_diff(const double b[])
    {
#if VECTOR3_USES_AVX
        vec = sub4(load3(b+3), load3(b));
#else
        XX = b[3] - b[0];
        YY = b[4] - b[1];
        ZZ = b[5] - b[2];
#endif
    }
    
    /// load difference: X = a[0] - b[0]; Y = a[1] - b[1]; Z = a[2] - b[2]
    void load_diff(const float a[], const float b[])
    {
        XX = a[0] - b[0];
        YY = a[1] - b[1];
        ZZ = a[2] - b[2];
#if VECTOR3_USES_AVX
        vec[3] = 0;
#endif
    }
    
    /// load difference: X = a[0] - b[0]; Y = a[1] - b[1]; Z = a[2] - b[2]
    void load_diff(const double a[], const double b[])
    {
#if VECTOR3_USES_AVX
        vec = sub4(load3(a), load3(b));
#else
        XX = a[0] - b[0];
        YY = a[1] - b[1];
        ZZ = a[2] - b[2];
#endif
    }

    /// copy coordinates to given array
    void store(float b[]) const
    {
        b[0] = (float)XX;
        b[1] = (float)YY;
        b[2] = (float)ZZ;
    }
    
    /// copy coordinates to given array
    void store(double b[]) const
    {
#if VECTOR3_USES_AVX
        store3(b, vec);
#else
        b[0] = (double)XX;
        b[1] = (double)YY;
        b[2] = (double)ZZ;
#endif
    }
    
    /// add content to given address
    void add_to(real b[]) const
    {
#if VECTOR3_USES_AVX
        assert_true(vec[3] == 0);
        storeu4(b, add4(vec, loadu4(b)));
#else
        b[0] += XX;
        b[1] += YY;
        b[2] += ZZ;
#endif
    }
    
    /// add content scaled by `alpha` to given address
    void add_to(real alpha, real b[]) const
    {
#if VECTOR3_USES_AVX
        assert_true(vec[3] == 0);
        storeu4(b, add4(mul4(set4(alpha), vec), loadu4(b)));
#else
        b[0] += alpha * XX;
        b[1] += alpha * YY;
        b[2] += alpha * ZZ;
#endif
    }
    
    /// add content `n` times to array `b` of size `ldd*n`
    void add_to(real b[], int n, int ldd) const
    {
        for ( int i = 0; i < n; ++i )
        {
            b[ldd*i  ] += XX;
            b[ldd*i+1] += YY;
            b[ldd*i+2] += ZZ;
        }
    }
    
    /// subtract to given address
    void sub_to(real b[]) const
    {
#if VECTOR3_USES_AVX
        assert_true(vec[3] == 0);
        storeu4(b, sub4(loadu4(b), vec));
#else
        b[0] -= XX;
        b[1] -= YY;
        b[2] -= ZZ;
#endif
    }
    
    /// subtract content scaled by `alpha` to given address
    void sub_to(real alpha, real b[]) const
    {
#if VECTOR3_USES_AVX
        assert_true(vec[3] == 0);
        storeu4(b, sub4(loadu4(b), mul4(set4(alpha), vec)));
#else
        b[0] -= alpha * XX;
        b[1] -= alpha * YY;
        b[2] -= alpha * ZZ;
#endif
    }
    
    /// set coordinates to zero
    void reset()
    {
#if VECTOR3_USES_AVX
        vec = setzero4();
#else
        XX = 0;
        YY = 0;
        ZZ = 0;
#endif
    }
    
    /// change coordinates
    void set(const real x, const real y, const real z)
    {
#if VECTOR3_USES_AVX
        vec = setr4(x, y, z, 0);
#else
        XX = x;
        YY = y;
        ZZ = z;
#endif
    }
    
    /// change signs of all coordinates
    void negate()
    {
#if VECTOR3_USES_AVX
        vec = flipsign4(vec);
#else
        XX = -XX;
        YY = -YY;
        ZZ = -ZZ;
#endif
    }
    
    //------------------------------------------------------------------
    
    /// the square of the standard norm
    real normSqr() const
    {
        return XX*XX + YY*YY + ZZ*ZZ;
    }
    
    /// the square of the standard norm, minus TT*TT
    real normSqrSubtracted(const real& TT) const
    {
        return ( XX*XX + YY*YY ) + ( ZZ*ZZ - TT*TT );
    }

    /// the square of the norm
    friend real normSqr(Vector3 const& V)
    {
        return V.normSqr();
    }

    /// the standard norm = std::sqrt(x^2+y^2+z^2)
    real norm() const
    {
        return std::sqrt(XX*XX + YY*YY + ZZ*ZZ);
    }
    
    /// the standard norm = std::sqrt(x^2+y^2+z^2)
    friend real norm(Vector3 const& V)
    {
        return V.norm();
    }
    
    /// the inversed magnitude = 1.0 / std::sqrt(x^2+y^2+z^2)
    real inv_norm() const
    {
        return 1 / std::sqrt(XX*XX + YY*YY + ZZ*ZZ);
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
    
    /// the 2D norm = x^2+z^2
    real normXZSqr() const
    {
        return XX*XX + ZZ*ZZ;
    }

    /// the 2D norm = std::sqrt(x^2+z^2)
    real normXZ() const
    {
        return std::sqrt(XX*XX + ZZ*ZZ);
    }

    /// the 2D norm = y^2+z^2
    real normYZSqr() const
    {
        return YY*YY + ZZ*ZZ;
    }

    /// the 2D norm = std::sqrt(y^2+z^2)
    real normYZ() const
    {
        return std::sqrt(YY*YY + ZZ*ZZ);
    }

    /// square of the distance between two points, equivalent to (a-b).normSqr()
    friend real distanceSqr(Vector3 const& a, Vector3 const& b)
    {
#if VECTOR3_USES_AVX
        assert_true(a.vec[3] == 0);
        assert_true(b.vec[3] == 0);
        return normsqr4(sub4(a.vec, b.vec))[0];
#else
        real x = a.XX - b.XX;
        real y = a.YY - b.YY;
        real z = a.ZZ - b.ZZ;
        return x*x + y*y + z*z;
#endif
    }

    /// distance between two points, equivalent to (a-b).norm()
    friend real distance(Vector3 const& a, Vector3 const& b)
    {
        return std::sqrt(distanceSqr(a,b));
    }

    /// absolute values: (|x|, |y|, |z|)
    Vector3 abs() const
    {
        return Vector3(abs_real(XX), abs_real(YY), abs_real(ZZ));
    }

    /// the infinite norm = max(|x|, |y|, |z|)
    real norm_inf() const
    {
        return std::max(std::max(abs_real(XX), abs_real(YY)), abs_real(ZZ));
    }
    
    /// true if no component is NaN
    bool valid() const
    {
#if VECTOR3_USES_AVX
        if ( vec[3] != 0 ) return false;
#endif
        return ( XX == XX ) & ( YY == YY ) & ( ZZ == ZZ );
    }
    
    /// true if some component is not zero
    bool is_not_zero() const
    {
        return ( XX != 0.0 ) | ( YY != 0.0 ) | ( ZZ != 0.0 );
    }

    /// scale to unit norm
    void normalize()
    {
#if VECTOR3_USES_AVX
        assert_true(vec[3] == 0);
        vec = normalize4(vec);
#else
        real s = norm();
        XX /= s;
        YY /= s;
        ZZ /= s;
#endif
    }

    /// scale to obtain norm `n`
    void normalize(const real n)
    {
#if VECTOR3_USES_AVX
        assert_true(vec[3] == 0);
        vec = normalize4(vec, n);
#else
        real s = n / norm();
        XX *= s;
        YY *= s;
        ZZ *= s;
#endif
    }

    /// returns vector parallel to argument of unit norm
    friend Vector3 normalize(Vector3 const& V)
    {
#if VECTOR3_USES_AVX
        assert_true(V.vec[3] == 0);
        return Vector3(normalize4(V.vec));
#else
        const real s = V.norm();
        return Vector3(V.XX/s, V.YY/s, V.ZZ/s);
#endif
    }

    /// returns the colinear vector of norm `n` (default 1.0)
    Vector3 normalized(const real n = 1.0) const
    {
#if VECTOR3_USES_AVX
        assert_true(vec[3] == 0);
        return Vector3(normalize4(vec, n));
#else
        real s = n / norm();
        return Vector3(s*XX, s*YY, s*ZZ);
#endif
    }

    //------------------------------------------------------------------
    /// returns a perpendicular vector, of comparable but unspecified norm
    Vector3 orthogonalB() const
    {
        if ( abs_real(XX) < abs_real(YY) )
        {
            if ( abs_real(XX) < abs_real(ZZ) )
                return Vector3(0.0, -ZZ,  YY); //XX is the smallest
            else
                return Vector3( YY, -XX, 0.0); //ZZ is the smallest
        }
        else
        {
            if ( abs_real(YY) < abs_real(ZZ) )
                return Vector3(-ZZ, 0.0,  XX); //YY is the smallest
            else
                return Vector3( YY, -XX, 0.0); //ZZ is the smallest
        }
    }

    /// returns a perpendicular vector, of comparable but unspecified norm
    /**
     Presumably branchless code inspired from:
     Stark, M. M., “Efficient Construction of Perpendicular Vectors without Branching”,
     Journal of Graphics Tools 14:1 (2009), 55-61.
     */
    Vector3 orthogonal() const
    {
        real ax = abs_real(XX);
        real ay = abs_real(YY);
        real az = abs_real(ZZ);
        
        // select axis for which vector component is the smallest:
        real x = ax - std::min(ay, az); // negative if XX is smallest
        real y = ay - std::min(az, ax); // negative if YY is smallest
        real z = std::min(x, y); // use z if x and y are positive (could use bitwise OR)

        // form cross product with axis (only one of x or y can be negative)
        return Vector3(sign_select(z, 0, YY) - sign_select(y, ZZ, 0),
                       sign_select(x, ZZ, 0) - sign_select(z, 0, XX),
                       sign_select(y, XX, 0) - sign_select(x, YY, 0));
    }
    
#if ( 1 )
    /// returns a perpendicular Vector, of norm `n`
    Vector3 orthogonal(const real n) const
    {
        real ax = abs_real(XX);
        real ay = abs_real(YY);
        real az = abs_real(ZZ);
        
        real x = ax - std::min(ay, az); // negative if XX is smallest
        real y = ay - std::min(az, ax); // negative if YY is smallest
        real z = std::min(x, y); // use z if x and y are positive

        // form cross product with axis (only one of x or y can be negative)
        ax = sign_select(z, 0, YY) - sign_select(y, ZZ, 0);
        ay = sign_select(x, ZZ, 0) - sign_select(z, 0, XX);
        az = sign_select(y, XX, 0) - sign_select(x, YY, 0);
        // normalize:
        real s = n / std::sqrt( ax*ax + ay*ay + az*az );
        return Vector3(ax*s, ay*s, az*s);
    }
#else
    /// returns a perpendicular vector, of norm `n`
    Vector3 orthogonal(const real n) const
    {
        if ( abs_real(XX) < abs_real(YY) )
        {
            if ( abs_real(XX) < abs_real(ZZ) )
            {
                // XX is the smallest component
                real s = n / std::sqrt(YY*YY+ZZ*ZZ);
                return Vector3(0.0, -s*ZZ, s*YY);
            }
            else
            {
                // ZZ is the smallest component
                real s = n / std::sqrt(XX*XX+YY*YY);
                return Vector3(s*YY, -s*XX, 0.0);
            }
        }
        else
        {
            if ( abs_real(YY) < abs_real(ZZ) )
            {
                // YY is the smallest component
                real s = n / std::sqrt(XX*XX+ZZ*ZZ);
                return Vector3(-s*ZZ, 0.0, s*XX);
            }
            else
            {
                // ZZ is the smallest component
                real s = n / std::sqrt(XX*XX+YY*YY);
                return Vector3(s*YY, -s*XX, 0.0);
            }
        }
    }
#endif
    /// returns a vector perpendicular to *this, close to `d` and of norm = `n`
    /**
     This removes the component of `n` parallel to *this,
     and will fail if `d` is parallel to *this
     */
    Vector3 orthogonal(Vector3 const& d, const real n) const
    {
        real s = dot(*this, d) / normSqr();
        return ( d - s * (*this) ).normalized(n);
    }
    
    /**
     Given: N = norm(this), (C, S) = random numbers in [-1, 1],
     @return a vector orthogonal to *this, of norm `N * ( C^2 + S^2 )`
     
     Derived from `Building an Orthonormal Basis, Revisited`,
     Tom Duff et al. Journal of Computer Graphics Techniques Vol. 6 N.1, 2017
     */
    Vector3 orthogonal(real N, real C, real S) const
    {
        assert_small(norm()/N - 1.0);
        N = std::copysign(N, ZZ);
        real a = YY / ( ZZ + N );
        real b = YY * a;
        real c = XX * a;
        return Vector3(S*c-C*(ZZ+b), S*(b-N)+C*c, S*YY+C*XX);
    }

    /**
     Set vectors 'E' and 'F' to build an orthonormal basis (this, E, F),
     assuming that 'norm(*this) == 1'
     
     From `Building an Orthonormal Basis, Revisited`,
     Tom Duff et al. Journal of Computer Graphics Techniques Vol. 6 N.1, 2017
     */
    void orthonormal(real E[3], real F[3]) const
    {
#if 0
        if ( abs_real(1.0 - normSqr()) > 0.01 )
        {
            // this should not happen...
            E = orthogonal(1);
            F = cross(*this, E).normalized();
            std::clog << "rescued orthonormal(" << to_string() << ")\n";
            return;
        }
#else
        assert_small(1.0 - normSqr());
#endif
        real s = std::copysign(real(1.0), ZZ);
        // optimized version by Marc B. Reynolds
        const real a = YY / ( ZZ + s );
        const real b = YY * a;
        const real c = XX * a;
        // below normSqr(F) = normSqr(this) + a*a*(normSqr(this)-s*s)
        E[0] = -ZZ - b;
        E[1] = c;
        E[2] = XX;
        F[0] = s * c;
        F[1] = s * b - 1;
        F[2] = s * YY;
        // if you do not mind an inverted basis, you may use F.set(c, b-s, YY);
        //printf("orthonormal %+9.6f %+9.6f %+9.6f\n", dot(*this, E), dot(*this, F), dot(E, F));
    }
    
    /**
     Set 'E' and 'F' to build a basis (this, E, F), with norm(E) = norm(F) = N
     assuming that 'norm(*this) == 1'
     
     From `Building an Orthonormal Basis, Revisited`,
     Tom Duff et al. Journal of Computer Graphics Techniques Vol. 6 N.1, 2017
     */
    void orthonormal(Vector3& E, Vector3& F, real N) const
    {
#if 0
        if ( abs_real(1.0 - normSqr()) > 0.01 )
        {
            // this should not happen...
            E = orthogonal(1);
            F = cross(*this, E).normalized();
            std::clog << "rescued orthonormal(" << to_string() << ")\n";
            return;
        }
#else
        assert_small(1.0 - normSqr());
#endif
        real s = std::copysign(real(1.0), ZZ);
        // optimized version by Marc B. Reynolds
        real nY = N * YY;
        real a = nY / ( ZZ + s );
        real b = YY * a;
        real c = XX * a;
        // below normSqr(F) = normSqr(this) + a*a*(normSqr(this)-s*s)
        E.set(-N * ZZ - b, c, N * XX);
        F.set(s * c, s * b - N, s * nY);
        // if you do not mind an inverted basis, use F.set(c, b - s * N, nY);
        //printf("orthonormal %+9.6f %+9.6f %+9.6f :", dot(*this, E), dot(*this, F), dot(E, F));
        //printf(" %+9.6f %+9.6f %+9.6f\n", normSqr(), dot(E, E), dot(F, F));
    }

    /**
     Set 'E' and 'F' to build a basis (this, E, F), with norm(E) = norm(F) = N
     assuming that 'norm(*this) == N'
     
     Derived from `Building an Orthonormal Basis, Revisited`,
     Tom Duff et al. Journal of Computer Graphics Techniques Vol. 6 N.1, 2017
     */
    void orthonormal(real N, Vector3& E, Vector3& F) const
    {
        assert_small(norm()/N - 1.0);
        real s = std::copysign(real(1.0), ZZ);
        // optimized version by Marc B. Reynolds
        real a = YY / ( ZZ + s * N );
        real b = YY * a;
        real c = XX * a;
        // below normSqr(F) = normSqr(this) + a*a*(normSqr(this)-s*s)
        E.set(-ZZ-b, c, XX);
        F.set(s*c, s*b-N, s*YY);
        // if you do not mind an inverted basis, use F.set(c, b - s*N, YY);
    }


    /**
     Set 'E' and 'F' to build a basis (this, E, F), with norm(E) = norm(F) = NoL*L
     assuming that 'norm(*this) == L'
     
     Derived from `Building an Orthonormal Basis, Revisited`,
     Tom Duff et al. Journal of Computer Graphics Techniques Vol. 6 N.1, 2017
   */
    void orthonormal(real L, Vector3& E, Vector3& F, real NoL) const
    {
        assert_small(norm()/L - 1.0);
        real s = std::copysign(real(1.0), ZZ);
        // optimized version by Marc B. Reynolds
        real nY = NoL * YY;
        real a = nY / ( ZZ + s * L );
        real b = YY * a;
        real c = XX * a;
        // below normSqr(F) = normSqr(this) + a*a*(normSqr(this)-s*s)
        E.set(-NoL*ZZ-b, c, NoL*XX);
        F.set(s*c, s*b-L*NoL, s*nY);
        // if you do not mind an inverted basis, use F.set(c, b - s*N, YY);
    }
    
    /// rotate `vec` around `*this`, by angle defined by cosine and sine
    /**
     It is assumed that norm(*this)==1
     The result is a Vector orthogonal to *this, of norm `C^2 + S^2`
     */
    Vector3 rotateOrtho(Vector3 const& vec, real C, real S)
    {
        // set two orthogonal vector to 'd' defining an orientated basis
        Vector3 ex, ey;
        orthonormal(ex, ey);
        // compute coordinates of `vec` in the reference frame (ex, ey):
        real x = dot(vec, ex);
        real y = dot(vec, ey);
        // normalization factor:
        real n = 1.0 / std::sqrt( x * x + y * y );
        x = x * n;
        y = y * n;
        // rotated vector:
        return ( C * x - S * y ) * ex + ( S * x + C * y ) * ey;
    }
    
    /// convert from cartesian to spherical coordinates ( r, theta, phi )
    Vector3 spherical() const
    {
        return Vector3(std::sqrt(XX*XX+YY*YY+ZZ*ZZ),
                       std::atan2(YY, XX),
                       std::atan2(std::sqrt(XX*XX+YY*YY), ZZ));
    }
    
    /// convert from spherical to cartesian coordinates ( x, y, z )
    Vector3 cartesian() const
    {
        return Vector3(XX*std::cos(YY)*std::sin(ZZ),
                       XX*std::sin(YY)*std::sin(ZZ),
                       XX*std::cos(ZZ));
    }
    
    //------------------------------------------------------------------
        
    /// Calculate intermediate position = A + C * ( B - A )
    void interpolate(const float a[], const float C, const float b[])
    {
        XX = a[0] + C * ( b[0] - a[0] );
        YY = a[1] + C * ( b[1] - a[1] );
        ZZ = a[2] + C * ( b[2] - a[2] );
#if VECTOR3_USES_AVX
        vec[3] = 0;
#endif
    }

    /// Calculate intermediate position = A + C * ( B - A )
    void interpolate(const double a[], const double C, const double b[])
    {
#if VECTOR3_USES_AVX
        vec4 A = load3(a), B = load3(b);
        vec = fmadd4(set4(C), sub4(B, A), A);
#else
        XX = a[0] + C * ( b[0] - a[0] );
        YY = a[1] + C * ( b[1] - a[1] );
        ZZ = a[2] + C * ( b[2] - a[2] );
#endif
    }
    
    /// linear interpolation, returning a + x * b
    friend Vector3 interpolated(const Vector3& a, real C, const Vector3& b)
    {
        return Vector3(a.XX+C*b.XX, a.YY+C*b.YY, a.ZZ+C*b.ZZ);
    }

    /// Calculate intermediate position = A + C * ( B - A )
    static Vector3 interpolated(const float a[], const float C, const float b[])
    {
        return Vector3(a[0]+C*(b[0]-a[0]), a[1]+C*(b[1]-a[1]), a[2]+C*(b[2]-a[2]));
    }

    /// Calculate intermediate position = A + C * ( B - A )
    static Vector3 interpolated(const double a[], const double C, const double b[])
    {
#if VECTOR3_USES_AVX
        vec4 A = load3(a), B = load3(b);
        return Vector3(fmadd4(set4(C), sub4(B, A), A));
#else
        return Vector3(a[0]+C*(b[0]-a[0]), a[1]+C*(b[1]-a[1]), a[2]+C*(b[2]-a[2]));
#endif
    }

    //------------------------------------------------------------------

    /// addition of two vectors
    friend Vector3 operator +(Vector3 const& a, Vector3 const& b)
    {
#if VECTOR3_USES_AVX
        return Vector3(add4(a.vec, b.vec));
#else
        return Vector3(a.XX+b.XX, a.YY+b.YY, a.ZZ+b.ZZ);
#endif
    }
    
    /// subtraction of two vectors
    friend Vector3 operator -(Vector3 const& a, Vector3 const& b)
    {
#if VECTOR3_USES_AVX
        return Vector3(sub4(a.vec, b.vec));
#else
        return Vector3(a.XX-b.XX, a.YY-b.YY, a.ZZ-b.ZZ);
#endif
    }
    
    /// unary + operator does nothing
    friend Vector3 operator +(Vector3 const& b)
    {
        return b;
    }
    
    /// opposition of a vector
    friend Vector3 operator -(Vector3 const& b)
    {
#if VECTOR3_USES_AVX
        return Vector3(flipsign4(b.vec));
#else
        return Vector3(-b.XX, -b.YY, -b.ZZ);
#endif
    }
    
    /// returns the element-by-element product
    Vector3 e_mul(Vector3 const& b) const
    {
#if VECTOR3_USES_AVX
        return Vector3(mul4(vec, b.vec));
#else
        return Vector3(XX*b.XX, YY*b.YY, ZZ*b.ZZ);
#endif
    }

    /// returns the element-by-element division
    Vector3 e_div(Vector3 const& b) const
    {
        return Vector3(XX/b.XX, YY/b.YY, ZZ/b.ZZ);
    }
    
    /// returns a vector with each element squared
    Vector3 e_squared() const
    {
#if VECTOR3_USES_AVX
        return Vector3(mul4(vec, vec));
#else
        return Vector3(XX*XX, YY*YY, ZZ*ZZ);
#endif
    }
    
    /// returns sum of all coordinates
    real e_sum() const
    {
        return XX + YY + ZZ;
    }
    
    /// returns min(x, y, z)
    real e_min() const
    {
        return std::min(std::min(XX, YY), ZZ);
    }
    
    /// returns max(x, y, z)
    real e_max() const
    {
        return std::max(std::max(XX, YY), ZZ);
    }
    
    /// returns the element-by-element minimum
    Vector3 e_min(Vector3 const& v) const
    {
        return Vector3(std::min(XX, v.XX), std::min(YY, v.YY), std::min(ZZ, v.ZZ));
    }
    
    /// returns the element-by-element maximum
    Vector3 e_max(Vector3 const& v) const
    {
        return Vector3(std::max(XX, v.XX), std::max(YY, v.YY), std::max(ZZ, v.ZZ));
    }
    

    /// cross product of two vectors
    friend Vector3 cross(Vector3 const& a, Vector3 const& b)
    {
#if VECTOR3_USES_AVX && defined(__AVX2__)
        assert_true((a.vec[3] == 0) & (b.vec[3] == 0));
        return Vector3(cross4(a.vec, b.vec));
#else
        return Vector3(a.YY * b.ZZ - a.ZZ * b.YY,
                       a.ZZ * b.XX - a.XX * b.ZZ,
                       a.XX * b.YY - a.YY * b.XX);
#endif
    }

    /// scalar product of two vectors
    friend real dot(Vector3 const& a, Vector3 const& b)
    {
#if 0 //VECTOR3_USES_AVX
        vec2 axy = load2(&a.XX), bxy = load2(&b.XX);
        vec2 az = load1(&a.ZZ), bz = load1(&b.ZZ);
        vec2 d = fmadd2(az, bz, mul2(axy, bxy));
        return add2(d, permute2(d))[0];
#else
        return a.XX * b.XX + a.YY * b.YY + a.ZZ * b.ZZ;
#endif
    }
    
    /// multiplication by scalar
    friend Vector3 operator *(Vector3 const& a, const real s)
    {
#if VECTOR3_USES_AVX
        return Vector3(mul4(a.vec, set4(s)));
#else
        return Vector3(s*a.XX, s*a.YY, s*a.ZZ);
#endif
    }
    
    /// mutiplication by scalar
    friend Vector3 operator *(const real s, Vector3 const& a)
    {
#if VECTOR3_USES_AVX
        return Vector3(mul4(set4(s), a.vec));
#else
        return Vector3(s*a.XX, s*a.YY, s*a.ZZ);
#endif
    }
    
    /// division by scalar
    friend Vector3 operator /(Vector3 const& a, const real s)
    {
#if VECTOR3_USES_AVX
        return Vector3(div4(a.vec, set4(s)));
#else
        return Vector3(a.XX/s, a.YY/s, a.ZZ/s);
#endif
    }
    
    /// addition of another vector
    void operator +=(Vector3 const& b)
    {
#if VECTOR3_USES_AVX
        vec = add4(vec, b.vec);
#else
        XX += b.XX;
        YY += b.YY;
        ZZ += b.ZZ;
#endif
    }
    
    /// subtraction of another vector
    void operator -=(Vector3 const& b)
    {
#if VECTOR3_USES_AVX
        vec = sub4(vec, b.vec);
#else
        XX -= b.XX;
        YY -= b.YY;
        ZZ -= b.ZZ;
#endif
    }
    
    /// multiplication by a scalar
    void operator *=(const real s)
    {
#if VECTOR3_USES_AVX
        vec = mul4(vec, set4(s));
#else
        XX *= s;
        YY *= s;
        ZZ *= s;
#endif
    }
    
    /// division by a scalar
    void operator /=(const real s)
    {
#if VECTOR3_USES_AVX
        vec = div4(vec, set4(s));
#else
        XX /= s;
        YY /= s;
        ZZ /= s;
#endif
    }
    
    //------------------------------------------------------------------
    
    /// equality test
    friend bool operator ==(Vector3 const& a, Vector3 const& b)
    {
        return ( a.XX==b.XX  &&  a.YY==b.YY  &&  a.ZZ==b.ZZ );
    }
    
    /// non-equality test
    friend bool operator !=(Vector3 const& a, Vector3 const& b)
    {
        return ( a.XX!=b.XX  ||  a.YY!=b.YY  ||  a.ZZ!=b.ZZ );
    }
    
    //------------------------------------------------------------------
    
    /// output
    void print(std::ostream& os) const
    {
        const int w = (int)os.width();
        os << std::showpos << XX << " ";
        os << std::setw(w) << YY << " ";
        os << std::setw(w) << ZZ << std::noshowpos;
    }

    /// output using width 'w' and precision 'p'
    void print(std::ostream& os, int w, int p) const
    {
        os.precision(p);
        os.setf(std::ios::showpos);
        os << std::setw(w) << XX << " ";
        os << std::setw(w) << YY << " ";
        os << std::setw(w) << ZZ;
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
        fprintf(out, "  %+9.3f %+9.3f %+9.3f", XX, YY, ZZ);
    }
    
    /// print to a file, surrounded by parenthesis
    void pprint(FILE * out = stdout) const
    {
        fprintf(out, "( %+9.3f %+9.3f %+9.3f )", XX, YY, ZZ);
    }
    
    /// print, followed by a new line
    void println(FILE * out = stdout) const
    {
        fprintf(out, "  %+9.3f %+9.3f %+9.3f\n", XX, YY, ZZ);
    }
    
    //------------------------------------------------------------------
    
    /// add a random component in [-s, s] to each coordinate
    void addRand(real s);
    
    
    /// a vector of norm n, orthogonal to *this, assuming `norm(*this)==1`
    Vector3 randOrthoU(real n) const;
    
    /// a vector of norm <= n, orthogonal to *this, assuming `norm(*this)==1`
    Vector3 randOrthoB(real n) const;
    
    
    /// Vector with random independent coordinates in [0,+1]
    static Vector3 randP();
    
    /// Vector with random independent coordinates in [0,+n]
    static Vector3 randP(real n);
    
    /// Vector with random independent coordinates in [-1,+1]
    static Vector3 randS();
    
    /// Vector with random independent coordinates in [-1/2,+1/2]
    static Vector3 randH();
    
    /// Vector with random independent coordinates in [-n,+n]
    static Vector3 randS(real n);
    
    
    /// random Vector of norm = 1; sampling is uniform
    static Vector3 randU();
    
    /// return a random vector of norm = n; sampling is uniform
    static Vector3 randU(real n);
    
    
    /// return a random vector of norm <= 1; sampling is uniform
    static Vector3 randB();
    
    /// return a random vector of norm <= n; sampling is uniform
    static Vector3 randB(real n);
    
    
    /// return a random vector with Normally distributed coordinates ~ N(0,n)
    static Vector3 randG(real n);
    
};

//-------------------------- associated global functions -----------------------

/// stream input operator
std::istream& operator >> (std::istream&, Vector3&);

/// output operator
inline std::ostream& operator << (std::ostream& os, Vector3 const& arg)
{
    std::ios::fmtflags fgs = os.flags();
    arg.print(os);
    os.setf(fgs);
    return os;
}

#endif
