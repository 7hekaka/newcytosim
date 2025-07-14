// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef CYMMATH_H
#define CYMMATH_H

#include "real.h"
#include <cmath>
#include <sstream>
#include <cstdint>
#include <iomanip>


#ifndef M_PI
/// Ratio of a circle's circumference to its diameter
constexpr real M_PI=3.14159265358979323846264338327950288;
#endif

#ifndef M_E
constexpr real M_E=2.7182818284590452354;
#endif


/// simple mathematical functions, mostly templated
namespace Cymath
{
    /// minimum of three arguments
    template <typename T>
    inline const T& min(const T& a, const T& b, const T& c)
    {
        return std::min(a, std::min(b, c));
    }
    
    /// maximum of three arguments
    template <typename T>
    inline const T& max(const T& a, const T& b, const T& c)
    {
        return std::max(a, std::max(b, c));
    }
    
    /// minimum of four arguments
    template <typename T>
    inline const T& min(const T& a, const T& b, const T& c, const T& d)
    {
        return std::min(std::min(a,b), std::min(c,d));
    }
    
    /// maximum of four arguments
    template <typename T>
    inline const T& max(const T& a, const T& b, const T& c, const T& d)
    {
        return std::max(std::max(a,b), std::max(c,d));
    }
    
    /// sort in ascending order using min() and max() functions
    template <typename T>
    inline void mm_sort(T& a, T& b)
    {
        T i = a;
        a = std::min(a, b);
        b = std::max(i, b);
    }
    
    /// sort in ascending order
    template <typename T>
    inline void mm_sort(T& a, T& b, T& c)
    {
        mm_sort(a, b);
        mm_sort(b, c);
        mm_sort(a, b);
    }

    /// sort in ascending order
    template <typename T>
    inline void mm_sort(T& a, T& b, T& c, T& d)
    {
        //4 inputs [[1 2][3 4][1 3][2 4][2 3]]
        mm_sort(a, b);
        mm_sort(c, d);
        mm_sort(a, c);
        mm_sort(b, d);
        mm_sort(b, c);
    }

    /// sort in ascending order
    template <typename T>
    inline void mm_sort(T& a, T& b, T& c, T& d, T& e)
    {
        //5 inputs [[1 2][3 4][1 3][2 5][1 2][3 4][2 3][4 5][3 4]]
        mm_sort(a, b);
        mm_sort(c, d);
        mm_sort(a, c);
        mm_sort(b, e);
        mm_sort(a, b);
        mm_sort(c, d);
        mm_sort(b, c);
        mm_sort(d, e);
        mm_sort(c, d);
    }

    /// sort in ascending order
    template <typename T>
    inline void mm_sort(T& a, T& b, T& c, T& d, T& e, T& f)
    {
        //6 inputs [[1 2][3 4][5 6][1 3][2 5][4 6][1 2][3 4][5 6][2 3][4 5][3 4]]
        mm_sort(a, b);
        mm_sort(c, d);
        mm_sort(e, f);
        mm_sort(a, c);
        mm_sort(b, e);
        mm_sort(d, f);
        mm_sort(a, b);
        mm_sort(c, d);
        mm_sort(e, f);
        mm_sort(b, c);
        mm_sort(d, e);
        mm_sort(c, d);
    }

    /// return index of the arguments that is the smallest: {0, 1, 2}
    template <typename T>
    inline int arg_min(const T& a, const T& b, const T& c)
    {
        return 2*int( c < std::min(b,a) ) | int( b < std::min(a,c) );
    }
    
    /// return index of the arguments that is the largest: {0, 1, 2}
    template <typename T>
    inline int arg_max(const T& a, const T& b, const T& c)
    {
        return 2*int( c > std::max(b,a) ) | int( b > std::max(a,c) );
    }

    /// return index of the arguments that is the smallest: {0, 1, 2, 3}
    template <typename T>
    inline int arg_min(const T& a, const T& b, const T& c, const T& d)
    {
        T ab = std::min(a, b);
        T cd = std::min(c, d);
        return 3*int( d < std::min(c,ab) ) | 2*int( c < std::min(d,ab) ) | int( b < std::min(a,cd) );
    }
    
    /// return index of the arguments that is the largest: {0, 1, 2, 3}
    template <typename T>
    inline int arg_max(const T& a, const T& b, const T& c, const T& d)
    {
        T ab = std::max(a, b);
        T cd = std::max(c, d);
        return 3*int( d > std::max(c,ab) ) | 2*int( c > std::max(d,ab) ) | int( b > std::max(a,cd) );
    }


    /**
     Set vectors 'x' and 'y' to make an orthonormal basis (x, y, z)
     
     Building an Orthonormal Basis, Revisited
     Tom Duff et al. Journal of Computer Graphics Techniques Vol. 6 N.1, 2017
     */
    template < typename FLOAT >
    inline void orthonormal(const FLOAT z[3], FLOAT x[3], FLOAT y[3])
    {
        const FLOAT s = std::copysign((FLOAT)1, z[2]);
#if ( 1 )
        // optimized version by Marc B. Reynolds
        const FLOAT a = z[1] / ( z[2] + s );
        const FLOAT b = z[1] * a;
        const FLOAT c = z[0] * a;
        x[0] = -z[2] - b;
        x[1] = c;
        x[2] = z[0];
        y[0] = s * c;
        y[1] = s * b - 1;
        y[2] = s * z[1];
#else
        // original code from Tom Duff et al.
        const FLOAT a = -1 / ( a[2] + s );
        const FLOAT b = a[0] * a[1] * a;
        x[0] = 1 + s * z[0] * z[0] * a;
        x[1] = s * b;
        x[2] = -s * z[0];
        y[0] = b;
        y[1] = s + z[1] * z[1] * a;
        y[2] = -z[1];
#endif
    }
    
    /// square of a number
    template <typename T> 
    inline T square(const T& a)
    {
        return a * a;
    }
    
    /// cube of a number
    template <typename T>
    inline T cube(const T& a)
    {
        return a * a * a;
    }
    
    /// power of `a` by positive integer exponent `n`
    /** This should be equivalent to std::pow(a, n) */
    template <typename T>
    inline T power_int(const T& a, unsigned n)
    {
        T x = a;
        T y = 1;
        while ( n )
        {
            if ( n & 1 )
                y = y * x;
            x = x * x;
            n = n >> 1;
        }
        return y;
    }
    
    ///power of `a` by signed integer exponent `n`
    template <typename T>
    inline T power(const T& a, const int n)
    {
        if ( n < 0 )
            return power_int(1.0/a, -n);
        return power_int(a, n);
    }
    
    
    ///square of distance between two vectors in dimension `dim`
    template <int dim, typename T>
    inline T distanceSqr(const T a[], const T b[])
    {
        T x = a[0] - b[0];
        T n = x * x;
        for( int i = 1; i < dim; ++i )
        {
            x = a[i] - b[i];
            n += x * x;
        }
        return n;
    }
    
    ///usual distance between two vectors of dimension `dim`
    template <int dim, typename T>
    inline T distance(const T a[], const T b[])
    {
        T x = a[0] - b[0];
        T n = x * x;
        for( int i = 1; i < dim; ++i )
        {
            x = a[i] - b[i];
            n += x * x;
        }
        return std::sqrt(n);
    }

    /// return the usual base-10 representation of number `x`
    template <typename T>
    std::string repr(T const& x, int width, unsigned precision)
    {
        std::ostringstream oss;
        oss.precision(precision);
        oss << std::setw(width) << std::fixed << x;
        return oss.str();
    }
    
    //------------------------------------------------------------------------------
#pragma mark - Function used for periodic boundary conditions
    /*
    /// used for periodic boundary conditions:
    inline void fold(real& x, const real p)
    {
        while ( x >  p ) x -= p+p;
        while ( x < -p ) x += p+p;
    }
*/
#ifdef WIN32
    
    //this is needed under windows:
    inline real remainder(const real a, const real b)
    {
        real p = std::floor( 0.5 + a / b );
        return a - p * b;
    }
    
    inline real round(real x)
    {
        if ( x < 0 )
            return -std::floor(0.5-x);
        else
            return  std::floor(0.5+x);
    }

#endif
    
    //------------------------------------------------------------------------------
#pragma mark -
    
    ///extract a 10-decimal digit form a number:
    /** 1st digit is really the first one, as index here do not start at zero */
    template <typename T> 
    inline int digit(T x, const int p)
    {
        for ( int q=1; q < p; ++q )
            x /= 10;
        return x % 10;
    }
    
    ///copy bytes
    inline void copyBytes( void * dest, const void * src, const unsigned cnt)
    {
        for ( size_t i = 0; i < cnt; ++i )
            ((char*)dest)[i] = ((char*)src)[i];
    }
    

    //------------------------------------------------------------------------------
    
    /// return smallest power of 2 that is greater or equal to `x`
    inline unsigned next_power(unsigned x)
    {
        if ( x > 0 )
        {
            --x;
            x |= x >> 1;
            x |= x >> 2;
            x |= x >> 4;
            x |= x >> 8;
            x |= x >> 16;
        }
        return x+1;
    }
    
    /// return smallest power of 2 that is greater or equal to `x`
    inline size_t next_power(size_t x)
    {
        if ( x > 0 )
        {
            --x;
            x |= x >> 1;
            x |= x >> 2;
            x |= x >> 4;
            x |= x >> 8;
            x |= x >> 16;
            x |= x >> 32;
        }
        return x+1;
    }

    /// number of '1' bits in a 32-bits integer (Charlie Gordon & Don Clugston)
    /** Should use Intel SIMD instruction POPCNT */
    inline unsigned int count_bits(uint32_t v)
    {
        v = v - ((v >> 1) & 0x55555555);
        v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
        return (((v + (v >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;
    }
    
    
    /// number of '1' bits, from: http://graphics.stanford.edu/~seander/bithacks.html
    /** Works up to 128 bits */
    template <typename T>
    unsigned int count_bits2(T v)
    {
        v = v - ((v >> 1) & (T)~(T)0/3);
        v = (v & (T)~(T)0/15*3) + ((v >> 2) & (T)~(T)0/15*3);
        v = (v + (v >> 4)) & (T)~(T)0/255*15;
        return (T)(v * ((T)~(T)0/255)) >> (sizeof(v) - 1) * 8;
    }
    
    
    /// set binary representation of integer 'val' with 0/1s, using 'wid' characters
    template < typename T >
    void binary_representation(char str[], size_t len, size_t wid, T val, char end = 0)
    {
        char * s = str;
        size_t i = 0;
        for ( ; val && i+1 < len; ++i )
        {
            *str++ = ( val & 1 ? '1' : '0' );
            val >>= 1;
        }
        for ( ; i+1 < len && i < wid; ++i )
            *str++ = '0';
        // swap order of characters in string:
        char t, * e = str - 1;
        while ( s < e )
        {
            t = *s;
            *s++ = *e;
            *e-- = t;
        }
        *str = end;
    }

    
    /* Find the nullpoint of a monotonously increasing function */
    inline real find_root(real (*func)(real), real a, real b)
    {
        if ( func(a) > 0 ) std::swap(a,b);
        for ( int i = 0; i < 16; ++i )
        {
            real m = 0.5 * ( a + b );
            real v = func(m);
            if ( v < 0 ) a = m;
            else b = m;
        }
        return 0.5 * ( a + b );
    }
}


#endif //#ifdef SMATH_H
