// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University
/**
 @file real.h
 SYNOPSIS: we define and use "real" to be able to easily change
 the floating point precision, depending on the application.
 REAL_EPSILON is a lower limit of the precision achieved.
 */

#ifndef REAL_H
#define REAL_H


#include <cmath>
#include <cfloat>
#include <cstring>  // memset
#include <algorithm>
#include <new>

/**
 It is possible to select double or single precision throughout Cytosim
 by editing this file.
 
 Calculations might be faster in single precision, but the iterative solver used
 in Meca::solve() (the conjugate-gradient method) may fail in adverse conditions.
 Much of the code was optimized for double precision but not for single precision.
 
 It is safer, and STRONGLY ADVISED therefore, to use double precision!
*/
#define REAL_IS_DOUBLE 1


#if REAL_IS_DOUBLE
   /// real is an alias to double
   typedef double real;
   constexpr real REAL_EPSILON = 128 * DBL_EPSILON;
#else
   /// real is an alias to float
   typedef float real;
   constexpr real REAL_EPSILON = 128 * FLT_EPSILON;
#endif


//----------------------------ALLOCATION----------------------------------------

/// add consistency checks to new_real() and free_real()
#define CHECK_ALLOCATIONS 0


#if CHECK_ALLOCATIONS
#  include <set>
#  include <cstdio>
#  include "backtrace.h"
static std::set<void*> allocations;
#endif


/// return a number greater or equal to 's' that is a multiple of 4
inline static size_t chunk_real(size_t cnt)
{
    // align to 4 doubles (of size 8 bytes), hence 32 bytes
    constexpr size_t chunk = 32 / sizeof(real);
    // return a multiple of chunk greater than 's'
    // this bit trickery works if chunk is a pure power of 2
    return ( cnt + chunk - 1 ) & ~( chunk - 1 );
}


/// allocate a new array to hold `size` real scalars
/** The returned pointer is aligned to a 64 byte boundary */
inline static real* new_real(size_t cnt)
{
    void* ptr = nullptr;
    /*
     We need to align to 4 doubles (of size 8 bytes), hence 32 bytes
     Allocating to 64 bytes matches the cache boundary on most CPUs
     */
    if ( posix_memalign(&ptr, 64, cnt*sizeof(real)) )
        throw std::bad_alloc();
    real* res = (real*)ptr;
    //printf("%p = new_real(%lu)  %lu\n", ptr, cnt, ((uintptr_t)ptr&63));
#if CHECK_ALLOCATIONS
    allocations.insert(res);
#endif
#if ( 0 )
    /*
     Allocated memory can be filled with signalling NaN, to catch access
     to uninitialized data, using the option '-fp-trap-all=divzero,invalid'
     from the intel compiler, or with GCC:
     feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
     */
    real n = std::numeric_limits<real>::signaling_NaN();
    for ( size_t u = 0; u < cnt; ++u )
        res[u] = n;
#endif
    return res;
}


/// release an array of reals allocated by `new_real`
inline static void free_real(void* ptr)
{
#if CHECK_ALLOCATIONS
    if ( ptr )
    {
        auto i = allocations.find(ptr);
        if ( i == allocations.end() )
        {
            printf("unallocated free_real(%p)\n", ptr);
            print_backtrace();
        }
        else
            allocations.erase(i);
        //printf("free_real(%p)\n", ptr);
    }
#endif
    free(ptr);
}


/// copy `cnt` real scalars from `src` to `dst`
inline static void copy_real(size_t cnt, real const* src, real * dst)
{
#if ( 0 )
    memcpy(dst, src, cnt*sizeof(real));
#else
    //#pragma vector unaligned
    for ( size_t u = 0; u < cnt; ++u )
        dst[u] = src[u];
#endif
}

#if REAL_IS_DOUBLE
/// copy `cnt` real scalars from `src` to `dst`
inline static void copy_real(size_t cnt, real const* src, float * dst)
{
    for ( size_t u = 0; u < cnt; ++u )
        dst[u] = (float)src[u];
}
#else
/// copy `cnt` real scalars from `src` to `dst`
inline static void copy_real(size_t cnt, real const* src, double * dst)
{
    for ( size_t u = 0; u < cnt; ++u )
        dst[u] = (double)src[u];
}
#endif


/// set `cnt` values of the array `ptr` to 0 (zero).
inline static void zero_real(size_t cnt, real * ptr)
{
#if ( 1 )
    // this works because IEEE 754 '+0.0' is represented with all bits at zero
    memset(ptr, 0, cnt*sizeof(real));
#else
    #pragma ivdep
    for ( size_t u = 0; u < cnt; ++u )
        ptr[u] = 0.0;
#endif
}

//-------------------------------CONSTANTS--------------------------------------

#ifndef M_SQRT3
constexpr real M_SQRT3 = 1.7320508075688772935274463415059;
#endif

//----------------------------BRANCHLESS? CODE----------------------------------

/// square of the argument: `x * x`
inline static real square(const real x) { return x * x; }

/// cube of the argument: `x * x * x`
inline static real cube(const real x) { return x * x * x; }

/// return `neg` if `val < 0` and `pos` otherwise
inline static real sign_select(real const val, real const neg, real const pos)
{
    // this should be branchless, using a conditional-move instruction (CMOVBE)
    return ( val < 0 ? neg : pos );
}

/// sign of a 'real': -1 or +1; result is +1 if ( x == 0 ) and -1 if ( x == -0 )
inline static real sign_real(const real x)
{
#if REAL_IS_DOUBLE
    return std::copysign(1.0, x);
#else
    return std::copysign(1.0f, x);
#endif
}

/// absolute value of `x`
inline static real abs_real(const real x) { return std::fabs(x); }

/// minimum between `x` and `y`
inline static real min_real(const real x, const real y) { return std::min(x, y); }

/// maximum between `x` and `y`
inline static real max_real(const real x, const real y) { return std::max(x, y); }

/// clamp value 'x' within [i, s]
inline static real clamp_real(const real x, const real i, const real s)
{
    return std::max(i, std::min(x, s));
}

/// adjust 'x' to canonical image with period 'p':
inline static real fold_real(const real x, const real p)
{
    // using remainder() function for branchless code
    return std::remainder(x, p);
}

/// return max absolute difference of `a[i] - b[i]` for i in [0, cnt]
inline static real norm_inf(size_t cnt, real const* a, real const* b)
{
    real s = 0;
    for ( size_t u = 0; u < cnt; ++u )
        std::max(s, std::fabs( a[u] - b[u] ));
    return s;
}

//----------------------------------- DEBUG ------------------------------------

inline static size_t has_nan(size_t cnt, real const* ptr)
{
    size_t res = 0;
    for ( size_t i = 0; i < cnt; ++i )
        res += std::isnan(ptr[i]);
    return res;
}

#endif
