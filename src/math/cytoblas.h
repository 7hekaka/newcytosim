// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

/*
 * BLAS-style extensions for reals; should be included after 'blas.h'
 */

#ifndef CYTOBLAS_H
#define CYTOBLAS_H

#include <algorithm>
#include <cstdio>
#include <cmath>

#include "simd.h"
#include "simd_float.h"

namespace blas
{

/// this is the standard Euclidian norm
inline double nrm2(int N, const real* X)
{
    //using double precision to accumulate:
    return std::sqrt(blas::xdot(N, X, 1, X, 1));
}


#ifdef __INTEL_MKL__
/**
 axpby() performs Y <- alpha * X + beta * Y
 as provided by Intel Math Kernel Library
 */
void BLAS(axpby)(int*, real*, const real*, int*, real*, real*, int*);
inline void xaxpby(int N, real alpha, const real*X, int incX, real beta, real*Y, int incY)
{
    BLAS(axpby)(&N, &alpha, X, &incX, &beta, Y, &incY);
}
#else
/**
 axpby() performs Y <- alpha * X + beta * Y
 */
inline void xaxpby(int N, real alpha, const real* X, int incX, real beta, real* Y, int incY)
{
    if ( incX == 1  &&  incY == 1 )
    {
        for ( int i = 0; i < N; ++i )
            Y[i] = alpha * X[i] + beta * Y[i];
    }
    else
    {
        for ( int i = 0; i < N; ++i )
            Y[i*incY] = alpha * X[i*incX] + beta * Y[i*incY];
    }
}
#endif


/// calculates Y <- X + alpha * Y
inline void xpay(size_t N, const real* X, real alpha, real* Y)
{
    #pragma omp simd
    for ( size_t i = 0; i < N; ++i )
        Y[i] = alpha * Y[i] + X[i];
}


/// addition Y[] <- Y[] + X[], for array of size N
inline void xadd(size_t N, const real* X, real* Y)
{
    //xaxpy(N, 1.0, X, 1, Y, 1);
    #pragma omp simd
    for ( size_t i = 0; i < N; ++i )
        Y[i] = Y[i] + X[i];
}
    
/// subtraction Y[] <- Y[] - X[], for array of size N
inline void xsub(size_t N, const real* X, real* Y)
{
    //xaxpy(N, -1.0, X, 1, Y, 1);
    #pragma omp simd
    for ( size_t i = 0; i < N; ++i )
        Y[i] = Y[i] - X[i];
}


/**
 return the infinite norm of the vector

     int inx = ixamax(N, X, inc);
     return abs(X[inx]);
 
 */
inline real nrm8(const int N, const real* X, int inc)
{
#if ( 1 )
    size_t inx = blas::ixamax(N, X, inc);
    return abs_real(X[inx-1]);
#else
    if ( N == 0 )
        return 0;
    real u = abs_real(X[0]);
    for ( int i = 1; i < N; ++i )
        u = std::max(u, abs_real(X[i*inc]));
    return u;
#endif
}

    
inline real nrm8seq(const size_t cnt, const real* X)
{
    real res = abs_real(X[0]);
    #pragma omp simd
    for ( size_t i = 1; i < cnt; ++i )
        res = std::max(res, abs_real(X[i]));
    return res;
}

#if defined(__AVX__) && USE_SIMD

inline double nrm8(const size_t cnt, const double* ptr)
{
    //double const* adr = ptr;
    double const* end = ptr + cnt;
    const vec4 sign = {-0.0, -0.0, -0.0, -0.0};
    vec4 u = setzero4();
    #pragma nounroll
    while ( ptr+12 <= end )
    {
        vec4 a = andnot4(sign, load4(ptr));
        vec4 b = andnot4(sign, load4(ptr+4));
        vec4 c = andnot4(sign, load4(ptr+8));
        u = max4(max4(u,a), max4(b,c));
        ptr += 12;
    }
    #pragma nounroll
    while ( ptr+8 <= end )
    {
        u = max4(u, andnot4(sign, load4(ptr)));
        ptr += 4;
    }
    vec2 w = max2(getlo(u), gethi(u));
    while ( ptr+2 <= end )
    {
        w = max2(w, andnot2(getlo(sign), load2(ptr)));
        ptr += 2;
    }
    double res = std::max(w[0], w[1]);
    for ( ; ptr < end; ++ptr )
        res = std::max(res, std::fabs(*ptr));
#if 0
    real z = std::fabs(adr[0]);
    for ( size_t i = 1; i < cnt; ++i )
        z = std::max(z, std::fabs(adr[i]));
    if ( z != res )
        printf("blas::nrm8 size %lu ERROR %f %f\n", cnt, z, res);
#endif
    return res;
}

inline float nrm8(const size_t siz, const float* ptr)
{
    float const* end = ptr + siz;
    vec8f u = setzero8f();
    while ( ptr+24 <= end )
    {
        vec8f a = abs8f(load8f(ptr));
        vec8f b = abs8f(load8f(ptr+8));
        vec8f c = abs8f(load8f(ptr+16));
        u = max8f(max8f(u,a), max8f(b,c));
        ptr += 24;
    }
    while ( ptr+8 <= end )
    {
        u = max8f(u, abs8f(load8f(ptr)));
        ptr += 8;
    }
    vec4f v = max4f(getlo4f(u), gethi4f(u));
    while ( ptr+4 <= end )
    {
        v = max4f(v, abs4f(load4f(ptr)));
        ptr += 4;
    }
    v = max4f(v, permute4f(v, 0b01));
    float res = v[0];
    while ( ptr < end )
        res = std::max(res, std::fabs(*ptr++));
    return res;
}

#elif defined(__ARM_NEON__) && USE_SIMD

inline double nrm8(const size_t cnt, const double* ptr)
{
    vec2 var{0.0, 0.0};
    double const* end = ptr + cnt;
    double const* stop = end - 2;
    for ( ; ptr < stop; ptr += 2 )
        var = max2(var, abs2(load2(ptr)));
    double res = std::max(var[0], var[1]);
    for ( ; ptr < end; ++ptr )
        res = std::max(res, std::fabs(*ptr));
    return res;
}

inline float nrm8(const size_t cnt, const float* ptr)
{
    vec4f var{0.0, 0.0, 0.0, 0.0};
    float const* end = ptr + cnt;
    float const* stop = end - 4;
    for ( ; ptr < stop; ptr += 4 )
        var = max4f(var, abs4f(load4f(ptr)));
    float res = std::max(std::max(var[0], var[1]), std::max(var[2], var[3]));
    for ( ; ptr < end; ++ptr )
        res = std::max(res, std::fabs(*ptr));
    return res;
}

#else
inline real nrm8(const size_t N, const real* X)
{
    real r = 0;
    #pragma omp simd
    for ( size_t i = 0; i < N; ++i )
        r = std::max(r, abs_real(X[i]));
    return r;
}
#endif

    
/**
 return the infinite norm of the difference between two vectors
 */
inline real difference(const size_t N, const real* X, const real* Y)
{
    if ( N == 0 )
        return 0;
    real u = abs_real(X[0]-Y[0]);
    for ( size_t i = 1; i < N; ++i )
        u = std::max(u, abs_real(X[i]-Y[i]));
    return u;
}


/**
Set N values of `X` to value `alpha`
 */
inline void xfill(const int N, real alpha, real* X)
{
    for ( int u = 0; u < N; ++u )
        X[u] = alpha;
}

/**
 Set N values of `X` to value `alpha`
*/
inline void xfill(const int N, real alpha, real* X, const int inc)
{
    for ( int u = 0; u < N*inc; u += inc )
        X[u] = alpha;
}

}


#endif
