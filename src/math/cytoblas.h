// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

/*
 * BLAS-style non-standard mathematical functions
 */

#ifndef CYTOBLAS_H
#define CYTOBLAS_H

#include <algorithm>
#include <cstdio>
#include "simd.h"


namespace blas
{

#ifdef __INTEL_MKL__
    /**
     axpby() performs Y <- alpha * X + beta * Y
     This routine is not part of BLAS, but is provided by Intel Math Kernel Library
     */
    void BLAS(axpby)(int*, real*, const real*, int*, real*, real*, int*);
    inline void xaxpby(int N, real alpha, const real*X, int incX, real beta, real*Y, int incY)
    {
        BLAS(axpby)(&N, &alpha, X, &incX, &beta, Y, &incY);
    }
#else
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
    #pragma ivdep
    #pragma vector always
    for ( size_t i = 0; i < N; ++i )
        Y[i] = alpha * Y[i] + X[i];
}


/// addition Y[] <- Y[] + X[], for array of size N
inline void add(size_t N, const real* X, real* Y)
{
    //xaxpy(N, 1.0, X, 1, Y, 1);
    #pragma ivdep
    #pragma vector always
    for ( size_t i = 0; i < N; ++i )
        Y[i] = Y[i] + X[i];
}
    
/// subtraction Y[] <- Y[] - X[], for array of size N
inline void sub(size_t N, const real* X, real* Y)
{
    //xaxpy(N, -1.0, X, 1, Y, 1);
    #pragma ivdep
    #pragma vector always
    for ( size_t i = 0; i < N; ++i )
        Y[i] = Y[i] - X[i];
}


/**
 return the infinite norm of the vector

     int inx = ixamax(N, X, inc);
     return fabs(X[inx]);
 
 */
inline real nrm8(const size_t N, const real* X, int inc)
{
#if ( 1 )
    int inx = blas::ixamax(N, X, inc);
    return std::abs(X[inx-1]);
#else
    if ( N == 0 )
        return 0;
    real u = std::abs(X[0]);
    for ( size_t i = 1; i < N; ++i )
        u = std::max(u, std::abs(X[i*inc]));
    return u;
#endif
}

    
inline real nrm8seq(const size_t siz, const real* X)
{
    real res = std::abs(X[0]);
#pragma ivdep
#pragma vector always
    for ( size_t i = 1; i < siz; ++i )
        res = std::max(res, std::abs(X[i]));
    return res;
}

#ifdef __AVX__

inline double nrm8(const size_t siz, const double* X)
{
    double const* ptr = X;
    double const* end = X + siz;
    double const* stop = end - 11;
    vec4 u = setzero4();
    #pragma nounroll
    while ( ptr < stop )
    {
        vec4 a = abs4(load4(ptr));
        vec4 b = abs4(load4(ptr+4));
        vec4 c = abs4(load4(ptr+8));
        u = max4(max4(u,a), max4(b,c));
        ptr += 12;
    }
    #pragma nounroll
    while ( ptr < end - 3 )
    {
        u = max4(u, abs4(load4(ptr)));
        ptr += 4;
    }
    vec2 v = getlo(max4(u, permute2f128(u, u, 0x01)));
    while ( ptr < end - 1 )
    {
        v = max2(v, abs2(load2(ptr)));
        ptr += 2;
    }
    v = max2(v, permute2(v, 0b01));
    double res = v[0];
    while ( ptr < end )
        res = std::max(res, std::abs(*ptr++));
#if 0
    real x = std::abs(X[0]);
    for ( size_t i = 1; i < siz; ++i )
        x = std::max(x, std::abs(X[i]));
    if ( x != res )
        printf("ERROR blas::nrm8 %f %f\n", x, res);
#endif
    return res;
}

    
inline float nrm8(const size_t siz, const float* X)
{
    float const* ptr = X;
    float const* end = X + siz;
    float const* stop = end - 24;
    vec8f u = setzero8f();
    while ( ptr <= stop )
    {
        vec8f a = abs8f(load8f(ptr));
        vec8f b = abs8f(load8f(ptr+8));
        vec8f c = abs8f(load8f(ptr+16));
        u = max8f(max8f(u,a), max8f(b,c));
        ptr += 24;
    }
    while ( ptr <= end - 8 )
    {
        u = max8f(u, abs8f(load8f(ptr)));
        ptr += 8;
    }
    vec4f v = _mm256_castps256_ps128(max8f(u, _mm256_permute2f128_ps(u, u, 0x01)));
    while ( ptr <= end - 4 )
    {
        v = max4f(v, abs4f(load4f(ptr)));
        ptr += 4;
    }
    v = max4f(v, permute4f(v, 0b01));
    float res = v[0];
    while ( ptr < end )
        res = std::max(res, std::abs(*ptr++));
    return res;
}

#else
inline real nrm8(const size_t N, const real* X)
{
    real r = 0;
    #pragma ivdep
    #pragma vector always
    for ( size_t i = 0; i < N; ++i )
        r = std::max(r, std::abs(X[i]));
    return r;
}
#endif

    
/**
 return the infinite norm of the difference between two vectors
 */
inline real max_diff(const size_t N, const real* X, const real* Y)
{
    if ( N == 0 )
        return 0;
    real u = std::abs(X[0] - Y[0]);
    for ( size_t i = 1; i < N; ++i )
        u = std::max(u, std::abs(X[i] - Y[i]));
    return u;
}


/**
Set N values of `X` to value `alpha`
 */
inline void xfill(const size_t N, real alpha, real* X)
{
    for ( size_t u = 0; u < N; ++u )
        X[u] = alpha;
}

/**
 Set N values of `X` to value `alpha`
*/
inline void xfill(const size_t N, real alpha, real* X, const int inc)
{
    for ( size_t u = 0; u < N; ++u )
        X[u*inc] = alpha;
}

}


#endif
