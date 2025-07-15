// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.
/*
 Testing Intel's Streaming SIMD on the vector dot product
 FJN, started August 2022
 
 To compile: c++ -O2 -mavx test.cc
 To generate assembly: c++ -S test.cc
 */

#include <cstdio>
#include "timer.h"
#include "simd.h"
#include "simd_float.h"

typedef double real;

const size_t CNT = 1<<14;

real vX[CNT], vY[CNT], vZ[CNT];

void init()
{
    for ( size_t i=0; i<CNT; ++i )
    {
        vX[i] = real(1)/real(CNT-i);
        vY[i] = real(CNT-i);
        vZ[i] = 1.0;
    }
}

// scalar version
real dot(const real* X, const real* Y)
{
    real d = 0;
    for ( size_t i = 0; i < CNT; ++i )
        d += X[i] * Y[i];
    return d;
}

#if defined(__SSE3__)
real dot_sse(const float* X, const float* Y)
{
    vec4f s = setzero4f();
    for ( size_t i = 0; i < CNT; i += 4 )
        s = add4f(s, mul4f( load4f(X+i), load4f(Y+i) ));
    _mm_empty();
    
    return s[0] + s[1] + s[2] + s[3];
}
#endif

#if USE_SIMD
real dot_SSE(const double* X, const double* Y)
{
    vec2 s = mul2(load2(X), load2(Y));
    for ( size_t i = 2; i < CNT; i += 2 )
        s = fmadd2(load2(X+i), load2(Y+i), s);
    _mm_empty();
    
    return s[0] + s[1];
}

real dot_SSEU(const double* X, const double* Y)
{
    vec2 s = mul2(load2(X), load2(Y));
    vec2 t = mul2(load2(X+2), load2(Y+2));
    for ( size_t i = 4; i < CNT; i += 4 )
    {
        s = fmadd2(load2(X+i), load2(Y+i), s);
        t = fmadd2(load2(X+i+2), load2(Y+i+2), t);
    }
    s = add2(s, t);
    _mm_empty();
    return s[0] + s[1];
}

real dot_SSEUU(const double* X, const double* Y)
{
    vec2 s = mul2(load2(X), load2(Y));
    vec2 t = mul2(load2(X+2), load2(Y+2));
    vec2 u = mul2(load2(X+4), load2(Y+4));
    vec2 v = mul2(load2(X+6), load2(Y+6));
    for ( size_t i = 8; i < CNT; i += 8 )
    {
        s = fmadd2(load2(X+i), load2(Y+i), s);
        t = fmadd2(load2(X+i+2), load2(Y+i+2), t);
        u = fmadd2(load2(X+i+4), load2(Y+i+4), u);
        v = fmadd2(load2(X+i+6), load2(Y+i+6), v);
    }
    s = add2(add2(s, t), add2(u, v));
    _mm_empty();
    return s[0] + s[1];
}
#endif

#ifdef __AVX__
real dot_AVX(const float* X, const float* Y)
{
    vec8f s = setzero8f();
    for ( size_t i = 0; i < CNT; i += 8 )
        s = add8f(s, mul8f( load8f(X+i), load8f(Y+i) ));
    _mm_empty();
    
    return s[0] + s[1] + s[2] + s[3];
}

real dot_AVX(const double* X, const double* Y)
{
    vec4 s = mul4(load4(X), load4(Y));
    for ( size_t i = 4; i < CNT; i += 4 )
        s = fmadd4(load4(X+i), load4(Y+i), s);
    _mm_empty();
    
    return s[0] + s[1] + s[2] + s[3];
}

real dot_AVXU(const double* X, const double* Y)
{
    vec4 v0 = mul4(load4(X), load4(Y));
    vec4 v1 = mul4(load4(X+4), load4(Y+4));
    vec4 v2 = mul4(load4(X+8), load4(Y+8));
    vec4 v3 = mul4(load4(X+12), load4(Y+12));
    
    for ( size_t i = 16; i < CNT; i += 16 )
    {
        v0 = fmadd4(load4(X+i   ), load4(Y+i   ), v0);
        v1 = fmadd4(load4(X+i+4 ), load4(Y+i+4 ), v1);
        v2 = fmadd4(load4(X+i+8 ), load4(Y+i+8 ), v2);
        v3 = fmadd4(load4(X+i+12), load4(Y+i+12), v3);
    }
    
    vec4 s = add4(add4(v0, v1), add4(v2, v3));
    _mm_empty();
    
    return esum4(s)[0];
}
#endif

#ifdef __ARM_NEON__
#include <arm_neon.h>

real dot_neon(const float* X, const float* Y)
{
    float32x4_t sum = { 0 };
    #pragma clang loop interleave_count(2)
    for ( int i = 0; i < CNT; i += 4 )
        sum = vfmaq_f32(sum, *(float32x4_t*)(X+i), *(float32x4_t*)(Y+i));
    return sum[0] + sum[1] + sum[2] + sum[3];
}

real dot_neon(const double* X, const double* Y)
{
    float64x2_t sum = { 0 };
    #pragma clang loop interleave_count(2)
    for ( int i = 0; i < CNT; i += 2 )
        sum = vfmaq_f64(sum, *(float64x2_t*)(X+i), *(float64x2_t*)(Y+i));
    return sum[0] + sum[1];
}

#endif


void run(real (*func)(const real*, const real*), const char name[], const size_t REP)
{
    real a = 0, b = 0, c = 0, d = 0;
    real e = 0, f = 0, g = 0, h = 0;
    fprintf(stderr, "%10s : ", name);
    tick();
    
    for ( size_t i=0; i<REP; ++i )
    {
        a = (*func)(vX, vY);
        b = (*func)(vY, vZ);
        c = (*func)(vX, vZ);
        d = (*func)(vX, vY);
        e = (*func)(vY, vZ);
        f = (*func)(vZ, vX);
        g = (*func)(vZ, vY);
        h = (*func)(vX, vY);
    }
    
    real s = ((a + b) + (c + d)) + ((e + f) + (g + h));
    fprintf(stderr, " %16f :  %8.0f ms\n", s, tock(8));
}


int main(int argc, char * argv[])
{
    const size_t REP = 1<<12;
    init();
    run(dot,  "scalar", REP);
#if USE_SIMD
    run(dot_SSE, "SSE", REP);
    run(dot_SSEU, "SSEU", REP);
    run(dot_SSEUU, "SSEUU", REP);
#endif
#ifdef __AVX__
    run(dot_AVX, "AVX", REP);
    run(dot_AVXU, "AVXU", REP);
#endif
#ifdef __ARM_NEON__
    run(dot_neon, "neon", REP);
#endif
}

