// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#include "real.h"
#include "random.h"
#include <cstdio>
#include <bitset>
#include <cstring>
#include <iostream>
#include <climits>
#include "timer.h"
#include <random>

/**
 Fill array `vec[]` with Gaussian values ~ N(0,1).
 the size of `vec` should be a multiple of 2, and sufficient to hold `end-src` values
 @Return the number of values that were stored in `vec`
 */
real * makeGaussians_(real dst[], size_t cnt, const int32_t src[])
{
    int32_t const*const end = src + cnt;
    while ( src < end )
    {
        real x = src[0] * TWO_POWER_MINUS_31;
        real y = src[1] * TWO_POWER_MINUS_31;
#if 1
        if ( std::abs(x) + std::abs(y) >= M_SQRT2 )
        {
            constexpr real S = M_SQRT1_2 + 1;
            // subtract corner and scale to recover a square of size sqrt(1/2)
            real cx = S * x - std::copysign(S, x);
            real cy = S * y - std::copysign(S, y);
            // apply rotation, scaling by sqrt(2): x' = y + x;  y' = y - x
            x = cy + cx;
            y = cy - cx;
        }
#endif
        real w = x * x + y * y;
        if (( w <= 1 ) & ( 0 < w ))
        {
            w = std::sqrt( std::log(w) / ( -0.5 * w ) );
            dst[0] = w * x;
            dst[1] = w * y;
            dst += 2;
        }
        src += 2;
    }
    return dst;
}


real * makeGaussians_std(real dst[], size_t cnt, const int32_t src[])
{
    std::random_device rd{};
    std::mt19937 gen{rd()};
    
    std::normal_distribution<> distribution{0,1};
    for ( size_t i = 0; i < cnt; ++i )
        dst[i] = distribution(gen);
    return dst + cnt;
}


real * makeExponentials_(real dst[], size_t cnt, const int32_t src[])
{
    for ( size_t i = 0; i < cnt; ++i )
    {
        real x = std::fabs(src[i]);
        dst[i] = -std::log(1 - x * TWO_POWER_MINUS_31);
    }
    return dst + cnt;
}


template < typename T >
void print_gaussian(size_t cnt, T const* vec)
{
    for ( size_t i = 0; i < cnt; )
    {
        for ( int k = 0; k < 4; ++k )
        {
            printf(" :");
            for ( int j = 0; j < 8 && i < cnt; ++j )
                printf(" %8.4f", vec[i++]);
        }
        printf("\n");
    }
}

template < typename REAL >
void check_gaussian(size_t cnt, REAL* vec)
{
    size_t nan = 0;
    double off = 1; // assumed mean
    double avg = 0, var = 0;
    for ( size_t i = 0; i < cnt; ++i )
    {
        if ( std::isnan(vec[i]) )
            ++nan;
        else
        {
            avg += vec[i];
            var += ( vec[i] - off ) * ( vec[i] - off );
        }
    }
    avg /= cnt;
    var = ( var - square( avg - off ) * cnt ) / ( cnt - 1 );
    // covariance of odd and even numbers:
    double cov = 0;
    for ( size_t i = 1; i < cnt; i += 2 )
    {
        if ( !std::isnan(vec[i]) )
            cov += ( vec[i-1] - avg ) * ( vec[i] - avg );
    }
    cnt -= nan;
    cov /= ( cnt / 2 );
    printf("%6lu + %6lu NaNs: avg %7.4f var %7.4f cov %7.4f ", cnt, nan, avg, var, cov);
}


//------------------------------------------------------------------------------
#pragma mark -
#include "random.h"

#include "simd.h"
#include "simd_math.h"
#include "../math/random_simd.cc"

#if defined(__AVX__)

#include "simd_float.h"
#include "simd_math.h"


/* Absolute error bounded by 1e-5 for normalized inputs
   Returns a finite number for +inf input
   Returns -inf for nan and <= 0 inputs.
   Continuous error.
 By Jacques-Henri Jourdan
 */
inline float logapprox(float val)
{
    union { float f; int32_t i; } valu;
    float exp, addcst, x;
    valu.f = val;
    exp = valu.i >> 23;
    /* -89.970756366f = -127 * log(2) + constant term of polynomial below. */
    /* -88.0296919311f = -127 * log(2) */
    /* -1.94106443489f = constant term */
    valu.i = (valu.i & 0x7FFFFF) | 0x3F800000;
    x = valu.f;
    
    /* Generated in Sollya using:
     > f = remez(log(x)-(x-1)*log(2),
     [|1,(x-1)*(x-2), (x-1)*(x-2)*x, (x-1)*(x-2)*x*x,
     (x-1)*(x-2)*x*x*x|], [1,2], 1, 1e-8);
     > plot(f+(x-1)*log(2)-log(x), [1,2]);
     > f+(x-1)*log(2)
     */
    /*
    float res = x * (3.529304993f + x * (-2.461222105f + x * (1.130626167f +
    x * (-0.288739945f + x * 3.110401639e-2f))))
    + ( -89.970756366f + 0.6931471805f*exp );
    */
    float cst = -89.970756366f + 0.6931471805f*exp;
    float res = 3.529304993f + x * (-2.461222105f + x * (1.130626167f +
                               x * (-0.288739945f + x * 3.110401639e-2f)));
    res = x * res + cst;
    return ( val > 0 ) ? res : -(float)INFINITY;
}


/// use this to check the log()
static real* check_log(real dst[], size_t cnt, const __m256i src[])
{
    const vec8f eps = set8f(0x1p-31f); //TWO_POWER_MINUS_31
    __m256i const* end = src + cnt;
    real * d = dst;
    while ( src < end )
    {
        // generate random floats in ]0, 1]:
        vec8f n = sub8f(set8f(1.0f), mul8f(eps, abs8f(cvt8if(load8i(src++)))));
        ++src; //vec8f j = mul8f(eps, cvt8if(load8i(src++)));
        vec8f x = n;
        vec8f y = logapprox8f(n);
#if REAL_IS_DOUBLE
        // convert 16 single-precision values
        store4d(d   , getlo4f(x));
        store4d(d+4 , gethi4f(x));
        store4d(d+8 , getlo4f(y));
        store4d(d+12, gethi4f(y));
#else
        // store 16 single-precision values
        store8f(d  , x);
        store8f(d+8, y);
#endif
        d += 16;
    }
    return dst+8*cnt;
}

#endif


template < float* (*FUNC)(float*, size_t, const int32_t*) >
void runGaussian(sfmt_t& sfmt, const char str[], int cnt)
{
    float vec[SFMT_N32] = { 0 };
    tick();
    for ( int i = 0; i < cnt; ++i )
    {
        sfmt_gen_rand_all(&sfmt);
        FUNC(vec, SFMT_N32, (int32_t*)sfmt.state);
    }
    float* end = FUNC(vec, SFMT_N32, (int32_t*)sfmt.state);
    printf("%-12s %5.2f :", str, tock(cnt>>10));
    check_gaussian(end-vec, vec);
    print_gaussian(std::min(end-vec, 16l), vec);
}

template < double* (*FUNC)(double*, size_t, const int32_t*) >
void runGaussian(sfmt_t& sfmt, const char str[], int cnt)
{
    double vec[SFMT_N32] = { 0 };
    tick();
    for ( int i = 0; i < cnt; ++i )
    {
        sfmt_gen_rand_all(&sfmt);
        FUNC(vec, SFMT_N32, (int32_t*)sfmt.state);
    }
    double* end = FUNC(vec, SFMT_N32, (int32_t*)sfmt.state);
    printf("%-12s %5.2f :", str, tock(cnt>>10));
    check_gaussian(end-vec, vec);
    print_gaussian(std::min(end-vec, 16l), vec);
}

/**
 Tests different implementation for speed
 */
void test_gaussian(int cnt)
{
    printf("test_gaussian --- %lu bytes real --- %s\n", sizeof(real), __VERSION__);
    sfmt_t sfmt;
    sfmt_init_gen_rand(&sfmt, time(nullptr));

    tick();
    for ( int i = 0; i < cnt; ++i )
        sfmt_gen_rand_all(&sfmt);
    printf("RNG.refill   %5.2f\n", tock(cnt>>10));
    //print(vec, end);
    
    runGaussian<makeGaussians_std>(sfmt, "GaussSTD", cnt);
    runGaussian<makeGaussians_>(sfmt, "Gauss", cnt);
#if USE_SIMD
    runGaussian<makeGaussians_SIMD>(sfmt, "GaussSSE", cnt);
#endif

#if defined(__AVX__)
    runGaussian<makeGaussians_AVXBM>(sfmt, "GaussAVXBM", cnt);
    runGaussian<makeGaussians_AVX1>(sfmt, "GaussAVX1", cnt);
    runGaussian<makeGaussians_AVX2>(sfmt, "GaussAVX2", cnt);
    runGaussian<makeExponentialsAVX>(sfmt, "Expon.AVX", cnt);
    if ( 0 )
    {
        printf("Approximate logarithm:\n");
        real vec[SFMT_N32] = { 0 };
        check_log(vec, SFMT_N256, (__m256i*)sfmt.state);
        print_gaussian(std::min(SFMT_N256, 32), vec);
    }
#endif
    runGaussian<makeExponentials_>(sfmt, "Exponential", cnt);
}

/**
 Prints many Gaussian distributecd random numbers
 */
void print_gaussian(int cnt)
{
    for ( int i = 0; i < cnt; ++i )
        printf("%10.5f\n", RNG.gauss());
}


int main(int argc, char* argv[])
{
    int mode = 0;
    RNG.seed();

    if ( argc > 1 )
        mode = atoi(argv[1]);

    if ( mode )
        print_gaussian(1<<14);
    else
        test_gaussian(1<<18);
}

