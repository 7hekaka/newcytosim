// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#include "real.h"
#include <cstdio>
#include <bitset>
#include <cstring>
#include <iostream>
#include <climits>
#include <random>
#include <new>

#include "timer.h"
#include "assert_macro.h"
#include "SFMT.h"

#include "simd.h"
#include "simd_float.h"
#include "simd_math.h"

std::random_device device;

template < typename T >
void print_floats(size_t cnt, T const* vec)
{
    for ( size_t i = 0; i < cnt; )
    {
        for ( int k = 0; k < 4; ++k )
        {
            if ( i >= cnt ) break;
            printf(" :");
            for ( int j = 0; j < 4; ++j )
                printf(" %10.3e", vec[i++]);
        }
    }
}


#define TWO_POWER_MINUS_31 0x1p-31
#define TWO_POWER_MINUS_32 0x1p-32
/**
 Fill array `vec[]` with Gaussian values ~ N(0,1).
 the size of `vec` should be a multiple of 2, and sufficient to hold `end-src` values
 @Return the number of values that were stored in `vec`
 */
real * makeGaussians_(real dst[], size_t cnt, const uint32_t arg[])
{
    int32_t const* src = reinterpret_cast<const int32_t*>(arg);
    int32_t const*const end = src + cnt;
    while ( src < end )
    {
        real x = src[0] * TWO_POWER_MINUS_31;
        real y = src[1] * TWO_POWER_MINUS_31;
#if 1
        if ( std::fabs(x) + std::fabs(y) >= M_SQRT2 )
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


real * makeGaussians_std(real dst[], size_t cnt, const uint32_t[])
{
    std::mt19937 gen{device()};
    
    std::normal_distribution<> distribution{0,1};
    for ( size_t i = 0; i < cnt; ++i )
        dst[i] = distribution(gen);
    return dst + cnt;
}

// defined in random.cc:
extern float * makeExponentials(float dst[], size_t cnt, const uint32_t src[]);

/** Using -log(R) where R is randomly distributed in [0, 1[ */
float * makeExponentials_(float dst[], size_t cnt, const uint32_t src[])
{
    const float alpha(TWO_POWER_MINUS_32);
    for ( size_t i = 0; i < cnt; ++i )
    {
        float x = static_cast<float>(src[i]);
        float y = alpha * x + alpha;
        assert_true( y > 0.f );
        dst[i] = -logf(y);
    }
    return dst + cnt;
}


template < typename REAL >
void check_gaussian(size_t num, REAL* vec, double off)
{
    size_t cnt = 0;
    size_t nan = 0;
    double avg = 0, var = 0;
    for ( size_t i = 0; i < num; ++i )
    {
        if ( std::isnan(vec[i]) )
            ++nan;
        else
        {
            ++cnt;
            double v = vec[i] - off;
            avg += v;
            var += v * v;
        }
    }
    avg = avg / cnt;
    var = ( var - square(avg) * cnt ) / ( cnt - 1 );
    avg += off;
    // covariance of odd and even numbers:
    double cov = 0;
    for ( size_t i = 1; i < num; i += 2 )
    {
        if ( vec[i] == vec[i] )
            cov += ( vec[i-1] - off ) * ( vec[i] - off );
    }
    cov /= ( cnt / 2 );
    printf("%9lu + %6lu NaNs: avg %7.4f var %7.4f cov %7.4f ", cnt, nan, avg, var, cov);
}


//------------------------------------------------------------------------------
#pragma mark -

#if defined(__AVX__)

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


/// use this to check the approximate log() function
static real* check_log(real dst[], size_t cnt, const uint32_t arg[])
{
    const vec8i * src = (vec8i*)arg;
    const vec8i * end = src + cnt / 8;

    const vec8f eps = set8f(0x1p-31f); //TWO_POWER_MINUS_31
    real * d = dst;
    while ( src < end )
    {
        // generate random floats in ]0, 1]:
        vec8f n = sub8f(set8f(1.0f), mul8f(eps, abs8f(load8if(src++))));
        ++src; //vec8f j = mul8f(eps, load8if(src++));
        vec8f x = n;
        vec8f y = log_approx8f(n);
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


template < typename REAL >
void run(sfmt_t& sfmt, size_t rep, const char str[], REAL* (*FUNC)(REAL*, size_t, const uint32_t*), double off = 0)
{
    void * ptr = nullptr;
    if ( 0 == posix_memalign(&ptr, 32, rep*sizeof(REAL)*SFMT_N32) )
    {
        REAL * mem = (REAL*)ptr;
        REAL * vec = mem;
        tick();
        for ( size_t i = 0; i < rep; ++i )
        {
            sfmt_gen_rand_all(&sfmt);
            vec = FUNC(vec, SFMT_N32, (uint32_t*)sfmt.state);
        }
        printf("\n%-12s %7.2f :", str, tock());
        check_gaussian(vec-mem, mem, off);
        
        size_t n = std::min(vec-mem, 8L);
        std::sort(mem, vec);
        print_floats(n, mem);
        print_floats(n, vec-n);
        
        size_t z = 0;
        for ( REAL * p = mem; p < vec; ++p )
            z += ( *p == 0 );
        printf(" (%lu zeros)", z);

        std::free(mem);
    }
}


template < typename REAL >
void scan(size_t chunk, size_t rep, REAL* (*FUNC)(REAL*, size_t, const uint32_t*))
{
    void * src = nullptr;
    void * dst = nullptr;

    int e = posix_memalign(&dst, 32, sizeof(REAL)*chunk);
    int f = posix_memalign(&src, 32, 4*chunk);
    if ( !e && !f )
    {
        uint32_t * rnd = (uint32_t*)src;
        
        for ( int r = 0; r < rep; ++r )
        {
            REAL * vec = (REAL*)dst;
            for ( uint32_t i = 0; i < chunk; ++i )
                rnd[i] = ~(r*chunk+i);
            FUNC(vec, chunk, rnd);
            for ( uint32_t i = 0; i < chunk; ++i )
                printf("%.4e ", vec[i]);
            printf("\n");
        }
        std::free(src);
        std::free(dst);
    }
}



int main(int argc, char* argv[])
{
    const size_t CNT = 1024*64;
    printf("test_gaussian --- %lu bytes real --- %s", sizeof(real), __VERSION__);
    sfmt_t sfmt alignas(32);
    sfmt_init_gen_rand(&sfmt, device());

    tick();
    for ( size_t i = 0; i < CNT; ++i )
        sfmt_gen_rand_all(&sfmt);
    printf("\nsfmt.refill  %7.2f", tock());
    //print(vec, end);
    
    run(sfmt, CNT, "Gauss.STD", makeGaussians_std);
    run(sfmt, CNT, "Gauss", makeGaussians_);

#if defined(__AVX__)
    run(sfmt, CNT, "GauBM.AVX", makeGaussiansBM_AVX);
    run(sfmt, CNT, "Gauss.AVX", makeGaussians_AVX);
#endif
    
    double one(1);
    run(sfmt, CNT, "Exponential", makeExponentials, one);
    run(sfmt, CNT, "Exponential_", makeExponentials_, one);
    //printf("\nSCAN makeExponentials_:\n");
    //scan(16, 64, makeExponentials_);
    printf("\n");
}

