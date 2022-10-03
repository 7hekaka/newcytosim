// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "random.h"
#include <cstdio>
#include <bitset>
#include <cstring>
#include <iostream>
#include <climits>
#include "timer.h"

template < typename T >
void print_bits(FILE* f, const T& val, char spc)
{
    unsigned char * ptr = (unsigned char*) & val;
    for ( int i = sizeof(T)-1; i >= 0; --i)
    {
        unsigned char byte = ptr[i];
        for ( int i = 0; i < CHAR_BIT; ++i )
        {
            putc('0' + (1 & (byte>>(CHAR_BIT-1))), f);
            byte <<= 1;
        }
        if ( spc ) putc(spc, f);
    }
    putc('\n', f);
}


void speed_test()
{
    const size_t cnt = 1 << 30;
    tick();
    uint32_t u = 10;
    for (size_t j=0; j<cnt; ++j)
    {
        u = RNG.pint32(1024);
        RNG.pint32(u);
    }
    printf("int %5.2f\n", tock(cnt));
}


void test_int()
{
    for (int j=0; j<8; ++j)
    {
        for (int k=0; k<8; ++k)
            printf(" %12u", RNG.pint32());
        printf("\n");
    }
    printf("\n");
    
    for (int j=0; j<8; ++j)
    {
        for (int k=0; k<8; ++k)
            printf(" %+12i", RNG.sint32());
        printf("\n");
    }
    printf("\n");

    for (int j=0; j<8; ++j)
    {
        for (int k=0; k<32; ++k)
            printf(" %2u", RNG.pint32(100));
        printf("\n");
    }
    printf("\n");
    
    for (int j=0; j<8; ++j)
    {
        for (int k=0; k<32; ++k)
            printf(" %2u", RNG.pint32_fair(100));
        printf("\n");
    }
    printf("\n");

    for (int j=0; j<8; ++j)
    {
        for (int k=0; k<32; ++k)
            printf(" %2u", RNG.pint32_slow(99));
        printf("\n");
    }
    printf("\n");
}


void silly_test()
{
    const uint32_t up = 1 << 30;
    
    const uint32_t cnt = 1 << 24;
    uint32_t hit = 0;
    
    for (uint32_t j=0; j<cnt; ++j)
        hit += ( RNG.pint32() < up );

    printf(" prob( pint32() < 1^30 ) = %f\n", hit/(float)cnt);
}


/**
 This assumes IEEE Standard 754 Floating point numbers
32 bits: 1 for sign, 8 for exponents, 23 for fraction
 */
float convertFix(uint32_t x)
{
    constexpr uint32_t FRAC  = 0x7FFFFFU;
    constexpr uint32_t EXPON = 127 << 23;
    uint32_t res = EXPON | ( x & FRAC );
    return *((float*)&res) - 1.0;
}


void testbits()
{
    const int SCALE = 2;
    for ( int i=0; i <= SCALE; ++i )
    {
        float x = i / float(SCALE);
        printf(" %f :", x);
        print_bits(stdout, x, 0);
        // x = -ii / float(SCALE);
        // printf("%f :", x);
        // print_bits(stdout, x, 0);
    }
    
    for ( int i=0; i < 16; ++i )
    {
        float y = convertFix(RNG.pint32());
        printf(" %f :", y);
        print_bits(stdout, y, ' ');
    }
}


#define TEST test
void test_test( const real prob, const size_t MAX )
{
    int cnt = 0, a, b, c;
    for ( size_t jj=0; jj < MAX; ++jj )
    {
        a = RNG.TEST(prob) + RNG.TEST(prob) + RNG.TEST(prob) + RNG.TEST(prob);
        b = RNG.TEST(prob) + RNG.TEST(prob) + RNG.TEST(prob) + RNG.TEST(prob);
        c = RNG.TEST(prob) + RNG.TEST(prob) + RNG.TEST(prob) + RNG.TEST(prob);
        cnt += a + b + c;
    }
    printf("prob = %f measured = %f cnt = %i\n", prob, cnt / double(12*MAX), cnt);
}

void test_RNG(const size_t MAX)
{
    for ( size_t jj=0; jj < MAX; ++jj )
    {
        RNG.preal();RNG.preal();RNG.preal();RNG.preal();RNG.preal();
        RNG.preal();RNG.preal();RNG.preal();RNG.preal();RNG.preal();
    }
}


void test_real()
{
    for (int j=0; j<8; ++j)
    {
        for (int k=0; k<8; ++k)
            printf(" %10f", RNG.sreal());
        printf("\n");
    }

    printf("\n");
    for (int j=0; j<8; ++j)
    {
        for (int k=0; k<8; ++k)
            printf(" %10f", RNG.preal());
        printf("\n");
    }
    
    printf("\n");
    for (int j=0; j<8; ++j)
    {
        for (int k=0; k<8; ++k)
            printf(" %10f", RNG.shalf());
        printf("\n");
    }

    printf("\npfloat:     ");
    float x;
    for ( int kk=0; kk < 10; ++kk )
    {
        x = RNG.pfloat();
        printf(" %+f", x);
    }

    printf("\nsfloat:     ");
    for ( int kk=0; kk < 10; ++kk )
    {
        x = RNG.sfloat();
        printf(" %+f", x);
    }
    
    double d;
    printf("\npdouble:    ");
    for ( int kk=0; kk < 10; ++kk )
    {
        d = RNG.pdouble();
        printf(" %+f", d);
    }

    printf("\nsdouble:    ");
    for ( int kk=0; kk < 10; ++kk )
    {
        d = RNG.sdouble();
        printf(" %+f", d);
    }

    printf("\nsflip:      ");
    for ( int kk=0; kk < 10; ++kk )
    {
        d = RNG.sflip();
        printf(" %+f", d);
    }
    printf("\n");
}

//==========================================================================

void test_uniform(size_t cnt)
{
    const double off = 0.5;
    double avg = 0, var = 0;
    for ( size_t i = 0; i < cnt; ++i )
    {
        real x = RNG.sreal() - off;
        real y = RNG.sreal() - off;
        real z = RNG.sreal() - off;
        real t = RNG.sreal() - off;
        avg += x + y + z + t;
        var += x*x + y*y + z*z + t*t;
    }
    cnt *= 4;
    avg /= cnt;
    var = ( var - square(avg) * cnt ) / real(cnt-1);
    printf("UNIFORM      avg = %.12e   var = %.12e\n", avg+off, var);
}


void test_gauss(size_t CNT)
{
    size_t cnt = 0;
    double avg = 0, var = 0;
    const size_t n_max = 1<<6;
    real vec[n_max] = { 0 };
    for ( size_t i = 0; i < CNT; ++i )
    {
        size_t n = RNG.pint32(n_max);
        RNG.gauss_set(vec, n);
        cnt += n;
        for ( size_t u = 0; u < n; ++u )
        {
            avg += vec[u];
            var += vec[u] * vec[u];
        }
    }
    avg /= cnt;
    var = ( var - square(avg) * cnt ) / real(cnt-1);
    printf("GAUSSIAN     avg = %.12e   var = %.12e\n", avg, var);
}


void test_prob()
{
    size_t avg = 0;
    size_t cnt = 1 << 28;
    for ( size_t i = 0; i < cnt; ++i )
        avg += RNG.flip_8th();

    printf("8th      prob = %.6f\n", avg/(double)cnt);
}


void test_exponential(size_t cnt)
{
    const double off = 1;
    double avg = 0, var = 0;
    for ( size_t i = 0; i < cnt; ++i )
    {
        real x = RNG.exponential() - off;
        real y = RNG.exponential() - off;
        real z = RNG.exponential() - off;
        real t = RNG.exponential() - off;
        avg += x + y + z + t;
        var += x*x + y*y + z*z + t*t;
    }
    cnt *= 4;
    avg /= cnt;
    var = ( var - square(avg) * cnt ) / real(cnt-1);
    printf("EXPONENTIAL  avg = %.12e   var = %.12e\n", avg+off, var);
}


void test_poisson(size_t sup)
{
    for ( size_t n = 0; n < sup; ++n )
    {
        int x = (int)(RNG.gauss() * std::sqrt(n) + n);
        printf("%10lu %9i %9i %9i\n", n, RNG.poisson_knuth(n), RNG.poisson(n), x);
    }
}


//==========================================================================


/**
 Fill array `vec[]` with Gaussian values ~ N(0,1).
 the size of `vec` should be a multiple of 2, and sufficient to hold `end-src` values
 @Return the number of values that were stored in `vec`
 */
template < typename REAL >
REAL * makeGaussians_(REAL dst[], size_t cnt, const int32_t src[])
{
    int32_t const*const end = src + cnt;
    while ( src < end )
    {
        REAL x = REAL(src[0]) * TWO_POWER_MINUS_31;
        REAL y = REAL(src[1]) * TWO_POWER_MINUS_31;
#if 1
        if ( std::abs(x) + std::abs(y) >= M_SQRT2 )
        {
            constexpr REAL S = M_SQRT1_2 + 1;
            // subtract corner and scale to recover a square of size sqrt(1/2)
            REAL cx = S * x - std::copysign(S, x);
            REAL cy = S * y - std::copysign(S, y);
            // apply rotation, scaling by sqrt(2): x' = y + x;  y' = y - x
            x = cy + cx;
            y = cy - cx;
        }
#endif
        REAL w = x * x + y * y;
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

template < typename REAL >
REAL * makeExponentials_(REAL dst[], size_t cnt, const int32_t src[])
{
    for ( size_t i = 0; i < cnt; ++i )
    {
        REAL x = std::fabs(static_cast<REAL>(src[i]));
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


#if defined(__AVX__)

#include "simd.h"
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

/*
 Absolute error bounded by 1e-5 for normalized inputs
   Returns a finite number for +inf input
   Returns -inf for nan and <= 0 inputs.
   Continuous error.
 By Jacques-Henri Jourdan, SIMD by Francois Nedelec
 */
inline vec8f log_approx8f(__m256 xxx)
{
    // masks:
    const __m256 mant = _mm256_castsi256_ps(_mm256_set1_epi32(0x007FFFFF));
    const __m256 expo = _mm256_castsi256_ps(_mm256_set1_epi32(0x3F800000));
    // polynomial coefficients
    const __m256 a = _mm256_set1_ps(+3.529304993f);
    const __m256 b = _mm256_set1_ps(-2.461222105f);
    const __m256 c = _mm256_set1_ps(+1.130626167f);
    const __m256 d = _mm256_set1_ps(-0.288739945f);
    const __m256 e = _mm256_set1_ps(+3.110401639e-2f);
    const __m256 f = _mm256_set1_ps(-89.970756366f);
    const __m256 g = _mm256_set1_ps(0.6931471805f);
    // used to clear negative / NaN arguments:
    __m256 invalid = _mm256_cmp_ps(xxx, _mm256_setzero_ps(), _CMP_NGT_UQ);
    // extract exponent:
#if defined(__AVX2__)
    __m256 cst = _mm256_cvtepi32_ps(_mm256_srli_epi32(_mm256_castps_si256(xxx), 23));
#else
    vec4f h = cvt4if(_mm_srli_epi32(_mm_castps_si128(gethi4f(xxx)), 23));
    vec4f l = cvt4if(_mm_srli_epi32(_mm_castps_si128(getlo4f(xxx)), 23));
    vec8f cst = cat44f(l, h);
#endif
    cst = _mm256_add_ps(_mm256_mul_ps(cst, g), f);
    // clear exponents:
    xxx = _mm256_or_ps(expo, _mm256_and_ps(mant, xxx));
    // calculate polynom:
#if 0 && defined(__FMA__)
    __m256 tmp = _mm256_fmadd_ps(xxx, e, d);
    tmp = _mm256_fmadd_ps(xxx, tmp, c);
    tmp = _mm256_fmadd_ps(xxx, tmp, b);
    tmp = _mm256_fmadd_ps(xxx, tmp, a);
    tmp = _mm256_fmadd_ps(xxx, tmp, cst);
#else
    __m256 tmp = _mm256_add_ps(_mm256_mul_ps(xxx, e), d);
    tmp = _mm256_add_ps(_mm256_mul_ps(xxx, tmp), c);
    tmp = _mm256_add_ps(_mm256_mul_ps(xxx, tmp), b);
    tmp = _mm256_add_ps(_mm256_mul_ps(xxx, tmp), a);
    tmp = _mm256_add_ps(_mm256_mul_ps(xxx, tmp), cst);
#endif
    // clear negative arguments:
    return _mm256_or_ps(tmp, invalid);
}


#if ( 0 )
/**
 Another implementation of logapproxf(float) found online...
 */
__m256 logf_app(__m256 val)
{
    // constants
    const __m256 mant = _mm256_castsi256_ps(_mm256_set1_epi32(0x007FFFFF)); // mantissa mask
    const __m256 expo = _mm256_castsi256_ps(_mm256_set1_epi32(0x3F800000)); // exponent mask
    // polynomial coefficients
    const __m256 a = _mm256_set1_ps(+3.529304993f);
    const __m256 b = _mm256_set1_ps(-2.461222105f);
    const __m256 c = _mm256_set1_ps(+1.130626167f);
    const __m256 d = _mm256_set1_ps(-0.288739945f);
    const __m256 e = _mm256_set1_ps(+3.110401639e-2f);
    const __m256 f = _mm256_set1_ps(-89.970756366f);
    const __m256 log2 = _mm256_set1_ps(M_LN2);      // log(2)

    // mask out anything <= 0
    __m256 invalid = _mm256_cmp_ps(val, _mm256_setzero_ps(), _CMP_NGT_UQ);

    // extract exponents
#if defined(__AVX2__)
    __m256 exp = _mm256_cvtepi32_ps(_mm256_srli_epi32((__m256i)val, 23));
#else
    __m128 hi = _mm256_extractf128_ps(val, 1);
    __m128 lo = _mm256_extractf128_ps(val, 0);
    hi = _mm_cvtepi32_ps(_mm_srli_epi32(_mm_castps_si128(hi), 23));
    lo = _mm_cvtepi32_ps(_mm_srli_epi32(_mm_castps_si128(lo), 23));
    __m256 exp = _mm256_set_m128(hi, lo);
#endif

    // clear exponent to 0 (+127)
    val = _mm256_and_ps(val, mant);
    val = _mm256_or_ps (val, expo);

    // horner's rule to evaluate polynomial
    __m256 res = e;
    res = _mm256_add_ps(_mm256_mul_ps(res,val), d);
    res = _mm256_add_ps(_mm256_mul_ps(res,val), c);
    res = _mm256_add_ps(_mm256_mul_ps(res,val), b);
    res = _mm256_add_ps(_mm256_mul_ps(res,val), a);
    res = _mm256_add_ps(_mm256_mul_ps(res,val), f);
    // add in exponent contribution
    res = _mm256_add_ps(_mm256_mul_ps(exp, log2), res);
    return _mm256_or_ps(res, invalid);
}
#endif


/// use this to check the log()
static real* check_log(real dst[], size_t cnt, const __m256i src[])
{
    const vec8f eps = set8f(0x1p-31f); //TWO_POWER_MINUS_31
    __m256i const* end = src + cnt;
    real * d = dst;
    while ( src < end )
    {
        // generate random floats in ]0, 1]:
        vec8f n = sub8f(set8f(1.0f), mul8f(eps, abs8f(cvt8if(load8si(src++)))));
        ++src; //vec8f j = mul8f(eps, cvt8if(load8si(src++)));
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

#include "random_simd.cc"

#endif


template < float* (*FUNC)(float*, size_t, const int32_t*) >
void runGaussian(sfmt_t& sfmt, const char str[], int cnt)
{
    float flt[SFMT_N32] = { 0 };
    tick();
    for ( int i = 0; i < cnt; ++i )
    {
        sfmt_gen_rand_all(&sfmt);
        FUNC(flt, SFMT_N32, (int32_t*)sfmt.state);
    }
    float* end = FUNC(flt, SFMT_N32, (int32_t*)sfmt.state);
    printf("%-12s %5.2f :", str, tock(cnt));
    check_gaussian(end-flt, flt);
    print_gaussian(std::min(end-flt, 16l), flt);
}


#if defined(__AVX__)
template < real* (*FUNC)(real*, size_t, const __m256i*) >
void runGaussian(sfmt_t& sfmt, const char str[], int cnt)
{
    real *end, vec[SFMT_N32];
    for ( int i = 0; i < SFMT_N32; ++i )
        vec[i] = 0;
    tick();
    for ( int i = 0; i < cnt; ++i )
    {
        sfmt_gen_rand_all(&sfmt);
        end = FUNC(vec, SFMT_N256, (__m256i*)sfmt.state);
    }
    printf("%-12s %5.2f :", str, tock(cnt));
    check_gaussian(end-vec, vec);
    print_gaussian(std::min(end-vec, 16l), vec);
}
#endif

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
    printf("RNG.refill   %5.2f\n", tock(cnt));
    //print(vec, end);
    
    runGaussian<makeGaussians_>(sfmt, "Gauss_", cnt);
#if defined(__AVX__)
    real *end, vec[SFMT_N32] = { 0 };
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


//==========================================================================
int main(int argc, char* argv[])
{
    int mode = 4;
    RNG.seed();

    if ( argc > 1 )
        mode = atoi(argv[1]);
    real rate = 1;
    if ( argc > 2 )
        rate = strtod(argv[2], 0);

    switch ( mode )
    {
        case 0:
            test_poisson(1024);
            test_prob();
            break;
            
        case 1:
            test_exponential(0x1p26);
            test_uniform(0x1p20);
            test_gauss(0x1p20);
            break;

        case 3:
            for ( int kk=0; kk < 11; ++kk )
                test_test(rate*kk, 5000000);
            break;
            
        case 4:
            printf("sizeof(uint32_t) = %lu\n", sizeof(uint32_t));
            test_int();
            test_real();
            break;
            
        case 5:
            speed_test();
            break;
            
        case 6:
            silly_test();
            break;
            
        case 7:
            test_gaussian(1<<18);
            break;
            
        case 8:
            print_gaussian(1<<14);
            break;
    }
}

