// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "random.h"
#include <cstdio>
#include <bitset>
#include <cstring>
#include <iostream>
#include <climits>
#include "tictoc.h"

using namespace TicToc;

template < typename T >
void print_bits(FILE* f, const T& val, char spc = 0)
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
    tic();
    uint32_t u = 10;
    for (size_t j=0; j<cnt; ++j)
    {
        u = RNG.pint32(1024);
        RNG.pint32(u);
    }
    printf("int %5.2f\n", toc(cnt));
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
        print_bits(stdout, x);
        // x = -ii / float(SCALE);
        // printf("%f :", x);
        // print_bits(stdout, x);
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

void test_uniform()
{
    size_t cnt = 1<<28;
    real avg = 0;
    real var = 0;
    for ( size_t i = 0; i < cnt; ++i )
    {
        real x = RNG.sreal();
        real y = RNG.sreal();
        real z = RNG.sreal();
        real t = RNG.sreal();
        avg += x + y + z + t;
        var += x*x + y*y + z*z + t*t;
    }
    cnt *= 4;
    avg /= (real)cnt;
    var = var/(real)cnt - avg * avg;
    printf("UNIFORM      avg = %.12e   var = %.12e\n", avg, var);
}


void test_gauss()
{
    printf("Gauss\n");
    size_t cnt = 0;
    real avg = 0;
    real var = 0;
    const size_t n_max = 1<<6;
    real vec[n_max] = { 0 };
    for ( size_t i = 0; i < 10000000; ++i )
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
    avg /= (real)cnt;
    var = var/(real)cnt - avg * avg;
    printf("GAUSS      avg = %.12e   var = %.12e\n", avg, var);

}


void test_prob()
{
    size_t avg = 0;
    size_t cnt = 1 << 28;
    for ( size_t i = 0; i < cnt; ++i )
        avg += RNG.flip_8th();

    printf("8th      prob = %.6f\n", avg/(double)cnt);
}


void test_exponential()
{
    size_t cnt = 1 << 29;
    real avg = 0;
    real var = 0;
    for ( size_t i = 0; i < cnt; ++i )
    {
        real x = RNG.exponential();
        real y = RNG.exponential();
        real z = RNG.exponential();
        real t = RNG.exponential();
        avg += x + y + z + t;
        var += x*x + y*y + z*z + t*t;
    }
    cnt *= 4;
    avg /= (real)cnt;
    var = var/(real)cnt - avg * avg;
    printf("EXPONENTIAL  avg = %.12e   var = %.12e\n", avg, var);
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
//test 3 methods to generate a random event time, when the rate varies in time
// F. Nedelec, Oct 2005

//this is our standard method: 64s CPU
int method1(const int maxTime, const real rate[])
{
    for ( int ii=0; ii<maxTime; ++ii )
    {
        if (RNG.test(rate[ii])) return ii;
    }
    return maxTime;
}

//this is 'exact' and very slow: 370s CPU (an exponential at each step!)
int method2(const int maxTime, const real rate[])
{
    for ( int ii=0; ii<maxTime; ++ii )
    {
        if ( RNG.preal() < -std::expm1(-rate[ii]) )
            return ii;
    }
    return maxTime;
}

//this is exact, and the fastest method: 10s CPU!
int method3(const int maxTime, const real rate[])
{
    real T = -std::log( RNG.preal() );
    for ( int ii=0; ii<maxTime; ++ii )
    {
        T -= rate[ii];
        if ( T < 0 ) return ii;
    }
    return maxTime;
}


int testGillespie(const int method)
{
    //test new idea for gillespie with changing rate (Oct 2005)
    const int maxTime = 200;
    real rate[maxTime];
    for ( int ii=0; ii<maxTime; ++ii )
        rate[ii] = ( ii % 10 ) / 30.0;
    
    int bins[3][maxTime+1];
    for ( int ii=0; ii<=maxTime; ++ii )
    {
        bins[0][ii] = 0;
        bins[1][ii] = 0;
        bins[2][ii] = 0;
    }
    
    const int nbSamples = 1000000;
    const int subSamples = 10;
    int result;
    switch( method )
    {
        case 0:
            for ( int ii=0; ii<nbSamples; ++ii )
            {
                bins[0][ method1(maxTime, rate) ]++;
                bins[1][ method2(maxTime, rate) ]++;
                bins[2][ method3(maxTime, rate) ]++;
            }
            break;
            
        case 1:
            printf("method 1:");
            for ( int ii=0; ii<nbSamples; ++ii )
                for ( int jj=0; jj<subSamples; ++jj )
                    result = method1(maxTime, rate);
            return result;
            
        case 2:
            printf("method 2:");
            for ( int ii=0; ii<nbSamples; ++ii )
                for ( int jj=0; jj<subSamples; ++jj )
                    result = method2(maxTime, rate);
            return result;
            
        case 3:
            printf("method 3:");
            for ( int ii=0; ii<nbSamples; ++ii )
                for ( int jj=0; jj<subSamples; ++jj )
                    result = method3(maxTime, rate);
            return result;
    }
    
    FILE* file = fopen("test.out", "w");
    for ( int ii=0; ii<=maxTime; ++ii )
        fprintf(file, "%4i   %6i %6i %6i\n", ii, bins[0][ii], bins[1][ii], bins[2][ii]);
    fclose(file);
    return 0;
}


//==========================================================================


/**
 Fill array `vec[]` with Gaussian values ~ N(0,1).
 the size of `vec` should be a multiple of 2, and sufficient to hold `end-src` values
 @Return the number of values that were stored in `vec`
 */
real * gauss_fill_0(real dst[], const int32_t src[], int32_t const*const end)
{
    while ( src < end )
    {
        real x = src[0] * TWO_POWER_MINUS_31;
        real y = src[1] * TWO_POWER_MINUS_31;
        real w = x * x + y * y;
        if (( w <= 1 ) & ( 0 < w ))
        {
            w = std::sqrt( std::log(w) / ( -0.5 * w ) );
            *dst++ = w * x;
            *dst++ = w * y;
        }
        src += 2;
    }
    return dst;
}

#if defined(__AVX__)

#include "simd.h"
#include "simd_float.h"

/**
 This packs an array by removing 'nan' values. The order of the list is irrelevant.
 This implementation is quite poor and takes much more time than the calculation
 of the Gaussian to start with, so there must be a better way to proceed!
 */
template < typename T >
T * remove_nans(T * s, T * e)
{
    while ( s < e )
    {
        --e;
        // find the next `nan` going upward:
        while ( *s == *s )
        {
            if ( ++s > e )
                return s;
        }
        // skip `nan` values going downward:
        while ( *e != *e )
        {
            if ( --e <= s )
                return s;
        }
        // copy number over:
        *s++ = *e;
    }
    return s;
}


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

#if 1
// convenience functions for defining constant vectors
inline static __m256 MM256_INT32(int32_t arg) { return _mm256_castsi256_ps(_mm256_set1_epi32(arg)); }
inline static __m256 MM256_FLOAT(float arg) { return _mm256_set1_ps(arg); }
#endif

/*
 Absolute error bounded by 1e-5 for normalized inputs
   Returns a finite number for +inf input
   Returns -inf for nan and <= 0 inputs.
   Continuous error.
 By Jacques-Henri Jourdan, SIMD by Francois Nedelec
 */
inline static vec8f logapprox8f(__m256 xxx)
{
    // masks:
    const __m256 mant = MM256_INT32(0x007fffff);
    const __m256 expo = MM256_INT32(0x3f800000);
    // polynomial coefficients
    const __m256 a = MM256_FLOAT(+3.529304993f);
    const __m256 b = MM256_FLOAT(-2.461222105f);
    const __m256 c = MM256_FLOAT(+1.130626167f);
    const __m256 d = MM256_FLOAT(-0.288739945f);
    const __m256 e = MM256_FLOAT(+3.110401639e-2f);
    const __m256 f = MM256_FLOAT(-89.970756366f);
    const __m256 g = MM256_FLOAT(0.6931471805f);
    // used to clear negative / NaN arguments:
    __m256 invalid = _mm256_cmp_ps(xxx, _mm256_setzero_ps(), _CMP_NGT_US);
    // extract exponent:
    __m256 cst = _mm256_cvtepi32_ps(_mm256_srli_epi32(_mm256_castps_si256(xxx), 23));
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


/**
 Another implementation of logapproxf(float) found online...
 */
__m256 logf_app(__m256 val)
{
    // constants
    const __m256 mant = MM256_INT32(0x007FFFFF); // mantissa mask
    const __m256 expo = MM256_INT32(0x3F800000); // exponent mask
    // polynomial coefficients
    const __m256 a = MM256_FLOAT(+3.529304993f);
    const __m256 b = MM256_FLOAT(-2.461222105f);
    const __m256 c = MM256_FLOAT(+1.130626167f);
    const __m256 d = MM256_FLOAT(-0.288739945f);
    const __m256 e = MM256_FLOAT(+3.110401639e-2f);
    const __m256 f = MM256_FLOAT(-89.970756366f);
    const __m256 log2 = MM256_FLOAT(M_LN2);      // log(2)

    // mask out anything <= 0
    __m256 invalid = _mm256_cmp_ps(val, _mm256_setzero_ps(), _CMP_NGT_UQ);

    // extract exponents
    __m256 exp = _mm256_cvtepi32_ps(_mm256_srli_epi32((__m256i)val, 23));

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

/*
 // generate random floats in [1, 2]:
 const vec8f mant = MM256_INT32(0x807fffff);
 const vec8f expo = MM256_INT32(0x3f800000);
 const vec8f ninf = MM256_INT32(0xff800000); // -INFINITY
 x = or8f(expo, and8f(mant, cast8f(load8si(src++))));
 y = or8f(expo, and8f(mant, cast8f(load8si(src++))));
*/

/**
 Calculates Gaussian-distributed, single precision random number,
 using SIMD AVX instructions
 Array `dst` should be able to hold as many 32-bit numbers as `src`.
 if 'real==float', for 256 bits of input, this produces ~64*PI bits of numbers.
 if 'real==double', this produces more output bits than input!

 The function used to calculate logarithm on SIMD data is part of the
 Intel SVML library, and is provided by the Intel compiler.

 F. Nedelec 02.01.2017
 */
real * gauss_fill_AVX(real dst[], const __m256i src[], __m256i* end)
{
    const vec8f eps = set8f(TWO_POWER_MINUS_31);
    const vec8f half = set8f(-0.5);
    
    real * d = dst;
    while ( src < end )
    {
        // generate random floats in [-1, 1]:
        vec8f x = mul8f(eps, cvt8i(load8si(src++)));
        vec8f y = mul8f(eps, cvt8i(load8si(src++)));
        vec8f n = add8f(mul8f(x,x), mul8f(y,y));
        //w = std::sqrt( -2 * std::log(n) / n );
        n = sqrt8f(div8f(logapprox8f(n), mul8f(half, n)));
        //n = rsqrt8f(div8f(mul8f(half, n), log8f(n)));
        //n = rsqrt8f(mul8f(mul8f(half, n), rcp8f(log8f(n))));
        x = mul8f(n, x);
        y = mul8f(n, y);
        // place corresponding X and Y values next to each other:
        vec8f t = unpacklo8f(x, y);
        y = unpackhi8f(x, y);
        // swap the NaNs to move them to 'y':
        x = blend8f(t, y, isnan8f(t));
        y = blend8f(y, t, isnan8f(t));
#if 0   // can push more NaNs out by swapping within the lanes:
        // swap the NaNs to move them to 'y':
        t = permute44f(x);
        x = blend8f(t, y, isnan8f(t)); d
        y = blend8f(y, t, isnan8f(t));
        // swap the NaNs to move them to 'y':
        t = swap2f128(x);
        x = blend8f(t, y, isnan8f(t));
        y = blend8f(y, t, isnan8f(t));
        // swap the NaNs to move them to 'y':
        t = permute44f(x);
        x = blend8f(t, y, isnan8f(t));
        y = blend8f(y, t, isnan8f(t));
#endif
#if REAL_IS_DOUBLE
        // convert 16 single-precision values
        store4(d   , cvt4sd(getlo4f(x)));
        store4(d+4 , cvt4sd(getlo4f(y)));
        store4(d+8 , cvt4sd(gethi4f(x)));
        store4(d+12, cvt4sd(gethi4f(y)));
#else
        // store 16 single-precision values
        store8f(d  , x);
        store8f(d+8, y);
#endif
        d += 16;
    }
    return d;
    //return remove_nans(dst, d);
}


/// use this to check the log()
real * check_log(real dst[], const __m256i src[], __m256i* end)
{
    const vec8f eps = set8f(TWO_POWER_MINUS_31);
    real * d = dst;
    while ( src < end )
    {
        // generate random floats in [-1, 1]:
        vec8f x = mul8f(eps, cvt8i(load8si(src++)));
        vec8f y = mul8f(eps, cvt8i(load8si(src++)));
        vec8f n = add8f(mul8f(x,x), mul8f(y,y));
        x = log8f(y);
        y = logapprox8f(y);
#if REAL_IS_DOUBLE
        // convert 16 single-precision values
        store4(d   , cvt4sd(getlo4f(x)));
        store4(d+4 , cvt4sd(getlo4f(y)));
        store4(d+8 , cvt4sd(gethi4f(x)));
        store4(d+12, cvt4sd(gethi4f(y)));
#else
        // store 16 single-precision values
        store8f(d  , x);
        store8f(d+8, y);
#endif
        d += 16;
    }
    return d;
    //return remove_nans(dst, d);
}


#endif


template < typename T >
void print_gaussian(size_t cnt, T const* vec)
{
    for ( size_t i = 0; i < cnt; )
    {
        for ( int j = 0; j < 8 && i < cnt; ++j, ++i )
            printf(" %8.4f", vec[i]);
        printf(" :");
        for ( int j = 0; j < 8 && i < cnt; ++j, ++i )
            printf(" %8.4f", vec[i]);
        printf("\n");
    }
}

void check_gaussian(size_t cnt, real* vec)
{
    size_t nan = 0;
    real avg = 0, dev = 0;
    for ( size_t i = 0; i < cnt; ++i )
    {
        if ( vec[i] != vec[i] )
            ++nan;
        else
        {
            avg += vec[i];
            dev += vec[i] * vec[i];
        }
    }
    cnt -= nan;
    avg /= cnt;
    dev = dev / cnt - avg;
    printf("%6lu numbers + %lu NaNs: avg %7.4f dev %7.4f\n", cnt, nan, avg, dev);
}

void test_gaussian(int cnt)
{
    printf("test_gaussian --- %lu bytes real --- %s\n", sizeof(real), __VERSION__);
    int32_t * buf = (int32_t*)RNG.data();

    if ( 1 )
    {
        tic();
        for ( int i = 0; i < cnt; ++i )
            RNG.refill();
        printf("RNG.refill   %5.2f\n", toc(cnt));
        //print(vec, end);
    }
    if ( 1 )
    {
        real *end, vec[SFMT_N32] = { 0 };
        tic();
        for ( int i = 0; i < cnt; ++i )
        {
            end = gauss_fill_0(vec, buf, buf+SFMT_N32);
            RNG.refill();
        }
        printf("gauss0       %5.2f  :", toc(cnt));
        check_gaussian(end-vec, vec);
        //print_gaussian(end-vec, vec);
    }
    __m256i * mem = (__m256i*)buf;
    if ( 1 )
    {
        real *end, vec[SFMT_N32] = { 0 };
        tic();
        for ( int i = 0; i < cnt; ++i )
        {
            end = gauss_fill_AVX(vec, mem, mem+SFMT_N256);
            RNG.refill();
        }
        printf("gauss AVX    %5.2f  :", toc(cnt));
        check_gaussian(end-vec, vec);
        print_gaussian((end-vec)/8, vec);
    }
    if ( 1 )
    {
        printf("Approximate logarithm:\n");
        real *end, vec[SFMT_N32] = { 0 };
        end = check_log(vec, mem, mem+SFMT_N256);
        print_gaussian(std::min(end-vec, 32L), vec);
    }
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
            test_exponential();
            test_uniform();
            test_gauss();
            break;
    
        case 2:
            testGillespie(rate);
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
    }
    
    printf("done\n");
}

