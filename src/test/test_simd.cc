// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/*
 Testing Intel's Streaming SIMD
 FJN, started July 2013
 
 To compile: c++ -O2 -mavx test.cc
 To generate assembly: c++ -S test.cc
 */

#include <cstdio>
#include "timer.h"
#include <cstdint>
#include "vecprint.h"
#include "simd.h"
#include "simd_float.h"
#include "simd_print.h"

#define shuffle2(a,b,k)   _mm_shuffle_pd(a,b,k)

template <typename FLOAT>
void dump(size_t len, const FLOAT* vec)
{
    printf(": ");
    for ( size_t n = 0; n < len; ++n )
    {
        printf("%+5.2f ", vec[n]);
        if (( n & 3 ) == 3 )
            printf(" ");
    }
    printf("\n");
}

template <typename FLOAT>
void dump(size_t lin, size_t col, const FLOAT* vec)
{
    for ( size_t i = 0; i < lin; ++i )
    {
        printf("%lu: ", i);
        for ( size_t j = 0; j < col; ++j )
            printf("%+5.2f ", vec[i+j*col]);
        printf("\n");
    }
}


void test_swapSSE()
{
    vec2 a = setr2(1, 2);
    vec2 b = setr2(3, 4);
    dump(a, "a");
    dump(b, "b");
    
#ifdef __AVX__
    dump(_mm_permute_pd(b,0b00), "permute 0b00");
    dump(_mm_permute_pd(b,0b01), "permute 0b01");
    dump(_mm_permute_pd(b,0b10), "permute 0b10");
    dump(_mm_permute_pd(b,0b11), "permute 0b11");
#endif
#if defined(__SSE3__)
    dump(shuffle2(a,b,0b00), "0b00");
    dump(shuffle2(a,b,0b01), "0b01");
    dump(shuffle2(a,b,0b10), "0b10");
    dump(shuffle2(a,b,0b11), "0b11");
#endif
    dump(unpacklo2(a,b), "unpacklo");
    dump(unpackhi2(a,b), "unpackhi");
    dump(duplo2(a), "duplo2(a)");
    dump(duphi2(a), "duphi2(a)");
    dump(blend11(a,b), "blend11(a,b)");
    dump(catshift(a,b), "catshift(a,b)");
    dump(swap2(a), "swap2(a)");
}


void test_shuffle()
{
    printf("---------- shuffle\n");
    // two shuffles to perform a blend
    vec4f a { 1, 2, 3, 4 };
    vec4f b {-1,-2,-3,-4 };
    dump(a, "a  ");
    dump(b, "b  ");

#if defined(__SSE3__)
    dump(_mm_movehl_ps(a, b), "_mm_movehl_ps(a, b)");
    dump(_mm_movelh_ps(a, b), "_mm_movelh_ps(a, b)");
    dump(_mm_movehl_ps(b, a), "_mm_movehl_ps(b, a)");
    dump(_mm_movelh_ps(b, a), "_mm_movelh_ps(b, a)");

    dump(_mm_shuffle_ps(a, b, 0x4E), "_mm_shuffle_ps(a, b, 0x4E)");
    dump(_mm_shuffle_ps(a, b, 0xEE), "_mm_shuffle_ps(a, b, 0xEE)");
    dump(_mm_shuffle_ps(a, b, 0xE4), "_mm_shuffle_ps(a, b, 0xE4)");
    dump(_mm_shuffle_ps(a, b, 0xC4), "_mm_shuffle_ps(a, b, 0xC4)");

    dump(_mm_shuffle_ps(b, a, 0x24), "_mm_shuffle_ps(b, a, 0x24)");
    dump(_mm_shuffle_ps(a, b, 0x23), "_mm_shuffle_ps(a, b, 0x23)");
    
    vec4f x = _mm_shuffle_ps(a, b, 0xEE);
    dump(_mm_shuffle_ps(a, x, 0xC4), "blend13");
    dump(_mm_shuffle_ps(a, b, 0xE4), "blend22");
    vec4f y = _mm_shuffle_ps(a, b, 0x44);
    dump(_mm_shuffle_ps(y, b, 0xEC), "blend31");
#endif
}


void test_signselect()
{
    vec4f s { -1, +1, -1, +1 };
    vec4f n { -2, -2, -2, -2 };
    vec4f p { +3, +3, +3, +3 };
    vec4f x = signselect4f(s, n, p);
    dump(x, "signselect");
}


void test_float4()
{
    printf("---------- float4\n");
    vec4f x{ 1, 2, 3, 4 };
    vec4f y{-5,-6,-7,-8 };
    dump(x, "x");
    dump(y, "y");
    
    dump(catshift1f(x, y), "catshift1f(x, y)");
    dump(catshift2f(x, y), "catshift2f(x, y)");
    dump(catshift3f(x, y), "catshift3f(x, y)");

    dump(blend13f(x, y), "blend13f(x, y)");
    dump(blend22f(x, y), "blend22f(x, y)");
    dump(blend31f(x, y), "blend31f(x, y)");
}

//------------------------------------------------------------------------------
#pragma mark -

#ifdef __AVX__

/**
 make dst = { XYZ XYZ XYZ XYZ }
 from src = { XYZ? }
 */
inline void twinedup12(double const* src, double* dst)
{
    vec4 s = load3(src);
    vec4 p = swap2f128(s);
    vec4 h = shuffle4(s, p, 0b0001);
    vec4 d0 = blend31(s, h);
    vec4 d1 = blend22(h, p);
    vec4 d2 = shuffle4(p, s, 0b0100);
    store4(dst  , d0);
    store4(dst+4, d1);
    store4(dst+8, d2);
}


/**
 make
     dst = { XXX YYY ZZZ TTT }
 from
     src = { XYZT }
 */
inline void repeat12(double const* src, double* dst)
{
    vec4 s = load4(src);
    dump(s, "s");

    vec4 zx = swap2f128(s);
    vec4 xy = unpacklo4(s, s);
    vec4 yz = unpackhi4(s, s);
    
    store4(dst  , blend22(xy, zx));
    store4(dst+4, blend22(yz, xy));
    store4(dst+8, blend22(zx, yz));
}


/**
 make
     dst = { XYZ XYZ XYZ XYZ }
 from
     dX = { XXXX }
     dY = { YYYY }
     dZ = { ZZZZ }
 */
inline void twine12(double const* X, double const* Y, double const* Z, double* dst)
{
    vec4 sx = load4(X);
    vec4 sy = load4(Y);
    vec4 sz = load4(Z);
    dump(sx, "sx");
    dump(sy, "sy");
    dump(sz, "sz");

    vec4 zx = blend0101(sz, sx);
    zx = swap2f128(zx);
    vec4 xy = unpacklo4(sx, sy);
    vec4 yz = unpackhi4(sy, sz);
    
    store4(dst  , blend22(xy, zx));
    store4(dst+4, blend22(yz, xy));
    store4(dst+8, blend22(zx, yz));
}

/**
 make
     dX = { XXXX }
     dY = { YYYY }
     dZ = { ZZZZ }
 from src = { XYZ XYZ XYZ XYZ }
 */
inline void untwine12(double const* src, double* X, double* Y, double* Z)
{
    vec4 s0 = load4(src);
    vec4 s1 = load4(src+4);
    vec4 s2 = load4(src+8);
    dump(s0, "s0");
    dump(s1, "s1");
    dump(s2, "s2");

    vec4 zx = blend22(s2, s0);
    zx = swap2f128(zx);
    vec4 xy = blend22(s0, s1);
    vec4 yz = blend22(s1, s2);
    
    store4(X, blend0101(xy, zx));
    store4(Y, shuffle4(xy, yz, 0b0101));
    store4(Z, blend0101(zx, yz));
}


/**
 make
     dst = { X+X+X, Y+Y+Y, Z+Z+Z, ? }
 from
     src = { XYZ XYZ XYZ XYZ }
 */
inline void sumXXX(double const* src, double* dst)
{
    vec4 s0 = load4(src);
    vec4 s1 = load4(src+4);
    vec4 s2 = load4(src+8);
    
    vec4 h = shuffle4(blend31(s1, s0), s2, 0b0101);
    vec4 d3 = catshift2(s1, s2);
    vec4 d2 = shuffle4(s2, s1, 0b0101);
    vec4 d1 = swap2f128(h);

    vec4 sum = add4(add4(s0, d2), add4(d3, d1));
    store4(dst, sum);
}


/**
 make
     dst = { X+Y+Z, X+Y+Z, X+Y+Z, X+Y+Z }
 from
     src = { XYZ XYZ XYZ XYZ }
 */
inline void sumXYZ(double const* src, double* dst)
{
    vec4 s0 = load4(src);
    vec4 s1 = load4(src+4);
    vec4 s2 = load4(src+8);
    
    vec4 zx = blend22(s2, s0);
    zx = swap2f128(zx);
    vec4 xy = blend22(s0, s1);
    vec4 yz = blend22(s1, s2);
    
    vec4 d1 = blend0101(xy, zx);
    vec4 d2 = shuffle4(xy, yz, 0b0101);
    vec4 d3 = blend0101(zx, yz);

    vec4 sum = add4(d2, add4(d3, d1));
    store4(dst, sum);
}

void test_twine()
{
    printf("---------- twine\n");
    double dst[12] = { 0 };
    const double src[12] = { 1.1, 1.2, 1.3, 2.1, 2.2, 2.3, 3.1, 3.2, 3.3, 4.1, 4.2, 4.3 };
    double X[4] = { 1 }, Y[4] = { 2 }, Z[4] = { 3 };
    untwine12(src, X, Y, Z);
    dump(4, X);
    dump(4, Y);
    dump(4, Z);
    twine12(X, Y, Z, dst);
    printf("twine12  "); dump(12, dst);
    repeat12(X, dst);
    printf("repeat12 "); dump(12, dst);
    twinedup12(src, dst);
    printf("twinedup "); dump(12, dst);
    sumXYZ(src, dst);
    printf("sumXYZ "); dump(4, dst);
    sumXXX(src, dst);
    printf("sumXXX "); dump(4, dst);
}


/**
 Change the data stride from 3 to 4, specifically:
 make dst = { XYZ? XYZ? XYZ? XYZ? }
 from src = { XYZ XYZ XYZ XYZ }
 */
inline void destride3x4(double const* src, double* dst)
{
    vec4 s0 = load4(src);
    vec4 s1 = load4(src+4);
    vec4 s2 = load4(src+8);
    store4(dst, s0);
    
    vec4 d1 = catshift2(s0, s1);
    d1 = shuffle4(d1, s1, 0b0101);
    store4(dst+4 , d1);
    
    vec4 d2 = blend22(s2, s1);
    d2 = swap2f128(d2);
    store4(dst+8 , d2);
    
    vec4 d3 = swap2f128(s2);
    d3 = shuffle4(s2, d3, 0b0101);
    store4(dst+12, d3);
}


void test_stride()
{
    printf("---------- stride\n");
    const size_t CNT = 1024;
    double a[CNT*12];
    double b[CNT*16];
    
    for ( size_t n = 0; n < 12*CNT; ++n )
        a[n] = 1 + n % 12;
    
    #pragma omp simd
    for ( size_t n = 0; n < CNT; ++n )
        destride3x4(a+12*n, b+16*n);
    
    for ( size_t n = 0; n < 4*CNT; ++n )
        b[4*n+3] = 0;
    
    dump(12, a);
    dump(16, b);
}

//------------------------------------------------------------------------------
#pragma mark -

/* transforms packed symmetric storage
 src = { a b c d e f }
 into { a b c 0 }, { b d e 0 }, { c e f 0 }
 */
void unpack_matrix()
{
    double src[] = { 1, 2, 3, 4, 5, 6 };
    
    vec2 zz = setzero2();
    vec2 ab = load2(src);
    vec2 cd = load2(src+2);
    vec2 ef = load2(src+4);
    
    vec2 m0h = blend11(cd, zz),  m0l = ab;
    vec2 m1h = blend11(ef, zz),  m1l = unpackhi2(ab, cd);
    vec2 m2h = catshift(ef, zz), m2l = unpacklo2(cd, ef);
    
    dump(cat22(m0l, m0h), "m0 ");
    dump(cat22(m1l, m1h), "m1 ");
    dump(cat22(m2l, m2h), "m2 ");
}


//------------------------------------------------------------------------------
#pragma mark -

void test_cat()
{
    printf("---------- cat\n");
    vec2 y{1.0, 2.0};
    vec2 x{3.0, 4.0};
    x = setr2(1.0, 2.0);
    y = setr2(3.0, 4.0);
    dump(x, "x");
    dump(y, "y");

    dump(cat22(y, x), "cat22(x, y)");
    dump(cat22(x, y), "cat22(y, x)");
}

void test_load()
{
    printf("---------- load\n");
    double mem[4] = { 1, 2, 3, 4 };
    vec4 x{1.0, 2.0, 3.0, 4.0};
    dump(x, "set ");
    
    vec4 y = load4(mem);
    dump(y, "load");
    
    // Attention: the upper 128-bits are in undefined state after a cast
    dump(cast4(load1(mem)), "cast4(load1)");
    dump(cast4(load2(mem)), "cast4(load2)");
    dump(load2Z(mem), "load2Z");

    dump(load3(mem), "load3");
    dump(cat22(load2(mem), load1(mem+2)), "cat22(load2, load1)");

    vec4 t = blend31(load4(mem), setzero4());
    dump(t, "blend(load4, zero)");


    vec2 u = load2(mem);
    vec4 a = broadcast2(mem);
    vec4 n = permute4(a, 0b1100);

    dump(u, "src");
    dump(n, "permute4(broadcast2)");
    dump(interleave2(u), "interleave2");
    dump(cat22(u,u), "cat22");
    dump(_mm256_set_m128d(u,u),    "set_m128d");
    dump(insertf128(cast4(u),u,1), "insertf128");
    dump(permute4(cat22(u,u), 0b1100), "permute4(cat22, 0b1100)");
    
    vec4 p{-1.0, -2.0, -3.0, -4.0};
    _mm_storeu_pd(mem, _mm256_castpd256_pd128(p));
    _mm_store_sd(mem+2,_mm256_extractf128_pd(p,1));
    dump(load4(mem), "load(store3)");
}

void test_load_float()
{
    printf("---------- load_float\n");
    float mem[4] = { 1, 2, 3, 4 };
    vec4f x{1.0, 2.0, 3.0, 4.0};
    dump(x, "set ");
    
    vec4f y = load4f(mem);
    dump(y, "load");
}

void test_broadcast()
{
    printf("---------- broadcast\n");
    double mem[4] = { 1, 2, 3, 4 };
    
    // using 4 loads
    dump(broadcast1(mem  ), "mem[0]");
    dump(broadcast1(mem+1), "mem[1]");
    dump(broadcast1(mem+2), "mem[2]");
    dump(broadcast1(mem+3), "mem[3]");
    
    // using 2 loads
    vec4 xy = broadcast2(mem);
    vec4 zt = broadcast2(mem+2);
    dump(duplo4(xy), " x");
    dump(duphi4(xy), " y");
    dump(duplo4(zt), " z");
    dump(duphi4(zt), " t");
    
    // using 1 load
    vec4 xyzt = load4(mem);
    xy = duplo2f128(xyzt);
    zt = duphi2f128(xyzt);
    dump(duplo4(xy), "X ");
    dump(duphi4(xy), "T ");
    dump(duplo4(zt), "Z ");
    dump(duphi4(zt), "T ");

    // using 2 permutes and 4 blends
    xyzt = load4(mem);
    vec4 u0 = unpacklo4(xyzt, xyzt);
    vec4 u1 = unpackhi4(xyzt, xyzt);
    vec4 v0 = swap2f128(u0);
    vec4 v1 = swap2f128(u1);
    dump(blend22(u0, v0), " x");
    dump(blend22(u1, v1), " y");
    dump(blend22(v0, u0), " z");
    dump(blend22(v1, u1), " t");

    // using 1 permute and 4 blends
    xyzt = load4(mem);
    vec4 ztxy = swap2f128(xyzt);
    u0 = unpacklo4(xyzt, xyzt);
    u1 = unpackhi4(xyzt, xyzt);
    v0 = unpacklo4(ztxy, ztxy);
    v1 = unpackhi4(ztxy, ztxy);
    //dump(u0, "u0");
    //dump(u1, "u1");
    //dump(v0, "v0");
    //dump(v1, "v1");
    dump(blend22(u0, v0), " x");
    dump(blend22(u1, v1), " y");
    dump(blend22(v0, u0), " z");
    dump(blend22(v1, u1), " t");

    // using 1 permute and 2 blends
    xyzt = load4(mem);
    ztxy = swap2f128(xyzt);
    xy = blend22(xyzt, ztxy);
    zt = blend22(ztxy, xyzt);
    dump(duplo4(xy), " x");
    dump(duphi4(xy), " y");
    dump(duplo4(zt), " z");
    dump(duphi4(zt), " t");
}

__m256i make_mask2(long i)
{
    switch( i )
    {
        case  0: return _mm256_setr_epi64x( 0, 0, 0, 0);
        case  1: return _mm256_setr_epi64x(-1, 0, 0, 0);
        case  2: return _mm256_setr_epi64x(-1,-1, 0, 0);
        case  3: return _mm256_setr_epi64x(-1,-1,-1, 0);
        default: return _mm256_setr_epi64x(-1,-1,-1,-1);
    }
}

__m256i make_mask(long i)
{
    vec4 v{0.5, 1.5, 2.5, 3.5};
    return _mm256_castpd_si256(cmp4(v, set4(i), _CMP_LT_OQ));
}

void test_store()
{
    printf("---------- store\n");
    double mem[4] = { 0, 0, 0, 0 };
    vec4 x{1.0, 2.0, 3.0, 4.0};
    dump(x, "value");
    for ( int i = 0; i < 5; ++i )
    {
        __m256i msk = make_mask(i);
        maskstore4(mem, msk, x);
        dump(load4(mem), "store");
    }
}


void test_swap1()
{
    printf("---------- swap1\n");
    vec4 a{ 1, 2, 3, 4};
    vec4 b{-1,-2,-3,-4};
    dump(a, "a = ");
    dump(b, "b = ");
    
    dump(permute4(a,0x05), "permute(a,0x05)");
    dump(permute4(b,0x05), "permute(b,0x05)");
    
    dump(permute2f128(a,b,0x20), "permute2f128(a,b,0x20)");
    dump(permute2f128(a,b,0x31), "permute2f128(a,b,0x31)");
    dump(permute2f128(a,b,0x01), "permute2f128(a,b,0x01)");
}


void test_rotate()
{
    printf("---------- rotate\n");
    vec4 a{ 1, 2, 3, 4 };
    vec4 b{-5,-6,-7,-8 };
    dump(a, "a/x = ");
    dump(b, "b/y = ");
    
    dump(catshift1(a, b), "catshift1(a, b)");
    dump(catshift2(a, b), "catshift2(a, b)");
    dump(catshift3(a, b), "catshift3(a, b)");
}


void test_swap2()
{
    printf("---------- swap2\n");
    vec4 a{ 1,  2,  3,  4};
    vec4 b{-1, -2, -3, -4};
    dump(a, "a");
    dump(b, "b");
    
    dump(permute2f128(a,a,0x08), "permute2f128(a,a,0x08)");
    dump(permute2f128(a,a,0x21), "permute2f128(a,a,0x21)");
    dump(permute2f128(a,a,0x81), "permute2f128(a,a,0x81)");
    
    dump(permute2f128(a,a,0x28), "permute2f128(a,a,0x28)");
    dump(permute2f128(a,a,0x81), "permute2f128(a,a,0x81)");
}


/*
 How to transform 'xyz' into 'xxx', 'yyy' and 'zzz'
 */
void test_swap4()
{
    printf("---------- swap4\n");
    vec4 s{1, 2, 3, 4};

    {
        vec4 z = unpacklo4(s, s);
        vec4 u = unpackhi4(s, s);
        
        dump(s, "src");
        dump(z, "lo ");
        dump(u, "hi ");
        
        dump(permute2f128(z, z, 0x00), "permute2f128(z,z) 0x00");
        dump(permute2f128(u, u, 0x00), "permute2f128(u,u) 0x00");
        dump(permute2f128(z, z, 0x11), "permute2f128(z,z) 0x11");
        dump(permute2f128(u, u, 0x11), "permute2f128(u,u) 0x11");
    }
    {
        vec4 p = permute2f128(s, s, 0x01);
        vec4 l = blend22(s, p);
        vec4 u = blend22(p, s);
        vec4 x0 = unpacklo4(l,l);
        vec4 x1 = unpackhi4(l,l);
        vec4 x2 = unpacklo4(u,u);
        vec4 x3 = unpackhi4(u,u);

        dump(s, "src");
        dump(x0, "xxxx");
        dump(x1, "yyyy");
        dump(x2, "zzzz");
        dump(x3, "tttt");
    }
    {
        vec4 p = permute2f128(s, s, 0x01);
        vec4 l = blend22(s, p);
        vec4 u = blend22(p, s);
        vec4 z = unpacklo4(u,u);
        
        dump(s, "xyzt");
        dump(l, "xyxy");
        dump(z, "zzzz");
    }
}

void test_hadd()
{
    printf("---------- hadd\n");
    vec4 a{1, -1, 2, -2};
    vec4 b{3, -3, 4, -4};
    dump(a, "a");
    dump(b, "b");

    vec4 p = permute2f128(a, b, 0x20);
    vec4 q = permute2f128(a, b, 0x31);
    dump(p, "p ");
    dump(q, "q ");

    vec4 z = unpacklo4(p, q);
    vec4 u = unpackhi4(p, q);

    dump(z, "z ");
    dump(u, "u ");
}

void test_hsum()
{
    /* finally sum horizontally:
     s0 = { Y0 Y0 Y0 - }, s1 = { Y1 Y1 Y1 - }, s2 = { Y2 Y2 Y2 - }
     to { Y0+Y0+Y0, Y1+Y1+Y1, Y2+Y2+Y2, 0 }
     */
    vec4f s0{ 1, 2, 3, 0.1 };
    vec4f s1{ 2, 3, 1, 0.1 };
    vec4f s2{ 3, 1, 2, 0.3 };
    vec4f s3{ 0, 0, 0, 0 };

    printf("---------- hsum\n");
    dump(s0, "s0  ");
    dump(s1, "s1  ");
    dump(s2, "s2  ");
    dump(s3, "s3  ");
    
    s0 = clear4th(s0); // clear garbage
    s1 = clear4th(s1); // clear garbage
    s2 = clear4th(s2); // clear garbage

    s0 = add4f(unpacklo4f(s0, s1), unpackhi4f(s0, s1));
    s2 = add4f(unpacklo4f(s2, s3), unpackhi4f(s2, s3));
    s0 = add4f(catshift2f(s0, s2), blend22f(s0, s2));
    dump(s0, "sum ");
    
    dump(blend0010f(s2, s0), "blend0010");
    dump(blend22f(s0, s2), "blend22");
}


//------------------------------------------------------------------------------
#pragma mark -

void test_mat()
{
    printf("---------- mat\n");
    // matrix:
    vec4 m012{0, 1, 2, -1};
    vec4 m345{3, 4, 5, -1};
    vec4 m678{6, 7, 8, -1};
    dump(m012, "m012");
    dump(m345, "m345");
    dump(m678, "m678");
    
    // symmetrized matrix:
    vec4 z = shuffle4(m012, m345, 0b0011);
    vec4 u = catshift2(z, m678);
    dump(z, "z");
    dump(u, "u");
    
    vec4 m145 = blend22(z, m345);
    vec4 m258 = blend22(u, m678);
    dump(m145, "m145");
    dump(m258, "m258");
    
    // transposed matrix:
    vec4 m036 = blend31(u, shuffle4(m012, m345, 0b1000));
    vec4 m147 = blend0010(m145, permute4(u, 0b0101));
    
    dump(m036, "m036");
    dump(m147, "m147");
}


void test_transpose2()
{
    printf("---------- transpose2\n");
    vec4 m{1, 2, 3, 4};
    vec4 t = blend0110(m, permute4(permute2f128(m,m,0x01),0b1100));
    dump(m, "m");
    dump(t, "t");
#if defined(__AVX2__)
    vec4 s = permute4x64(m, 0b11011000);
    dump(s, "s");
    dump(permute4x64(m, 0x88), "permute4x64(m, 0x88)");
    dump(permute4x64(m, 0xD8), "permute4x64(m, 0xD8)");
    dump(permute4x64(m, 0xDD), "permute4x64(m, 0xDD)");
    dump(permute4x64(m, 0x50), "permute4x64(m, 0x50)");
#endif
}


void test_transpose3()
{
    printf("---------- transpose3\n");
    // matrix:
    vec4 m012{1, 2, 3, 0};
    vec4 m345{4, 5, 6, 0};
    vec4 m678{7, 8, 9, 0};

    dump(m012, "m012");
    dump(m345, "m345");
    dump(m678, "m678");
    
    // symmetrized matrix:
    vec4 z = shuffle4(m012, m345, 0b0011);
    vec4 u = catshift2(z, m678);
    vec4 t = shuffle4(m012, m345, 0b1000);
    dump(z, "z");
    dump(u, "u");
    dump(t, "t");
    dump(shuffle4(u, m345, 0b1110), "tmp");
    // transposed matrix:
    vec4 m036 = blend0010(t, u);
    vec4 m147 = blend22(z, shuffle4(u, m345, 0b1100));
    vec4 m258 = blend22(u, m678);
    
    dump(m036, "m036");
    dump(m147, "m147");
    dump(m258, "m258");
}

void test_transpose4()
{
    printf("---------- transpose4\n");
    // matrix:
    vec4 m0{ 0,  1,  2,  3};
    vec4 m1{ 4,  5,  6,  7};
    vec4 m2{ 8,  9, 10, 11};
    vec4 m3{12, 13, 14, 15};

    dump(m0, "m0");
    dump(m1, "m1");
    dump(m2, "m2");
    dump(m3, "m3");
    printf("\n");

    // transpose all 2x2 subblocks:
    vec4 u0 = unpacklo4(m0, m1);
    vec4 u1 = unpackhi4(m0, m1);
    vec4 u2 = unpacklo4(m2, m3);
    vec4 u3 = unpackhi4(m2, m3);

    dump(u0, "u0");
    dump(u1, "u1");
    dump(u2, "u2");
    dump(u3, "u3");
    printf("\n");
    
    // using 4 permutes
    vec4 t0 = permute2f128(u0, u2, 0x20);
    vec4 t1 = permute2f128(u1, u3, 0x20);
    vec4 t2 = permute2f128(u0, u2, 0x31);
    vec4 t3 = permute2f128(u1, u3, 0x31);

    // transposed matrix:
    dump(t0, "t0");
    dump(t1, "t1");
    dump(t2, "t2");
    dump(t3, "t3");
    printf("\n");

    // using 2 permutes and 4 blend
    vec4 x02 = permute2f128(u0, u2, 0x21);
    vec4 x13 = permute2f128(u1, u3, 0x21);
    dump(x02, "x02");
    dump(x13, "x13");
    t0 = blend22(u0, x02);
    t1 = blend22(u1, x13);
    t2 = blend22(x02, u2);
    t3 = blend22(x13, u3);

    // transposed matrix:
    dump(t0, "t0");
    dump(t1, "t1");
    dump(t2, "t2");
    dump(t3, "t3");
}


/* return transposed matrix
make dst = { XXXX YYYY ZZZZ TTTT }
from src = { XYZT XYZT XYZT XYZT }
*/
void transpose4x4(double const* src, double* dst)
{
    vec4 u0 = loadu4(src);
    vec4 u1 = loadu4(src+4);
    vec4 v2 = loadu4(src+8);
    vec4 v3 = loadu4(src+12);
    vec4 v0 = unpacklo4(u0, u1);
    vec4 v1 = unpackhi4(u0, u1);
    u0 = unpacklo4(v2, v3);
    u1 = unpackhi4(v2, v3);
    v2 = catshift2(v0, u0);
    v3 = catshift2(v1, u1);
    storeu4(dst   , blend22(v0, v2));
    storeu4(dst+4 , blend22(v1, v3));
    storeu4(dst+8 , blend22(v2, u0));
    storeu4(dst+12, blend22(v3, u1));
}


void test_swap7()
{
    printf("---------- swap7\n");
    vec4 s{1, 2, 3, 4};
    dump(s, "source");

    dump(permute4(s, 0b1010), "permute 0b1010");
    dump(permute4(s, 0b0101), "permute 0b0101");

    dump(shuffle4(s, s, 0b1100), "shuffle 0b1100");
    dump(shuffle4(s, s, 0b0011), "shuffle 0b0011");
    
    dump(permute4(s, 0b1100), "permute 0b1100");
    dump(permute4(s, 0b0011), "permute 0b0011");

    dump(permute4(s, 0b0000), "permute 0b0000");
    dump(permute4(s, 0b1111), "permute 0b1111");

    dump(unpacklo4(s, s), "unpacklo");
    dump(unpackhi4(s, s), "unpackhi");

#if defined(__AVX2__)
    dump(permute4x64(s, 0xDD), "permute4x64 0xDD");
    dump(permute4x64(s, 0x88), "permute4x64 0x88");
    dump(permute4x64(s, 0xD8), "permute4x64 0xD8");
    dump(permute4x64(s, 0xC9), "permute4x64 0xC9");
    dump(permute4x64(s, 0xD2), "permute4x64 0xD2"); // Z X Y T
    dump(permute4x64(s, 0xC9), "permute4x64 0xC9"); // Y Z X T
#endif
}

#endif

/* return transposed matrix
make dst = { XXXX YYYY ZZZZ TTTT }
from src = { XYZT XYZT XYZT XYZT }
*/
void transpose4x4(float const* src, float* dst)
{
    vec4f v0 = loadu4f(src);
    vec4f v1 = loadu4f(src+4);
    vec4f u2 = loadu4f(src+8);
    vec4f u3 = loadu4f(src+12);
    vec4f u0 = unpacklo4f(v0, v1);
    vec4f u1 = unpackhi4f(v0, v1);
    v0 = unpacklo4f(u2, u3);
    v1 = unpackhi4f(u2, u3);
    storeu4f(dst   , movelh4f(u0, v0));
    storeu4f(dst+4 , movehl4f(v0, u0));
    storeu4f(dst+8 , movelh4f(u1, v1));
    storeu4f(dst+12, movehl4f(v1, u1));
}

#if defined(__ARM_NEON__)
void transpose4x4neon(float const* src, float* dst)
{
    float32x4_t v0 = vld1q_f32(src);
    float32x4_t v1 = vld1q_f32(src+4);
    float32x4_t v2 = vld1q_f32(src+8);
    float32x4_t v3 = vld1q_f32(src+12);
    float32x4_t x = vtrn1q_f32(v0, v1);
    float32x4_t y = vtrn2q_f32(v0, v1);
    float32x4_t z = vtrn1q_f32(v2, v3);
    float32x4_t t = vtrn2q_f32(v2, v3);
#if 1
    vst1q_f32(dst   , vtrn1q_f64(x, z));
    vst1q_f32(dst+4 , vtrn1q_f64(y, t));
    vst1q_f32(dst+8 , vtrn2q_f64(x, z));
    vst1q_f32(dst+12, vtrn2q_f64(y, t));
#else
    vst1q_f32(dst   , movelh4f(x, z));
    vst1q_f32(dst+4 , movelh4f(y, t));
    vst1q_f32(dst+8 , movehl4f(z, x));
    vst1q_f32(dst+12, movehl4f(t, y));
#endif
}
#endif

void test_transpose2x2()
{
    printf("---------- transpose2x2\n");
    vec4f a{1,2,3,4};
    vec4f b{5,6,7,8};
    dump(a, "a");
    dump(b, "b");
    dump(movehl4f(a,b), "movehl4f(a,b)");
    dump(movelh4f(a,b), "movelh4f(a,b)");
}

void test_transpose4x4()
{
    //const double src[16] = { 1.1, 1.2, 1.3, 1.4, 2.1, 2.2, 2.3, 2.4, 3.1, 3.2, 3.3, 3.4, 4.1, 4.2, 4.3, 4.4 };
    const float src[16] = { 1.1, 2.1, 3.1, 4.1, 1.2, 2.2, 3.2, 4.2, 1.3, 2.3, 3.3, 4.3, 1.4, 2.4, 3.4, 4.4 };
    printf("---------- transpose4x4\n");
    {
        float dst[16] = { 0 };
        transpose4x4(src, dst);
        printf("float\n"); dump(4, 4, src);
        printf("transpose\n"); dump(4, 4, dst);
    }
#if defined(__ARM_NEON__)
    {
        float dst[16] = { 0 };
        transpose4x4neon(src, dst);
        printf("Transpose\n"); dump(4, 4, dst);
    }
#endif
#if defined(__AVX__)
    {
        double dst[16] = { 0 };
        //const double src[16] = { 1.1, 1.2, 1.3, 1.4, 2.1, 2.2, 2.3, 2.4, 3.1, 3.2, 3.3, 3.4, 4.1, 4.2, 4.3, 4.4 };
        const double src[16] = { 1.1, 2.1, 3.1, 4.1, 1.2, 2.2, 3.2, 4.2, 1.3, 2.3, 3.3, 4.3, 1.4, 2.4, 3.4, 4.4 };
        transpose4x4(src, dst);
        printf("double\n"); dump(4, 4, src);
        printf("transpose\n"); dump(4, 4, dst);
    }
#endif
}

#pragma mark - ARM NEON loads


void test_loads()
{
        printf("---------- loads\n");
    double src[] = { 1, 2, 3, 4, 5, 6 };
    
#if defined(__ARM_NEON__)
    float64x2x3_t tmp = vld3q_f64(src);
    dump(tmp.val[0], "vld3[0] ");
    dump(tmp.val[1], "vld3[1] ");
    dump(tmp.val[2], "vld3[2] ");
        
    float64x2_t one = vld1q_f64(src);
    dump(one, "vld1 ");

    float flt[] = { 1, 2, 3, 4, 5, 6 };
    float32x2_t two = vld1_f32(flt);
    dump(two, "vldf ");
#endif
}


    

int main(int argc, char * argv[])
{
    int i = 0;
    if ( argc > 1 )
        i = std::max(0, atoi(argv[1]));

    switch ( i )
    {
        case 0:
            test_swapSSE();
            test_loads();
            test_transpose2x2();
            test_transpose4x4();
            test_signselect();
            test_float4();
            break;
#ifdef __AVX__
        case 1:
            unpack_matrix();
            test_twine();
            test_stride();
            break;
        case 2:
            test_cat();
            test_load();
            test_load_float();
            test_broadcast();
            test_store();
            break;
        case 3:
            test_rotate();
            test_hadd();
            test_hsum();
            break;
        case 4:
            test_swap1();
            test_swap2();
            test_swap4();
            break;
        case 5:
            test_transpose2();
            test_transpose3();
            test_transpose4();
            break;
        case 6:
            test_shuffle();
            test_swap7();
            break;
#else
        default:
            printf("AVX was not enabled\n");
#endif
    }
}

