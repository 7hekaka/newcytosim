// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <sys/time.h>
#include <limits>
#include <cfenv>

#define DIM 3

#include "assert_macro.h"
#include "real.h"
#include "vector.h"
#include "random.h"
#include "vecprint.h"
#include "blas.h"
#include "lapack.h"
#include "xpttrf.h"
#include "cytoblas.h"
#include "simd.h"
#include "simd_float.h"
#include "simd_print.h"


/// keeping time using Intel's cycle counters
unsigned long long rdt = 0;
/// start timer
inline void tic() { rdt = __rdtsc(); }
/// stop timer and print time
inline double toc(double num) { return double(__rdtsc()-rdt)/num; }


/// type used for indexing
typedef size_t SIZE_T;

/// number of segments:
const SIZE_T NSEG = 127;
/// number of vertex coordinates
const SIZE_T NVAL = DIM * ( NSEG + 1 );
const SIZE_T ALOC = NVAL + 8;

const SIZE_T DISP = 24UL;


/// fill vector with signalling NANs
void nan_fill(SIZE_T cnt, real * ptr)
{
    real n = std::numeric_limits<real>::signaling_NaN();
    for ( SIZE_T u = 0; u < cnt; ++u )
        ptr[u] = n;
}


void new_nans(SIZE_T cnt, real*& ptr)
{
    free_real(ptr);
    ptr = new_real(cnt);
    nan_fill(cnt, ptr);
    //printf("%p = new_nans(%lu)\n", ptr, cnt);
}

void new_nans(SIZE_T cnt, real*& x, real*& y, real*& z)
{
    new_nans(cnt, x);
    new_nans(cnt, y);
    new_nans(cnt, z);
}

void randomize(SIZE_T cnt, real*& x, real*& y, real*& z, real mag)
{
    for ( SIZE_T i = 0; i < cnt; ++i )
    {
        x[i] = mag * RNG.sreal();
        y[i] = mag * RNG.sreal();
        z[i] = mag * RNG.sreal();
    }
}

void free_reals(real*& x, real*& y, real*& z)
{
    free_real(x); x = nullptr;
    free_real(y); y = nullptr;
    free_real(z); z = nullptr;
}

void free_reals(real*& x, real*& y, real*& z, real*& t)
{
    free_real(x); x = nullptr;
    free_real(y); y = nullptr;
    free_real(z); z = nullptr;
    free_real(t); t = nullptr;
}


#pragma STDC_FENV_ACCESS on
void print_fe_exceptions(const char* str)
{
    int n = std::fetestexcept(FE_ALL_EXCEPT);
    n &= ~FE_INEXACT;
    if ( n )
    {
        fprintf(stderr, " %s(", str);
        if (n&FE_DIVBYZERO)  fprintf(stderr, "DIVBYZERO ");
        if (n&FE_INEXACT)    fprintf(stderr, "INEXACT ");
        if (n&FE_INVALID)    fprintf(stderr, "INVALID ");
        if (n&FE_OVERFLOW)   fprintf(stderr, "OVERFLOW ");
        if (n&FE_UNDERFLOW)  fprintf(stderr, "UNDERFLOW ");
        if (n&FE_ALL_EXCEPT) fprintf(stderr, "UNKNOWN ");
        if ( std::feclearexcept(FE_ALL_EXCEPT) )
            fprintf(stderr, "unclear ");
        fprintf(stderr, "\b)");
    }
}



#ifdef __AVX__

inline __m256d loadc4(double const* ptr)
{
#if ( 0 )
    // we can check memory alignment here:
    uintptr_t x = ((uintptr_t)ptr) & 31;
    if ( x )
    {
        fprintf(stderr, "loading unaligned memory %i\n", x);
        return _mm256_loadu_pd(ptr);
    }
#endif
    return _mm256_load_pd(ptr);
}

#endif


inline void print_alignment(real const* ptr, const char msg[])
{
    fprintf(stderr, "%s %p alignment %lu\n", msg, ptr, (uintptr_t)ptr&31);
}

//------------------------------------------------------------------------------

/// vectors for one filament
real *pos_=nullptr, *dir_=nullptr, *ani_=nullptr;
real *lag_=nullptr, *force_=nullptr, *tmp_=nullptr;
real *diag_=nullptr, *upper_=nullptr;


inline void DPTTRF(SIZE_T nbs)
{
    int info = 0;
    alsatian_xpttrf(nbs, diag_, upper_, &info);
}

inline void DPTTS2(SIZE_T nbs)
{
    alsatian_xptts2(nbs, 1, diag_, upper_, lag_, 0);
}

//------------------------------------------------------------------------------

void setFilament(SIZE_T nbs, real seg, real persistence_length, real mag)
{
    nbs = std::min(nbs, NSEG);
    SIZE_T nbv = DIM * ( nbs + 1 );
    new_nans(nbv, pos_);
    new_nans(nbv, dir_);
    new_nans(nbv, lag_);
    new_nans(nbv, force_);

    real sigma = std::sqrt(2.0*seg/persistence_length);
    
    Vector pos(0,0,0);
    Vector dir = Vector::randU();
    
    pos.store(pos_);
    dir.store(dir_);
    for ( SIZE_T p = 1 ; p < nbs; ++p )
    {
        pos += seg * dir;
        pos.store(pos_+DIM*p);
        dir.store(dir_+DIM*p);
        //rotate dir in a random direction:
        real a = sigma * RNG.gauss();
        dir = std::cos(a) * dir + dir.randOrthoU(std::sin(a));
    }
    pos.store(pos_+DIM*nbs);
    
    for ( SIZE_T i = 0 ; i < nbv; ++i )
        force_[i] = mag * RNG.sreal();
    
    print_fe_exceptions("setFilament");
}


void setProjection(SIZE_T nbs)
{
    new_nans(nbs, diag_);
    new_nans(nbs, upper_);
    
    SIZE_T j = 0;
    for ( ; j < nbs-1; ++j )
    {
        const real* X = dir_ + DIM * j;
#if ( DIM == 2 )
        upper_[j] = -( X[0]*X[2] + X[1]*X[3] );
#else
        upper_[j] = -( X[0]*X[3] + X[1]*X[4] + X[2]*X[5] );
#endif
        diag_[j] = 2.0;
    }
    diag_[j] = 2.0;

    DPTTRF(nbs);
#if ( 0 )
    std::clog << "\nD "; VecPrint::print(std::clog, nbs, diag_, 3);
    std::clog << "\nU "; VecPrint::print(std::clog, nbs-1, upper_, 3); std::clog << '\n';
#endif
    print_fe_exceptions("setProjection");
}


void setAnisotropy(SIZE_T nbs)
{
    const SIZE_T sup = DIM * nbs;
    new_nans(sup+DIM, ani_);
    new_nans(sup+DIM, tmp_);

    // for the extremities, the direction of the nearby segment is used.
    for ( SIZE_T d = 0; d < DIM; ++d )
    {
        ani_[d]     = dir_[d];
        ani_[d+sup] = dir_[d+sup-DIM];
    }
    // for intermediate points, the directions of the flanking segments are averaged
    for ( SIZE_T p = DIM ; p < sup; ++p )
        ani_[p] = 0.5 * ( dir_[p-DIM] + dir_[p] );
    
    print_fe_exceptions("setAnisotroy");
}


void setRandom(int np, real * vec, real mag)
{
    for ( SIZE_T p = 0 ; p < ALOC; ++p )
        vec[p] = mag * RNG.sreal();
}

//------------------------------------------------------------------------------
#pragma mark - PROJECT UP

void projectForcesU_(SIZE_T nbs, const real* dif, const real* src, real* mul)
{
    #pragma vector unaligned
    for ( SIZE_T i = 0; i < nbs; ++i )
    {
        const real * X = src + DIM * i;
        const real * d = dif + DIM * i;
        mul[i] = d[0] * ( X[DIM  ] - X[0] )
               + d[1] * ( X[DIM+1] - X[1] )
#if ( DIM > 2 )
               + d[2] * ( X[DIM+2] - X[2] )
#endif
        ;
    }
}

void projectForcesU_PTR(SIZE_T nbs, const real* dif, const real* src, real* mul)
{
    const real *const end = mul + nbs;

    while ( mul < end )
    {
        *mul = dif[0] * ( src[DIM  ] - src[0] )
             + dif[1] * ( src[DIM+1] - src[1] )
#if ( DIM > 2 )
             + dif[2] * ( src[DIM+2] - src[2] )
#endif
        ;
        src += DIM;
        dif += DIM;
        ++mul;
    }
}

/**
 Perform first calculation needed by projectForces:
 */
inline void projectForcesU_TWO(SIZE_T nbs, const real* dif, const real* src, real* mul)
{
    real x3, x0 = src[0];
    real x4, x1 = src[1];
#if ( DIM >= 3 )
    real x5, x2 = src[2];
#endif
    src += DIM;
    real *const end = mul + nbs;
    
    //further optimization with manual loop-unrolling
    if ( nbs & 1 )
    {
        x3 = src[0];
        x4 = src[1];
#if ( DIM == 2 )
        mul[0] = dif[0] * (x3 - x0) + dif[1] * (x4 - x1);
#elif ( DIM >= 3 )
        x5 = src[2];
        mul[0] = dif[0] * (x3 - x0) + dif[1] * (x4 - x1) + dif[2] * (x5 - x2);
        x2 = x5;
#endif
        ++mul;
        src += DIM;
        dif += DIM;
        x0 = x3;
        x1 = x4;
    }
    
    while ( mul < end )
    {
        x3 = src[0];
        x4 = src[1];
#if ( DIM == 2 )
        mul[0] = dif[0] * (x3 - x0) + dif[1] * (x4 - x1);
#elif ( DIM >= 3 )
        x5 = src[2];
        mul[0] = dif[0] * (x3 - x0) + dif[1] * (x4 - x1) + dif[2] * (x5 - x2);
#endif
        
#if ( DIM == 2 )
        x0 = src[2];
        x1 = src[3];
        mul[1] = dif[2] * (x0 - x3) + dif[3] * (x1 - x4);
#elif ( DIM >= 3 )
        x0 = src[3];
        x1 = src[4];
        x2 = src[5];
        mul[1] = dif[3] * (x0 - x3) + dif[4] * (x1 - x4) + dif[5] * (x2 - x5);
#endif
        
        mul += 2;
        src += 2*DIM;
        dif += 2*DIM;
    }
    assert_true( mul == end );
}


#if REAL_IS_DOUBLE && defined(__SSE3__)

/**
 Perform first calculation needed by projectForces:
 tmp <- J * X
 */
void projectForcesU2D_SSE(SIZE_T nbs, const real* dif, const real* src, real* mul)
{
    real const*const end = mul - 1 + nbs;

    vec2 y, x = load2(src);
    while ( mul < end )
    {
        y = load2(src+2);
        vec2 a = mul2(sub2(y, x), load2(dif));
        x = load2(src+4);
        src += 4;
        vec2 b = mul2(sub2(x, y), load2(dif+2));
        dif += 4;
        //storeu2(mul, hadd2(a, b));
        storeu2(mul, add2(unpacklo2(a, b), unpackhi2(a, b)));
        mul += 2;
    }
    
    if ( mul < end+1 )
    {
        y = load2(src+2);
        vec2 a = mul2(sub2(y, x), load2(dif));
        //store1(mul, hadd2(a, a));
        store1(mul, add2(a, unpackhi2(a, a)));
    }
}

#endif

#if REAL_IS_DOUBLE && defined(__AVX__)

__m256i make_mask(long i)
{
    vec4 v{0.5, 1.5, 2.5, 3.5};
    return _mm256_castpd_si256(cmp4(v, set4(i), _CMP_LT_OQ));
}

/**
 Perform first calculation needed by projectForces:
 tmp <- J * X
 F. Nedelec, 9.12.2016
 */
void projectForcesU2D_AVX(SIZE_T nbs, const real* dif, const real* src, real* mul)
{
    real const*const end = mul - 7 + nbs;
    
    while ( mul < end )
    {
        vec4 a = mul4(sub4(loadu4(src+2 ), loadu4(src   )), loadc4(dif   ));
        vec4 b = mul4(sub4(loadu4(src+6 ), loadu4(src+4 )), loadc4(dif+4 ));
        vec4 c = mul4(sub4(loadu4(src+10), loadu4(src+8 )), loadc4(dif+8 ));
        vec4 d = mul4(sub4(loadu4(src+14), loadu4(src+12)), loadc4(dif+12));
        src += 16;
        dif += 16;
        vec4 p = permute2f128(a,b,0x20);
        vec4 q = permute2f128(a,b,0x31);
        vec4 s = permute2f128(c,d,0x20);
        vec4 t = permute2f128(c,d,0x31);
        store4(mul  , add4(unpacklo4(p, q), unpackhi4(p, q)));
        store4(mul+4, add4(unpacklo4(s, t), unpackhi4(s, t)));
        mul += 8;
    }

    while ( mul < end+4 )
    {
        vec4 a = mul4(sub4(loadu4(src+2), loadu4(src  )), loadc4(dif  ));
        vec4 b = mul4(sub4(loadu4(src+6), loadu4(src+4)), loadc4(dif+4));
        dif += 8;
        src += 8;
        vec4 p = permute2f128(a,b,0x20);
        vec4 q = permute2f128(a,b,0x31);
        store4(mul, add4(unpacklo4(p, q), unpackhi4(p, q)));
        mul += 4;
    }

    while ( mul < end+6 )
    {
        //mul[jj] = dif[0] * ( X[DIM] - X[0] ) + dif[1] * ( X[DIM+1] - X[1] )
        vec4 d = mul4(sub4(loadu4(src+2), loadu4(src)), load4(dif));
        src += 4;
        dif += 4;
        vec2 h = gethi(d);
        storeu2(mul, add2(unpacklo2(getlo(d),h), unpackhi2(getlo(d),h)));
        mul += 2;
    }
    
    if ( mul < end+7 )
    {
        vec2 a = mul2(sub2(load2(src+2), load2(src)), load2(dif));
        store1(mul, add2(a, permute2(a,1)));
    }
}


/**
 Attention: this does not check the boundaries and will write beyond the
 nbs-th point of tmp, which should be allocated accordingly.
 F. Nedelec, 11.01.2018
 */
void projectForcesU2D_AVY(SIZE_T nbs, const real* dif, const real* src, real* mul)
{
    real *const end = mul - 3 + nbs;
    
    // calculate the terms 4 by 4
    while ( mul < end )
    {
        vec4 a = mul4(sub4(loadu4(src+2), loadu4(src  )), loadc4(dif  ));
        vec4 b = mul4(sub4(loadu4(src+6), loadu4(src+4)), loadc4(dif+4));
        dif += 8;
        src += 8;
        vec4 p = permute2f128(a, b, 0x20);
        vec4 q = permute2f128(a, b, 0x31);
        store4(mul, add4(unpacklo4(p, q), unpackhi4(p, q)));
        mul += 4;
    }
}

#endif

//------------------------------------------------------------------------------
#pragma mark - PROJECT UP 3D

#if REAL_IS_DOUBLE && defined(__SSE3__)

void projectForcesU3D_SSE(SIZE_T nbs, const real* dif, const real* src, real* mul)
{
    const real *const end = mul + nbs - 1;
#if 1
    // unrolled bulk of the calculation, processing 4x3 scalars
    while ( mul < end-2 )
    {
        /*
         *mul = dif[0] * ( src[DIM  ] - src[0] )
              + dif[1] * ( src[DIM+1] - src[1] )
              + dif[2] * ( src[DIM+2] - src[2] );
         */
        vec2 s0 = mul2(load2(dif  ), sub2(loadu2(src+3), loadu2(src  )));
        vec2 s1 = mul2(load2(dif+2), sub2(loadu2(src+5), loadu2(src+2)));
        vec2 s2 = mul2(load2(dif+4), sub2(loadu2(src+7), loadu2(src+4)));
        vec2 ss = unpacklo2(s0, s2);
        vec2 tt = unpackhi2(s0, s2);

        vec2 t0 = mul2(load2(dif+ 6), sub2(loadu2(src+ 9), loadu2(src+ 6)));
        vec2 t1 = mul2(load2(dif+ 8), sub2(loadu2(src+11), loadu2(src+ 8)));
        vec2 t2 = mul2(load2(dif+10), sub2(loadu2(src+13), loadu2(src+10)));
        vec2 xx = unpacklo2(t0, t2);
        vec2 yy = unpackhi2(t0, t2);

        store2(mul, add2(tt, add2(s1, ss)));
        store2(mul+2, add2(yy, add2(t1, xx)));

        src += 12;
        dif += 12;
        mul += 4;
    }
#endif
    // bulk of the calculation, processing 2x3 scalars
    while ( mul < end )
    {
        /*
         *mul = dif[0] * ( src[DIM  ] - src[0] )
              + dif[1] * ( src[DIM+1] - src[1] )
              + dif[2] * ( src[DIM+2] - src[2] );
         */
        vec2 s0 = mul2(load2(dif  ), sub2(loadu2(src+3), loadu2(src  )));
        vec2 s1 = mul2(load2(dif+2), sub2(loadu2(src+5), loadu2(src+2)));
        vec2 s2 = mul2(load2(dif+4), sub2(loadu2(src+7), loadu2(src+4)));

        vec2 ss = unpacklo2(s0, s2);
        vec2 tt = unpackhi2(s0, s2);

        store2(mul, add2(tt, add2(s1, ss)));
        
        src += 6;
        dif += 6;
        mul += 2;
    }
    // process remaining vector, 3 scalars
    if ( mul <= end )
    {
         /*
          *mul = dif[0] * ( src[DIM  ] - src[0] )
               + dif[1] * ( src[DIM+1] - src[1] )
               + dif[2] * ( src[DIM+2] - src[2] );
          */
         vec2 s0 = mul2(load2(dif  ), sub2(loadu2(src+3), loadu2(src)));
         vec2 s1 = mul1(load1(dif+2), sub1(load1(src+5), load1(src+2)));
         vec2 ss = unpackhi2(s0, setzero2());

         store1(mul, add1(s0, add1(s1, ss)));
   }
}

#endif

#if REAL_IS_DOUBLE && defined(__AVX__)

inline void twine3x4(real const* X, real const* Y, real const* Z, real* dst)
{
    vec4 sx = load4(X);
    vec4 sy = load4(Y);
    vec4 sz = load4(Z);

    vec4 zx = blend4(sx, sz, 0b0101);
    zx = swap2f128(zx);
    vec4 xy = unpacklo4(sx, sy);
    vec4 yz = unpackhi4(sy, sz);
    
    storeu4(dst  , blend22(xy, zx));
    storeu4(dst+4, blend22(yz, xy));
    storeu4(dst+8, blend22(zx, yz));
}

/**
 make
     dX = { XXXX }
     dY = { YYYY }
     dZ = { ZZZZ }
 from src = { XYZ XYZ XYZ XYZ }
 */
inline void untwine4x3(real const* src, real* X, real* Y, real* Z)
{
    vec4 s0 = load4(src);
    vec4 s1 = load4(src+4);
    vec4 s2 = load4(src+8);

    vec4 zx = blend22(s2, s0);
    vec4 xy = blend22(s0, s1);
    zx = swap2f128(zx);
    vec4 yz = blend22(s1, s2);
    
    store4(X,   blend4(zx, xy, 0b0101));
    store4(Y, shuffle4(xy, yz, 0b0101));
    store4(Z,   blend4(zx, yz, 0b1010));
}


void projectForcesU3D_AVX(SIZE_T nbs, const real* dif, const real* src, real* mul)
{
    const real *const end = mul - 3 + nbs;
    // bulk of the calculation, processing 4x3 scalars
    while ( mul < end )
    {
        /*
         *mul = dif[0] * ( src[DIM  ] - src[0] )
              + dif[1] * ( src[DIM+1] - src[1] )
              + dif[2] * ( src[DIM+2] - src[2] );
         */
        vec4 s0 = mul4(load4(dif  ), sub4(loadu4(src+ 3), loadu4(src  )));
        vec4 s1 = mul4(load4(dif+4), sub4(loadu4(src+ 7), loadu4(src+4)));
        vec4 s2 = mul4(load4(dif+8), sub4(loadu4(src+11), loadu4(src+8)));

        vec4 zx = blend22(s2, s0);
        vec4 xy = blend22(s0, s1);
        zx = swap2f128(zx);
        vec4 yz = blend22(s1, s2);
        
        vec4 mm = shuffle4(xy, yz, 0b0101);
        zx = add4(blend4(zx, yz, 0b1010), blend4(zx, xy, 0b0101));
        store4(mul, add4(mm, zx));
        
        src += 12;
        dif += 12;
        mul += 4;
    }
    while ( mul < end+2 )
    {
        vec4 s0 = mul4(load4(dif  ), sub4(loadu4(src+3), loadu4(src)));
        vec2 s1 = mul2(load2(dif+4), sub2(loadu2(src+7), loadu2(src+4)));

        vec2 xy = getlo(s0);
        vec2 zx = gethi(s0);

        vec2 mm = gethilo2(xy, s1);
        zx = add2(blend11(zx, s1), blend11(xy, zx));
        store2(mul, add2(mm, zx));

        src += 6;
        dif += 6;
        mul += 2;
    }
    while ( mul < end+3 )
    {
        vec2 x = mul2(loadu2(dif), sub2(loadu2(src+3), loadu2(src)));
        vec2 z = mul1(load1(dif+2), sub1(load1(src+5), load1(src+2)));
        vec2 y = permute2(x, 0b1);
        
        store1(mul, add1(add1(x, y), z));
        
        src += 3;
        dif += 3;
        ++mul;
    }
    assert_true(mul==end+3);
}

#endif

template <void (*FUNC)(SIZE_T, const real*, const real*, real*)>
void testU(SIZE_T cnt, char const* str)
{
    real *x = nullptr, *y = nullptr, *z = nullptr;
    new_nans(ALOC, x, y, z);
    randomize(NVAL, x, y, z, 1.0);
    nan_fill(NVAL, lag_);

    FUNC(NSEG, dir_, force_, lag_);
    VecPrint::print(std::cout, std::min(DISP,NSEG), lag_);
    std::cout << " |";
    VecPrint::print(std::cout, 1, lag_+NSEG);
    tic();
    for ( SIZE_T i=0; i<cnt; ++i )
    {
        FUNC(NSEG, dir_, y, z);
        // check the code with unaligned memory:
        FUNC(NSEG, dir_, x+DIM, y);
        FUNC(NSEG, dir_, z, y);
    }
    printf("  %10s %5.2f\n", str, toc(cnt*NSEG));
    
    free_reals(x,y,z);
}


void testProjectionU(SIZE_T cnt)
{
    std::cout << "testProjection UP " << DIM << "D, " << NSEG << " segments\n";
    testU<projectForcesU_>(cnt,    " U_   ");
    testU<projectForcesU_PTR>(cnt, " U_PTR");
    testU<projectForcesU_TWO>(cnt, " U_TWO");
#if ( DIM == 2 ) && REAL_IS_DOUBLE && defined(__SSE__)
    testU<projectForcesU2D_SSE>(cnt, " U_SSE");
#endif
#if ( DIM == 2 ) && REAL_IS_DOUBLE && defined(__AVX__)
    testU<projectForcesU2D_AVX>(cnt, " U_AVX");
    testU<projectForcesU2D_AVY>(cnt, " U_AVY");
#endif
#if ( DIM == 3 ) && REAL_IS_DOUBLE && defined(__SSE__)
    testU<projectForcesU3D_SSE>(cnt, " U_SSE");
#endif
#if ( DIM == 3 ) && REAL_IS_DOUBLE && defined(__AVX__)
    testU<projectForcesU3D_AVX>(cnt, " U_AVX");
#endif
}


//------------------------------------------------------------------------------
#pragma mark - PROJECT DOWN

/**
 Perform second calculation needed by projectForces:
 Y <- X + Jt * tmp
 */
void projectForcesD_(SIZE_T nbs, const real* dif, const real* X, const real* mul, real* Y)
{
    for ( SIZE_T d = 0, e = DIM*nbs; d < DIM; ++d, ++e )
    {
        Y[d] = X[d] + dif[d    ] * mul[    0];
        Y[e] = X[e] - dif[e-DIM] * mul[nbs-1];
    }
    
    for ( SIZE_T jj = 1; jj < nbs; ++jj )
    {
        const SIZE_T kk = DIM*jj;
        Y[kk  ] = X[kk  ] + dif[kk  ] * mul[jj] - dif[kk-DIM  ] * mul[jj-1];
        Y[kk+1] = X[kk+1] + dif[kk+1] * mul[jj] - dif[kk-DIM+1] * mul[jj-1];
#if ( DIM > 2 )
        Y[kk+2] = X[kk+2] + dif[kk+2] * mul[jj] - dif[kk-DIM+2] * mul[jj-1];
#endif
    }
}


/**
 1xMUL 2xADD
 */
void projectForcesD_ADD(SIZE_T nbs, const real* dif, const real* X, const real* mul, real* Y)
{
    real a0 = X[0];
    real a1 = X[1];
#if ( DIM > 2 )
    real a2 = X[2];
#endif
    
    for ( SIZE_T jj = 0; jj < nbs; ++jj )
    {
        real b0 = dif[DIM*jj  ] * mul[jj];
        real b1 = dif[DIM*jj+1] * mul[jj];
#if ( DIM > 2 )
        real b2 = dif[DIM*jj+2] * mul[jj];
#endif
        Y[DIM*jj  ] = a0 + b0;
        Y[DIM*jj+1] = a1 + b1;
#if ( DIM > 2 )
        Y[DIM*jj+2] = a2 + b2;
#endif
        a0 = X[DIM*jj+ DIM   ] - b0;
        a1 = X[DIM*jj+(DIM+1)] - b1;
#if ( DIM > 2 )
        a2 = X[DIM*jj+(DIM+2)] - b2;
#endif
    }
    
    const SIZE_T ee = DIM * nbs;
    Y[ee  ] = a0;
    Y[ee+1] = a1;
#if ( DIM > 2 )
    Y[ee+2] = a2;
#endif
}


/// proceeding downward, this could be fused with downward part of DPTTS2
void projectForcesD_REV(SIZE_T nbs, const real* dif, const real* X, const real* mul, real* Y)
{
    const real * E = X + DIM * nbs;
    real a0 = E[0];
    real a1 = E[1];
#if ( DIM > 2 )
    real a2 = E[2];
#endif
    
    for ( SIZE_T i = nbs; i-- > 0; )
    {
        real m0 = dif[DIM*i  ] * mul[i];
        real m1 = dif[DIM*i+1] * mul[i];
#if ( DIM > 2 )
        real m2 = dif[DIM*i+2] * mul[i];
#endif
        Y[DIM*(i+1)  ] = a0 - m0;
        Y[DIM*(i+1)+1] = a1 - m1;
#if ( DIM > 2 )
        Y[DIM*(i+1)+2] = a2 - m2;
#endif
        a0 = X[DIM*i  ] + m0;
        a1 = X[DIM*i+1] + m1;
#if ( DIM > 2 )
        a2 = X[DIM*i+2] + m2;
#endif
    }
    
    Y[0] = a0;
    Y[1] = a1;
#if ( DIM > 2 )
    Y[2] = a2;
#endif
}


/// proceeding downward with pointers, this could be fused with downward part of DPTTS2
void projectForcesD_RIV(SIZE_T nbs, const real* dif, const real* X, const real* mul, real* Y)
{
    const real * E = X + DIM * nbs;
    real a0 = E[0];
    real a1 = E[1];
#if ( DIM > 2 )
    real a2 = E[2];
#endif

    mul += nbs-1;
    dif += DIM*(nbs-1);
    Y += DIM*nbs;
    
    while ( E > X )
    {
        E -= DIM;
        real m0 = dif[0] * mul[0];
        real m1 = dif[1] * mul[0];
#if ( DIM > 2 )
        real m2 = dif[2] * mul[0];
#endif
        --mul;
        dif -= DIM;
        Y[0] = a0 - m0;
        Y[1] = a1 - m1;
#if ( DIM > 2 )
        Y[2] = a2 - m2;
#endif
        Y -= DIM;

        a0 = E[0] + m0;
        a1 = E[1] + m1;
#if ( DIM > 2 )
        a2 = E[2] + m2;
#endif
    }
    
    Y[0] = a0;
    Y[1] = a1;
#if ( DIM > 2 )
    Y[2] = a2;
#endif
}



/**
 2xFMA
 */
void projectForcesD_FMA(SIZE_T nbs, const real* dif, const real* X, const real* mul, real* Y)
{
    real a0 = X[0];
    real a1 = X[1];
#if ( DIM > 2 )
    real a2 = X[2];
#endif
    
    for ( SIZE_T jj = 0; jj < nbs; ++jj )
    {
        real m = mul[jj];
        real d0 = dif[DIM*jj  ];
        real d1 = dif[DIM*jj+1];
#if ( DIM > 2 )
        real d2 = dif[DIM*jj+2];
#endif
        Y[DIM*jj  ] = a0 + d0 * m;
        Y[DIM*jj+1] = a1 + d1 * m;
#if ( DIM > 2 )
        Y[DIM*jj+2] = a2 + d2 * m;
#endif
        a0 = X[DIM*jj+ DIM   ] - d0 * m;
        a1 = X[DIM*jj+(DIM+1)] - d1 * m;
#if ( DIM > 2 )
        a2 = X[DIM*jj+(DIM+2)] - d2 * m;
#endif
    }
    
    const SIZE_T ee = DIM * nbs;
    Y[ee  ] = a0;
    Y[ee+1] = a1;
#if ( DIM > 2 )
    Y[ee+2] = a2;
#endif
}


/**
 Perform second calculation needed by projectForces:
 */
void projectForcesD_PTR(SIZE_T nbs, const real* dif,
                        const real* X, const real* mul, real* Y)
{
    real a0 = X[0];
    real a1 = X[1];
#if ( DIM > 2 )
    real a2 = X[2];
#endif
    X += DIM;

    const real* const end = mul + nbs;
    while ( mul < end )
    {
        real b0 = dif[0] * mul[0];
        real b1 = dif[1] * mul[0];
#if ( DIM > 2 )
        real b2 = dif[2] * mul[0];
#endif
        dif += DIM;
        Y[0] = a0 + b0;
        Y[1] = a1 + b1;
#if ( DIM > 2 )
        Y[2] = a2 + b2;
#endif
        Y += DIM;
        a0 = X[0] - b0;
        a1 = X[1] - b1;
#if ( DIM > 2 )
        a2 = X[2] - b2;
#endif
        X += DIM;
        ++mul;
    }
    
    Y[0] = a0;
    Y[1] = a1;
#if ( DIM > 2 )
    Y[2] = a2;
#endif
}

#if defined(__SSE3__) && ( DIM == 2 )

/**
 Perform second calculation needed by projectForces:
 */
void projectForcesD2D_SSE(SIZE_T nbs, const real* dif,
                          const real* src, const real* mul, real* dst)
{
    vec2 ss = load2(src);

    real const*const end = mul + nbs;
#if 0
    // unrolled
    while ( mul+3 < end )
    {
        vec2 d = mul2(load2(dif), loaddup2(mul));
        vec2 e = mul2(load2(dif+2), loaddup2(mul+1));
        dif += 4;
        mul += 2;
        store2(dst, add2(ss, d));
        store2(dst+2, add2(load2(src+2), sub2(e, d)));
        ss = sub2(load2(src+4), e);
        src += 4;
        dst += 4;
    }
#endif
    while ( mul < end )
    {
        vec2 d = mul2(load2(dif), loaddup2(mul));
        dif += 2;
        ++mul;
        store2(dst, add2(ss, d));
        ss = sub2(load2(src+2), d);
        src += 2;
        dst += 2;
    }
    store2(dst, ss);
}

#endif

#if defined(__AVX__) && ( DIM == 2 )


/**
 Perform second calculation needed by projectForces:
 Y <- X + Jt * tmp
 F. Nedelec, 9.12.2016
 */
void projectForcesD2D_AVX(SIZE_T nbs, const real* dif,
                          const real* src, const real* mul, real* dst)
{
    vec4 cc = setzero4();
    
    const bool odd = ( nbs & 1 );
    real const*const end = mul - odd + nbs;
    
    while ( mul < end )
    {
        vec4 t = broadcast2(mul);
        vec4 x = loadu4(src);
        mul += 2;
        vec4 m = permute4(t, 0b1100);
        vec4 d = mul4(m, load4(dif));
        dif += 4;
        vec4 n = twine2f128(cc,d);
        cc = d;
        vec4 z = add4(x, sub4(d, n));
        src += 4;
        storeu4(dst, z);
        dst += 4;
    }
    
    vec2 c = gethi(cc);
    
    if ( odd )
    {
        vec2 m = loaddup2(mul);
        vec2 x = mul2(m, load2(dif));
        dif += 2;
        vec2 z = add2(load2(src), sub2(x, c));
        storeu2(dst, z);
        c = x;
        dst += 2;
        src += 2;
    }
    
    vec2 z = sub2(load2(src), c);
    storeu2(dst, z);
}

#endif

/*
 3D workflow:
     Y[0] = X[0] + dif[0] * mul[0] - dif[-3] * mul[-1];
     Y[1] = X[1] + dif[1] * mul[0] - dif[-2] * mul[-1];

     Y[2] = X[2] + dif[2] * mul[0] - dif[-1] * mul[-1];
     Y[3] = X[3] + dif[3] * mul[1] - dif[0] * mul[0];

     Y[4] = X[4] + dif[4] * mul[1] - dif[1] * mul[0];
     Y[5] = X[5] + dif[5] * mul[1] - dif[2] * mul[0];
 */

#if REAL_IS_DOUBLE && defined(__SSE3__)
void projectForcesD3D_SSE(SIZE_T nbs, const real* dif,
                          const real* src, const real* mul, real* dst)
{
    real const*const end = mul - 1 + nbs;
    // first round involving 6 scalars = 3 SSE vectors
    if ( mul < end )
    {
        vec2 BC = loadu2(mul);
        vec2 BB = unpacklo2(BC, BC);
        vec2 CC = unpackhi2(BC, BC);
        vec2 AB = unpacklo2(setzero2(), BB);

        mul += 2;
        vec2 a0 = fmadd2(BB, load2(dif  ), loadu2(src  ));
        vec2 a1 = fmadd2(BC, load2(dif+2), loadu2(src+2));
        vec2 a2 = fmadd2(CC, load2(dif+4), loadu2(src+4));

        storeu2(dst  , a0);
        storeu2(dst+2, fnmadd2(AB, loaddup2(dif), a1));
        storeu2(dst+4, fnmadd2(BB, loadu2(dif+1), a2));
        dif += 6;
        dst += 6;
        src += 6;
    }
    else
    {
        // only 2 points overall
        vec2 BB = loaddup2(mul);
        vec2 AB = unpacklo2(setzero2(), BB);
        vec2 BC = unpacklo2(BB, setzero2());
        
        vec2 a0 = fmadd2(BB, load2(dif  ), loadu2(src  ));
        vec2 a1 = fmadd2(BC, load1(dif+2), loadu2(src+2));
        
        storeu2(dst  , a0);
        storeu2(dst+2, fnmadd2(AB, loaddup2(dif), a1));
        storeu2(dst+4, fnmadd2(BB, loadu2(dif+1), loadu2(src+4)));
        return;
    }
#if 0
    // unrolled processing 4 multipliers
    while ( mul < end-3 )
    {
        vec2 AB = loadu2(mul-1);
        vec2 CC = loaddup2(mul+1);
        vec2 DE = loadu2(mul+2);
        vec2 AA = unpacklo2(AB, AB);
        vec2 BB = unpackhi2(AB, AB);
        vec2 CD = unpacklo2(CC, DE);
        vec2 DD = unpacklo2(DE, DE);
        vec2 BC = unpackhi2(AB, CC);
        vec2 EE = unpackhi2(DE, DE);

        mul += 4;
        vec2 a0 = fnmadd2(loadu2(dif-3), AA, loadu2(src  ));
        vec2 a1 = fnmadd2(loadu2(dif-1), AB, loadu2(src+2));
        vec2 a2 = fnmadd2(loadu2(dif+1), BB, loadu2(src+4));
        
        vec2 b0 = fnmadd2(loadu2(dif+3), CC, loadu2(src+6));
        vec2 b1 = fnmadd2(loadu2(dif+5), CD, loadu2(src+8));
        vec2 b2 = fnmadd2(loadu2(dif+7), DD, loadu2(src+10));

        storeu2(dst  , fmadd2(load2(dif  ), BB, a0));
        storeu2(dst+2, fmadd2(load2(dif+2), BC, a1));
        storeu2(dst+4, fmadd2(load2(dif+4), CC, a2));
        
        storeu2(dst+ 6, fmadd2(load2(dif+6), DD, b0));
        storeu2(dst+ 8, fmadd2(load2(dif+8), DE, b1));
        storeu2(dst+10, fmadd2(load2(dif+10), EE, b2));
        dif += 12;
        dst += 12;
        src += 12;
    }
#endif
    // bulk work could be unrolled
    while ( mul < end )
    {
        vec2 AB = loadu2(mul-1);
        vec2 BC = loadu2(mul);
        vec2 AA = unpacklo2(AB, AB);
        vec2 BB = unpackhi2(AB, AB);
        vec2 CC = unpackhi2(BC, BC);

        mul += 2;
        vec2 a0 = fnmadd2(AA, loadu2(dif-3), loadu2(src  ));
        vec2 a1 = fnmadd2(AB, loadu2(dif-1), loadu2(src+2));
        vec2 a2 = fnmadd2(BB, loadu2(dif+1), loadu2(src+4));

        storeu2(dst  , fmadd2(BB, load2(dif  ), a0));
        storeu2(dst+2, fmadd2(BC, load2(dif+2), a1));
        storeu2(dst+4, fmadd2(CC, load2(dif+4), a2));
        dif += 6;
        dst += 6;
        src += 6;
    }
    if ( nbs & 1 )
    {
        // 2 points remaining = 6 scalars
        vec2 AA = loaddup2(mul-1);
        vec2 BB = loaddup2(mul);
        vec2 AB = unpacklo2(AA, BB);

        vec2 a0 = fnmadd2(AA, loadu2(dif-3), loadu2(src  ));
        vec2 a1 = fnmadd2(AB, loadu2(dif-1), loadu2(src+2));
        vec2 a2 = fnmadd2(BB, loadu2(dif+1), loadu2(src+4));

        storeu2(dst  , fmadd2(BB, load2(dif  ), a0));
        storeu2(dst+2, fmadd2(BB, load1(dif+2), a1));
        storeu2(dst+4, a2);
    }
    else
    {
        // 1 point remains = 3 scalars
        vec2 AA = loaddup2(mul-1);
        storeu2(dst, fnmadd2(AA, loadu2(dif-3), loadu2(src)));
        store1(dst+2, fnmadd1(AA, load1(dif-1), load1(src+2)));
    }
}
#endif

/*
 3D workflow:
 
     Y[0] = X[0] + dif[0] * mul[0] - dif[-3] * mul[-1];
     Y[1] = X[1] + dif[1] * mul[0] - dif[-2] * mul[-1];
     Y[2] = X[2] + dif[2] * mul[0] - dif[-1] * mul[-1];

     Y[3] = X[3] + dif[3] * mul[1] - dif[0] * mul[0];
     Y[4] = X[4] + dif[4] * mul[1] - dif[1] * mul[0];
     Y[5] = X[5] + dif[5] * mul[1] - dif[2] * mul[0];

     Y[6] = X[6] + dif[6] * mul[2] - dif[3] * mul[1];
     Y[7] = X[7] + dif[7] * mul[2] - dif[4] * mul[1];
     Y[8] = X[8] + dif[8] * mul[2] - dif[5] * mul[1];

     Y[9] = X[9] + dif[9] * mul[3] - dif[6] * mul[2];
     Y[A] = X[A] + dif[A] * mul[3] - dif[7] * mul[2];
     Y[B] = X[B] + dif[B] * mul[3] - dif[8] * mul[2];
 */

#if REAL_IS_DOUBLE && defined(__AVX__)

/*
 Ugly piece of code to harvest AVX power...
 FJN 18 and 19.04.2020
 */
void projectForcesD3D_AVX(SIZE_T nbs, const double* dir, const double* src, const double* mul, double* dst)
{
    const double* const end = mul - 3 + nbs;
    /*
     This follows the standard pattern defined below, except
     that the negative terms are not present on the first vector.
     This handles 12 scalars (4 vectors) in one round.
     */
    if ( mul < end )
    {
        vec4 m1 = broadcast1(mul  );
        vec4 m0 = blend31(setzero4(), m1);
        vec4 m2 = broadcast1(mul+1);
        vec4 p0 = blend31(m1, m2);
        vec4 p2 = broadcast1(mul+2);
        m1 = blend22(m1, m2);
        vec4 p1 = blend22(m2, p2);
        m2 = blend13(m2, p2);
        p2 = blend13(p2, broadcast1(mul+3));
        
        mul += 4;
        vec4 a0 = fmadd4(p0, load4(dir  ), loadu4(src  ));
        vec4 a1 = fmadd4(p1, load4(dir+4), loadu4(src+4));
        vec4 a2 = fmadd4(p2, load4(dir+8), loadu4(src+8));

        storeu4(dst  , fnmadd4(m0, broadcast1(dir), a0));
        storeu4(dst+4, fnmadd4(m1, loadu4(dir+1), a1));
        storeu4(dst+8, fnmadd4(m2, loadu4(dir+5), a2));
        dir += 12;
        dst += 12;
        src += 12;
    }
    /*
     This is where the bulk of the work is done, handing 12 scalars per pass.
     The loop can be unrolled.
     */
    while ( mul < end )
    {
        vec4 m0 = broadcast1(mul-1);
        vec4 m1 = broadcast1(mul  );
        m0 = blend31(m0, m1);
        vec4 m2 = broadcast1(mul+1);
        vec4 p0 = blend31(m1, m2);
        vec4 p2 = broadcast1(mul+2);
        m1 = blend22(m1, m2);
        vec4 p1 = blend22(m2, p2);
        m2 = blend13(m2, p2);
        p2 = blend13(p2, broadcast1(mul+3));
        
        mul += 4;
        vec4 a0 = fnmadd4(m0, loadu4(dir-3), loadu4(src  ));
        vec4 a1 = fnmadd4(m1, loadu4(dir+1), loadu4(src+4));
        vec4 a2 = fnmadd4(m2, loadu4(dir+5), loadu4(src+8));

        storeu4(dst  , fmadd4(p0, load4(dir  ), a0));
        storeu4(dst+4, fmadd4(p1, load4(dir+4), a1));
        storeu4(dst+8, fmadd4(p2, load4(dir+8), a2));

        dir += 12;
        dst += 12;
        src += 12;
    }
    /*
     We need to consider here multiple cases depending on how many vectors are
     left, since the positive terms are not present on the last vector.
     In addition the above code was not executed if there was less than 5 vectors
     in total, and the first vector is still a special case.
     We introduce: `mm` and `dd` for this reason, avoiding to load outside the
     valid data range.
     */
    vec4 mm, dd;
    if ( nbs < 4 )
    {
        mm = setzero4();
        dd = blend31(mm, broadcast1(dir));
    }
    else
    {
        mm = broadcast1(mul-1);
        dd = loadu4(dir-3);
    }
    switch ( mul - end )
    {
        case 0: {
            // 4 vectors remaining
            vec4 m1 = broadcast1(mul  );
            vec4 m0 = blend31(mm, m1);
            vec4 m2 = broadcast1(mul+1);
            vec4 p0 = blend31(m1, m2);
            vec4 p2 = broadcast1(mul+2);
            m1 = blend22(m1, m2);
            vec4 p1 = blend22(m2, p2);
            m2 = blend13(m2, p2);
            p2 = blend13(p2, setzero4());
            
            mul += 3;
            vec4 a0 = fmadd4(p0, load4(dir  ), loadu4(src  ));
            vec4 a1 = fmadd4(p1, load4(dir+4), loadu4(src+4));
            vec4 a2 = fmadd4(p2, broadcast1(dir+8), loadu4(src+8));
            a0 = fnmadd4(m0, dd, a0);
            a1 = fnmadd4(m1, loadu4(dir+1), a1);
            a2 = fnmadd4(m2, loadu4(dir+5), a2);
            storeu4(dst  , a0);
            storeu4(dst+4, a1);
            storeu4(dst+8, a2);
            //dir += 12; dst += 12; src += 12;
        } break;
        case 1: {
            // 3 vectors remaining
            vec4 m1 = broadcast1(mul  );
            vec4 m0 = blend31(mm, m1);
            vec4 m2 = broadcast1(mul+1);
            vec4 p0 = blend31(m1, m2);
            m1 = blend22(m1, m2);
            vec4 p1 = blend22(m2, setzero4());
            
            mul += 2;
            vec4 a0 = fmadd4(p0, load4(dir  ), loadu4(src  ));
            vec4 a1 = fmadd4(p1, load2Z(dir+4), loadu4(src+4));
            
            storeu4(dst  , fnmadd4(m0, dd, a0));
            storeu4(dst+4, fnmadd4(m1, loadu4(dir+1), a1));
             store1(dst+8, fnmadd2(getlo(m2), load1(dir+5), load1(src+8)));
            //dir += 9; dst += 9; src += 9;
        } break;
        case 2: {
            // 2 vectors remaining
            vec4 m1 = broadcast1(mul);
            vec4 m0 = blend31(mm, m1);
            vec4 p0 = blend31(m1, setzero4());
            
            mul += 1;
            vec4 a0 = fmadd4(p0, load3(dir), loadu4(src));
            
            storeu4(dst  , fnmadd4(m0, dd, a0));
            storeu2(dst+4, fnmadd2(getlo(m1), loadu2(dir+1), loadu2(src+4)));
            //dir += 6; dst += 6; src += 6;
        } break;
        case 3: {
            // 1 vector remaining
            //store3(dst, fnmadd4(mm, blend31(dd, setzero4()), load3(src)));
            storeu2(dst, fnmadd2(getlo(mm), getlo(dd), loadu2(src)));
            store1(dst+2, fnmadd1(getlo(mm), gethi(dd), load1(src+2)));
            //dir += 3; dst += 3; src += 3;
        } break;
        default:
            puts("unexpected case in projectForcesD3D_AVX!");
    }
    assert_true( mul == end+3 );
}

#endif

template < void (*FUNC)(SIZE_T, const real*, const real*, const real*, real*) >
void testD(SIZE_T cnt, char const* str)
{
    real *x = nullptr, *y = nullptr, *z = nullptr;
    new_nans(ALOC, x, y, z);
    randomize(NVAL, x, y, z, 1.0);
    nan_fill(NVAL, x);
    
    FUNC(NSEG, dir_, pos_, lag_, x);
    VecPrint::print(std::cout, std::min(DISP,NVAL), x);
    std::cout << " |";
    VecPrint::print(std::cout, DIM, x+NVAL);

    tic();
    for ( SIZE_T i=0; i<cnt; ++i )
    {
        FUNC(NSEG, dir_, x, lag_, z);
        // check the code with unaligned memory in 3D:
        FUNC(NSEG, dir_, y, lag_, x+DIM);
        FUNC(NSEG, dir_, z+DIM, lag_, y);
    }
    printf("  %10s %5.2f\n", str, toc(cnt*NSEG));
    
    free_reals(x,y,z);
}


void testProjectionD(SIZE_T cnt)
{
    projectForcesU_(NSEG, dir_, force_, lag_);
    std::cout << "testProjection DOWN " << DIM << "D,  " << NSEG << " segments\n";
    testD<projectForcesD_>(cnt,    " D_   ");
    testD<projectForcesD_ADD>(cnt, " D_ADD");
    testD<projectForcesD_REV>(cnt, " D_REV");
    testD<projectForcesD_RIV>(cnt, " D_RIV");
    testD<projectForcesD_FMA>(cnt, " D_FMA");
    testD<projectForcesD_PTR>(cnt, " D_PTR");
#if ( DIM == 2 ) && REAL_IS_DOUBLE && defined(__SSE__)
    testD<projectForcesD2D_SSE>(cnt, " D_SSE");
#endif
#if ( DIM == 2 ) && REAL_IS_DOUBLE && defined(__AVX__)
    testD<projectForcesD2D_AVX>(cnt, " D_AVX");
#endif
#if ( DIM == 3 ) && REAL_IS_DOUBLE && defined(__SSE__)
    testD<projectForcesD3D_SSE>(cnt, " D_SSE");
#endif
#if ( DIM == 3 ) && REAL_IS_DOUBLE && defined(__AVX__)
    testD<projectForcesD3D_AVX>(cnt, " D_AVX");
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Comparison against reference implementation


void projectForces(SIZE_T nbs, const real* X, real* Y)
{
    projectForcesU_(nbs, dir_, X, lag_);
    DPTTS2(nbs);
    projectForcesD_(nbs, dir_, X, lag_, Y);
}


#if REAL_IS_DOUBLE && defined(__SSE3__)
void projectForces_SSE(SIZE_T nbs, const real* X, real* Y)
{
#if ( DIM == 2 )
    projectForcesU2D_SSE(nbs, dir_, X, lag_);
#else
    projectForcesU3D_SSE(nbs, dir_, X, lag_);
#endif
    DPTTS2(nbs);
#if ( DIM == 2 )
    projectForcesD2D_SSE(nbs, dir_, X, lag_, Y);
#else
    projectForcesD3D_SSE(nbs, dir_, X, lag_, Y);
#endif
}
#endif


#if REAL_IS_DOUBLE && defined(__AVX__)
void projectForces_AVX(SIZE_T nbs, const real* X, real* Y)
{
#if ( DIM == 2 )
    projectForcesU2D_AVX(nbs, dir_, X, lag_);
#else
    projectForcesU3D_AVX(nbs, dir_, X, lag_);
#endif
    DPTTS2(nbs);
#if ( DIM == 2 )
    projectForcesD2D_AVX(nbs, dir_, X, lag_, Y);
#else
    projectForcesD3D_AVX(nbs, dir_, X, lag_, Y);
#endif
}
#endif


void checkProject()
{
    real *x = nullptr, *y = nullptr, *z = nullptr;

    for ( SIZE_T nbs = std::min(NSEG,(SIZE_T)11); nbs > 0; --nbs )
    {
        SIZE_T nbv = DIM * ( nbs + 1 );
        new_nans(nbv, x, y, z);
        setFilament(nbs, 0.1, 20.0, 1.0);
        setProjection(nbs);
        
        for ( SIZE_T i = 0; i < nbv; ++i )
            x[i] = RNG.sreal();
        //VecPrint::print(std::cout, nbv, dir_, 2); printf(" dir_\n");

        nan_fill(nbs, lag_);
        //projectForces(nbs, x, y);
        projectForcesU_(nbs, dir_, x, lag_);
        DPTTS2(nbs);
        projectForcesD_(nbs, dir_, x, lag_, y);
        printf("%2lu ", nbs);
        //VecPrint::print(std::cout, std::min(DISP,nbv), y);
        print_fe_exceptions("projectForces");
        printf(" SCA ");
        //printf("%2lu ", nbs); VecPrint::print(std::cout, std::min(DISP,nbs), lag_); printf(" lag_\n");

        nan_fill(nbs, lag_);
#if REAL_IS_DOUBLE && defined(__AVX__)
        projectForces_AVX(nbs, x, z);
        printf(" AVX ");
        //printf("%2lu ", nbs); VecPrint::print(std::cout, std::min(DISP,nbs), lag_); printf(" lag_\n");
        //printf("%2lu ", nbs); VecPrint::print(std::cout, std::min(DISP,nbv), z);
        print_fe_exceptions("projectForcesAVX");
#elif REAL_IS_DOUBLE && defined(__SSE3__)
        projectForces_SSE(nbs, x, z);
        printf(" SSE ");
        //printf("%2lu ", nbs); VecPrint::print(std::cout, std::min(DISP,nbs), lag_); printf(" lag_\n");
        //printf("%2lu ", nbs); VecPrint::print(std::cout, std::min(DISP,nbv), z);
        print_fe_exceptions("projectForces_SSE");
#else
        projectForces(nbs, x, z);
        printf(" ... ");
#endif
        real err = blas::max_diff(nbv, y, z);
        if ( abs_real(err) > 64*REAL_EPSILON )
            printf(" XXXX %e\n", err);
        else
            printf("  --> %e\n", err);
    }
    free_reals(x,y,z);
}


//------------------------------------------------------------------------------
#pragma mark - Fiber::projectForces()


void projectDPTTS(SIZE_T nbs, const real* X, real* Y)
{
    alsatian_xptts2(nbs, 1, diag_, upper_, lag_, NSEG);
}

void scaleTangentially(SIZE_T nbp, const real* X, const real* dir, real* Y)
{
    for ( SIZE_T p = 0; p < nbp; ++p )
    {
        real const* xxx = X   + DIM * p;
        real const* ddd = dir + DIM * p;
        real      * yyy = Y   + DIM * p;
#if ( DIM == 2 )
        real s = xxx[0] * ddd[0] + xxx[1] * ddd[1];
        yyy[0] = xxx[0] + s * ddd[0];
        yyy[1] = xxx[1] + s * ddd[1];
#elif ( DIM >= 3 )
        real s = xxx[0] * ddd[0] + xxx[1] * ddd[1] + xxx[2] * ddd[2];
        yyy[0] = xxx[0] + s * ddd[0];
        yyy[1] = xxx[1] + s * ddd[1];
        yyy[2] = xxx[2] + s * ddd[2];
#endif
    }
}

void scaleTangentiallyPTR(SIZE_T nbp, const real* src, const real* dir, real* dst)
{
    const real* const end = src + DIM * nbp;
    while ( src < end )
    {
#if ( DIM == 2 )
        real s = src[0] * dir[0] + src[1] * dir[1];
        dst[0] = src[0] + s * dir[0];
        dst[1] = src[1] + s * dir[1];
#elif ( DIM >= 3 )
        real s = src[0] * dir[0] + src[1] * dir[1] + src[2] * dir[2];
        dst[0] = src[0] + s * dir[0];
        dst[1] = src[1] + s * dir[1];
        dst[2] = src[2] + s * dir[2];
#endif
        src += DIM;
        dir += DIM;
        dst += DIM;
    }
}

void projectTangent(SIZE_T nbs, const real* X, real* Y)
{
    scaleTangentially(nbs+1, X, ani_, tmp_);
    projectForcesU_(nbs, dir_, tmp_, lag_);
    
    DPTTS2(nbs);

    projectForcesD_(nbs, dir_, X, lag_, Y);
    scaleTangentially(nbs+1, Y, ani_, Y);
}

void projectScale(SIZE_T nbs, const real* X, real* Y)
{
    scaleTangentially(nbs+1, X, ani_, Y);
    scaleTangentially(nbs+1, Y, ani_, Y);
}

template < void (*FUNC)(SIZE_T, const real*, real*) >
void timeProject(SIZE_T cnt, char const* str)
{
    real *x = nullptr, *y = nullptr, *z = nullptr;
    new_nans(ALOC, x, y, z);
    randomize(NVAL, x, y, z, 1.0);

    zero_real(ALOC, x);
    FUNC(NSEG, force_, x);
    VecPrint::print(std::cout, std::min(DISP,NVAL+2), x);

    tic();
    for ( SIZE_T ii=0; ii<cnt; ++ii )
    {
        FUNC(NSEG, x, y);
        // check the code with unaligned memory:
        FUNC(NSEG, y+2, z);
        FUNC(NSEG, z, x);
    }
    printf("  %10s %5.2f\n", str, toc(cnt*NSEG));
    free_reals(x,y,z);
}

//------------------------------------------------------------------------------
#pragma mark - Main

int main(int argc, char* argv[])
{
    const SIZE_T CNT = 1<<20;
    RNG.seed();
    checkProject();
    
    std::cout << __VERSION__ << "\n";
    if ( 1 )
    {
        setFilament(NSEG, 0.1, 20.0, 1.0);
        testProjectionU(CNT);
        testProjectionD(CNT);
    }
    if ( 1 )
    {
        setProjection(NSEG);
        setAnisotropy(NSEG);
        std::cout << "testProject " << DIM << "D, " << NSEG << " segments\n";
        timeProject<projectForces>(CNT,    " projF");
#if REAL_IS_DOUBLE && defined(__SSE3__)
        timeProject<projectForces_SSE>(CNT, " prSSE");
#endif
#if REAL_IS_DOUBLE && defined(__AVX__)
        timeProject<projectForces_AVX>(CNT, " prAVX");
#endif
        timeProject<projectDPTTS>(CNT,   " dptts");
        timeProject<projectTangent>(CNT, " projT");
        timeProject<projectScale>(CNT,   " scale");
    }
    
    free_reals(tmp_, lag_, force_, pos_);
    free_reals(dir_, ani_, diag_, upper_);
}
