// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <sys/time.h>

#define DIM 3

#include "assert_macro.h"
#include "real.h"
#include "vector.h"
#include "tictoc.h"
#include "random.h"
#include "vecprint.h"
#include "blas.h"
#include "lapack.h"
#include "dpttrf.h"
#include "cytoblas.h"
#include "simd.h"
#include "simd_print.h"


/// number of segments:
const size_t NBS = 6;
const size_t NCO = DIM * ( NBS + 1 );
const size_t ALOC = NCO + 8;

const size_t DISP = 24UL;


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

void setPoint(size_t n, Vector pos, Vector dir)
{
    pos.store(pos_+DIM*n);
    dir.store(dir_+DIM*n);
}

void setFilament(size_t np, real seg, real persistence_length)
{
    np = std::min(np, NBS+1);
    delete(pos_);
    delete(dir_);
    pos_ = new_real(NCO);
    dir_ = new_real(NCO);

    real sigma = sqrt(2.0*seg/persistence_length);
    
    Vector pos(0,0,0);
    Vector dir = Vector::randU();
    
    setPoint(0, pos, dir);
    for ( size_t p = 1 ; p < np; ++p )
    {
        pos += seg * dir;
        setPoint(p, pos, dir);
        //rotate dir in a random direction:
        real a = sigma * RNG.gauss();
        dir = cos(a) * dir + dir.randOrthoU(sin(a));
    }
}

void setProjection()
{
    delete(diag_);
    delete(upper_);
    diag_ = new_real(NBS);
    upper_ = new_real(NBS);
    
    size_t j = 0;
    for ( ; j < NBS-1; ++j )
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

    int info = 0;
    alsatian_xpttrf(NBS, diag_, upper_, &info);
}

void setAnisotropy()
{
    delete(ani_);
    ani_ = new_real(NCO);

    const size_t end = DIM * NBS;

    // for the extremities, the direction of the nearby segment is used.
    for ( size_t d = 0; d < DIM; ++d )
    {
        ani_[d]     = dir_[d];
        ani_[d+end] = dir_[d+end-DIM];
    }
    
    // for intermediate points, the directions of the flanking segments are averaged
    for ( size_t p = DIM ; p < end; ++p )
        ani_[p] = 0.5 * ( dir_[p-DIM] + dir_[p] );
}

void setRandom(int np, real * vec, real mag)
{
    for ( size_t p = 0 ; p < ALOC; ++p )
        vec[p] = mag * RNG.sreal();
}

void new_reals(real*& x, real*& y, real*& z, real mag)
{
    x = new_real(ALOC);
    y = new_real(ALOC);
    z = new_real(ALOC);
    
    for ( size_t ii=0; ii<ALOC; ++ii )
    {
        x[ii] = mag * RNG.sreal();
        y[ii] = mag * RNG.sreal();
        z[ii] = mag * RNG.sreal();
    }
}

void free_reals(real* x, real* y, real* z)
{
    free_real(x);
    free_real(y);
    free_real(z);
}

//------------------------------------------------------------------------------
#pragma mark - PROJECT UP

void projectForcesU_(size_t nbs, const real* dif, const real* src, real* mul)
{
    #pragma vector unaligned
    for ( size_t jj = 0; jj < nbs; ++jj )
    {
        const real * X = src + DIM * jj;
        const real * d = dif + DIM * jj;
        mul[jj] = d[0] * ( X[DIM  ] - X[0] )
                + d[1] * ( X[DIM+1] - X[1] )
#if ( DIM > 2 )
                + d[2] * ( X[DIM+2] - X[2] )
#endif
        ;
    }
}

void projectForcesU_PTR(size_t nbs, const real* dif, const real* src, real* mul)
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
inline void projectForcesU_TWO(size_t nbs, const real* dif, const real* src, real* mul)
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


#if REAL_IS_DOUBLE && defined __SSE3__

/**
 Perform first calculation needed by projectForces:
 tmp <- J * X
 */
void projectForcesU2D_SSE(size_t nbs, const real* dif, const real* src, real* mul)
{
    real const*const end = mul + nbs - 1;

    vec2 y, x = load2(src);
    while ( mul < end )
    {
        y = load2(src+2);
        src += 4;
        vec2 a = mul2(sub2(y, x), load2(dif));
        x = load2(src);
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
        //storelo(mul, hadd2(a, a));
        storelo(mul, add2(a, unpackhi2(a, a)));
    }
}

#endif

#ifdef __AVX__

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
void projectForcesU2D_AVX(size_t nbs, const real* dif, const real* src, real* mul)
{
    real const*const end = mul + nbs - 3;
    
#if ( 0 )
    while ( mul+4 < end )
    {
        vec4 a = mul4(sub4(loadu4(src+2 ), loadu4(src   )), loadc4(dif   ));
        vec4 b = mul4(sub4(loadu4(src+6 ), loadu4(src+4 )), loadc4(dif+4 ));
        vec4 c = mul4(sub4(loadu4(src+10), loadu4(src+8 )), loadc4(dif+8 ));
        vec4 d = mul4(sub4(loadu4(src+14), loadu4(src+12)), loadc4(dif+12));
        src += 16;
        dif += 16;
        store4(mul  , hadd4(permute2f128(a,b,0x20), permute2f128(a,b,0x31)));
        store4(mul+4, hadd4(permute2f128(c,d,0x20), permute2f128(c,d,0x31)));
        mul += 8;
    }
#endif

    while ( mul < end )
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

    while ( mul < end+2 )
    {
        //mul[jj] = dif[0] * ( X[DIM] - X[0] ) + dif[1] * ( X[DIM+1] - X[1] )
        vec4 d = mul4(sub4(loadu4(src+2), loadu4(src)), load4(dif));
        src += 4;
        dif += 4;
        vec2 h = gethi(d);
        storeu2(mul, add2(unpacklo2(getlo(d),h), unpackhi2(getlo(d),h)));
        mul += 2;
    }
    
    if ( mul < end+3 )
    {
        vec2 a = mul2(sub2(load2(src+2), load2(src)), load2(dif));
        storelo(mul, add2(a, permute2(a,1)));
    }
}


/**
 Attention: this does not check the boundaries and will write beyond the
 nbs-th point of tmp, which should be allocated accordingly.
 F. Nedelec, 11.01.2018
 */
void projectForcesU2D_AVY(size_t nbs, const real* dif, const real* src, real* mul)
{
    real *const end = mul + nbs;
    
    // calculate the terms 4 by 4
    while ( mul < end )
    {
        vec4 a = mul4(sub4(loadu4(src+2), loadu4(src  )), loadc4(dif  ));
        vec4 b = mul4(sub4(loadu4(src+6), loadu4(src+4)), loadc4(dif+4));
        dif += 8;
        src += 8;
        //store4(mul, hadd4(permute2f128(a,b,0x20), permute2f128(a,b,0x31)));
        vec4 p = permute2f128(a, b, 0x20);
        a = permute2f128(a, b, 0x31);
        store4(mul, add4(unpacklo4(p, a), unpackhi4(p, a)));
        mul += 4;
    }
}

#endif

//------------------------------------------------------------------------------
#pragma mark - PROJECT UP 3D

#if REAL_IS_DOUBLE && defined __AVX__

inline void twine3x4(real const* X, real const* Y, real const* Z, real* dst)
{
    vec4 sx = load4(X);
    vec4 sy = load4(Y);
    vec4 sz = load4(Z);

    vec4 zx = blend4(sx, sz, 0b0101);
    zx = permute2f128(zx, zx, 0x21);
    vec4 xy = unpacklo4(sx, sy);
    vec4 yz = unpackhi4(sy, sz);
    
    storeu4(dst  , blend4(xy, zx, 0b1100));
    storeu4(dst+4, blend4(yz, xy, 0b1100));
    storeu4(dst+8, blend4(zx, yz, 0b1100));
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

    vec4 zx = blend4(s0, s2, 0b0011);
    zx = permute2f128(zx, zx, 0x21);
    vec4 xy = blend4(s0, s1, 0b1100);
    vec4 yz = blend4(s1, s2, 0b1100);
    
    store4(X,   blend4(zx, xy, 0b0101));
    store4(Y, shuffle4(xy, yz, 0b0101));
    store4(Z,   blend4(zx, yz, 0b1010));
}


void projectForcesU3D_AVX(size_t nbs, const real* dif, const real* src, real* mul)
{
    const real *const end = mul + nbs - 3;
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

        vec4 zx = blend4(s0, s2, 0b0011);
        vec4 xy = blend4(s0, s1, 0b1100);
        zx = permute2f128(zx, zx, 0x21);
        vec4 yz = blend4(s1, s2, 0b1100);
        
        vec4 mm = shuffle4(xy, yz, 0b0101);
        mm = add4(mm, blend4(zx, xy, 0b0101));
        store4(mul, add4(mm, blend4(zx, yz, 0b1010)));
        
        src += 12;
        dif += 12;
        mul += 4;
    }
    while ( mul < end+2 )
    {
        /*
         *mul = dif[0] * ( src[DIM  ] - src[0] )
              + dif[1] * ( src[DIM+1] - src[1] )
              + dif[2] * ( src[DIM+2] - src[2] );
         */
        vec4 s0 = mul4(load4(dif  ), sub4(loadu4(src+3), loadu4(src)));
        vec4 s1 = mul4(load4(dif+4), sub4(broadcast1(src+7), loadu4(src+4)));

        vec4 xy = blend4(s0, s1, 0b1100);
        vec4 zx = permute2f128(s0, s0, 0x21);
        
        vec4 mm = shuffle4(xy, s1, 0b0101);
        mm = add4(mm, blend4(zx, xy, 0b0101));
        store2(mul, add4(mm, blend4(zx, s1, 0b1010)));
        
        src += 6;
        dif += 6;
        mul += 2;
    }
    while ( mul < end+3 )
    {
        /*
         *mul = dif[0] * ( src[DIM  ] - src[0] )
              + dif[1] * ( src[DIM+1] - src[1] )
              + dif[2] * ( src[DIM+2] - src[2] );
         */
        vec4 x = mul4(load3(dif), sub4(load3(src+3), loadu4(src)));
        vec4 z = permute2f128(x, x, 0x21);
        vec4 y = permute4(x, 0b0011);
        
        store1(mul, add4(add4(x, y), z));
        
        src += 3;
        dif += 3;
        ++mul;
    }
    assert_true(mul==end);
}

#endif

void testU(size_t cnt, void (*func)(size_t, const real*, const real*, real*), char const* str)
{
    real *x = nullptr, *y = nullptr, *z = nullptr;
    new_reals(x, y, z, 1.0);
    
    zero_real(ALOC, lag_);
    func(NBS, dir_, force_, lag_);
    VecPrint::print(std::cout, std::min(DISP,NBS+1), lag_);

    TicToc::tic();
    for ( size_t ii=0; ii<cnt; ++ii )
    {
        func(NBS, dir_, y, z);
        // check the code with unaligned memory:
        func(NBS, dir_, x+2, y+4);
        func(NBS, dir_, z+4, y+2);
    }
    TicToc::toc(str);
    
    free_reals(x,y,z);
}


void testProjectionU(size_t cnt)
{
    std::cout << "testProjection UP " << DIM << "D\n";
    testU(cnt, projectForcesU_,    " U_   ");
    testU(cnt, projectForcesU_PTR, " U_PTR");
    testU(cnt, projectForcesU_TWO, " U_TWO");
#if ( DIM == 2 ) && REAL_IS_DOUBLE && defined __SSE__
    testU(cnt, projectForcesU2D_SSE, " U_SSE");
#endif
#if ( DIM == 2 ) && REAL_IS_DOUBLE && defined __AVX__
    testU(cnt, projectForcesU2D_AVX, " U_AVX");
    testU(cnt, projectForcesU2D_AVY, " U_AVY");
#endif
#if ( DIM == 3 ) && REAL_IS_DOUBLE && defined __AVX__
    testU(cnt, projectForcesU3D_AVX, " U_AVX");
#endif
}


//------------------------------------------------------------------------------
#pragma mark - PROJECT DOWN

/**
 Perform second calculation needed by projectForces:
 Y <- X + Jt * tmp
 */
void projectForcesD_(size_t nbs, const real* dif, const real* X, const real* mul, real* Y)
{
    for ( size_t d = 0, e = DIM*nbs; d < DIM; ++d, ++e )
    {
        Y[d] = X[d] + dif[d    ] * mul[    0];
        Y[e] = X[e] - dif[e-DIM] * mul[nbs-1];
    }
    
    for ( size_t jj = 1; jj < nbs; ++jj )
    {
        const size_t kk = DIM*jj;
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
void projectForcesD_ADD(size_t nbs, const real* dif, const real* X, const real* mul, real* Y)
{
    real a0 = X[0];
    real a1 = X[1];
#if ( DIM > 2 )
    real a2 = X[2];
#endif
    
    for ( size_t jj = 0; jj < nbs; ++jj )
    {
        const size_t kk = DIM * jj;
        real b0 = dif[kk  ] * mul[jj];
        real b1 = dif[kk+1] * mul[jj];
#if ( DIM > 2 )
        real b2 = dif[kk+2] * mul[jj];
#endif
        Y[kk  ] = a0 + b0;
        Y[kk+1] = a1 + b1;
#if ( DIM > 2 )
        Y[kk+2] = a2 + b2;
#endif
        a0 = X[DIM+kk  ] - b0;
        a1 = X[DIM+kk+1] - b1;
#if ( DIM > 2 )
        a2 = X[DIM+kk+2] - b2;
#endif
    }
    
    const size_t ee = DIM * nbs;
    Y[ee  ] = a0;
    Y[ee+1] = a1;
#if ( DIM > 2 )
    Y[ee+2] = a2;
#endif
}



/**
 2xFMA
 */
void projectForcesD_FMA(size_t nbs, const real* dif, const real* X, const real* mul, real* Y)
{
    real a0 = X[0];
    real a1 = X[1];
#if ( DIM > 2 )
    real a2 = X[2];
#endif
    
    for ( size_t jj = 0; jj < nbs; ++jj )
    {
        const size_t kk = DIM * jj;
        real m = mul[jj];
        real d0 = dif[kk  ];
        real d1 = dif[kk+1];
#if ( DIM > 2 )
        real d2 = dif[kk+2];
#endif
        Y[kk  ] = a0 + d0 * m;
        Y[kk+1] = a1 + d1 * m;
#if ( DIM > 2 )
        Y[kk+2] = a2 + d2 * m;
#endif
        a0 = X[DIM+kk  ] - d0 * m;
        a1 = X[DIM+kk+1] - d1 * m;
#if ( DIM > 2 )
        a2 = X[DIM+kk+2] - d2 * m;
#endif
    }
    
    const size_t ee = DIM * nbs;
    Y[ee  ] = a0;
    Y[ee+1] = a1;
#if ( DIM > 2 )
    Y[ee+2] = a2;
#endif
}


/**
 Perform second calculation needed by projectForces:
 */
void projectForcesD_PTR(size_t nbs, const real* dif,
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

#if defined __SSE3__ && ( DIM == 2 )

/**
 Perform second calculation needed by projectForces:
 */
void projectForcesD2D_SSE(size_t nbs, const real* dif,
                          const real* src, const real* mul, real* dst)
{
    vec2 cc = load2(src);
    
    real const*const end = mul + nbs;
    while ( mul < end )
    {
        src += DIM;
        vec2 d = mul2(load2(dif), loaddup2(mul));
        ++mul;
        dif += DIM;
        store2(dst, add2(cc, d));
        dst += DIM;
        cc = sub2(load2(src), d);
    }
    store2(dst, cc);
}

#endif

#if defined __AVX__ && ( DIM == 2 )


/**
 Perform second calculation needed by projectForces:
 Y <- X + Jt * tmp
 F. Nedelec, 9.12.2016
 */
void projectForcesD2D_AVX(size_t nbs, const real* dif,
                          const real* src, const real* mul, real* dst)
{
    vec4 cc = setzero4();
    
    const bool odd = ( nbs & 1 );
    real const*const end = mul + nbs - odd;
    
    while ( mul < end )
    {
        vec4 t = broadcast2(mul);
        vec4 x = loadu4(src);
        mul += 2;
        vec4 m = permute4(t, 0b1100);
        vec4 d = mul4(m, load4(dif));
        dif += 4;
        vec4 n = permute2f128(cc,d,0x21);
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
 void projectForcesD_(size_t nbs, const real* dif, const real* X, const real* mul, real* Y)
 {
     for ( size_t d = 0, e = DIM*nbs; d < DIM; ++d, ++e )
     {
         Y[d] = X[d] + dif[d    ] * mul[    0];
         Y[e] = X[e] - dif[e-DIM] * mul[nbs-1];
     }
     
     for ( size_t jj = 1; jj < nbs; ++jj )
     {
         const size_t kk = DIM*jj;
         Y[kk  ] = X[kk  ] + dif[kk  ] * mul[jj] - dif[kk-DIM  ] * mul[jj-1];
         Y[kk+1] = X[kk+1] + dif[kk+1] * mul[jj] - dif[kk-DIM+1] * mul[jj-1];
         Y[kk+2] = X[kk+2] + dif[kk+2] * mul[jj] - dif[kk-DIM+2] * mul[jj-1];
     }
 
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

         Y[C] = X[C] + dif[C] * mul[4] - dif[9] * mul[3];
         Y[D] = X[D] + dif[D] * mul[4] - dif[A] * mul[3];
         Y[E] = X[E] + dif[E] * mul[4] - dif[B] * mul[3];
      }
 }
 */

#if REAL_IS_DOUBLE && defined __AVX__

/*
 Ugly piece of code to harvest AVX power...
 FJN 18 and 19.04.2020
 */
void projectForcesD3D_AVX(size_t nbs, const real* dif, const real* src, const real* mul, real* dst)
{
    const real* const end = mul + nbs - 3;
    /*
     This follows the standard pattern defined below, except
     that the negative terms are not present on the first vector.
     This handles 12 scalars (4 vectors) in one round.
     */
    if ( mul < end )
    {
        vec4 m1 = broadcast1(mul  );
        vec4 m0 = blend4(setzero4(), m1, 0b1000);
        vec4 m2 = broadcast1(mul+1);
        vec4 p0 = blend4(m1, m2, 0b1000);
        vec4 p2 = broadcast1(mul+2);
        m1 = blend4(m1, m2, 0b1100);
        vec4 p1 = blend4(m2, p2, 0b1100);
        m2 = blend4(m2, p2, 0b1110);
        p2 = blend4(p2, broadcast1(mul+3), 0b1110);
        
        mul += 4;
        vec4 a0 = fmadd4(p0, load4(dif  ), loadu4(src  ));
        vec4 a1 = fmadd4(p1, load4(dif+4), loadu4(src+4));
        vec4 a2 = fmadd4(p2, load4(dif+8), loadu4(src+8));

        storeu4(dst  , fnmadd4(m0, broadcast1(dif), a0));
        storeu4(dst+4, fnmadd4(m1, loadu4(dif+1), a1));
        storeu4(dst+8, fnmadd4(m2, loadu4(dif+5), a2));
        dif += 12;
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
        m0 = blend4(m0, m1, 0b1000);
        vec4 m2 = broadcast1(mul+1);
        vec4 p0 = blend4(m1, m2, 0b1000);
        vec4 p2 = broadcast1(mul+2);
        m1 = blend4(m1, m2, 0b1100);
        vec4 p1 = blend4(m2, p2, 0b1100);
        m2 = blend4(m2, p2, 0b1110);
        p2 = blend4(p2, broadcast1(mul+3), 0b1110);
        
        mul += 4;
        vec4 a0 = fmadd4(p0, load4(dif  ), loadu4(src  ));
        vec4 a1 = fmadd4(p1, load4(dif+4), loadu4(src+4));
        vec4 a2 = fmadd4(p2, load4(dif+8), loadu4(src+8));

        storeu4(dst  , fnmadd4(m0, loadu4(dif-3), a0));
        storeu4(dst+4, fnmadd4(m1, loadu4(dif+1), a1));
        storeu4(dst+8, fnmadd4(m2, loadu4(dif+5), a2));
        dif += 12;
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
        dd = broadcast1(dif);
    }
    else
    {
        mm = broadcast1(mul-1);
        dd = loadu4(dif-3);
    }
    switch ( mul - end )
    {
        case 0: {
            // 4 vectors remaining
            vec4 m1 = broadcast1(mul  );
            vec4 m0 = blend4(mm, m1, 0b1000);
            vec4 m2 = broadcast1(mul+1);
            vec4 p0 = blend4(m1, m2, 0b1000);
            vec4 p2 = broadcast1(mul+2);
            m1 = blend4(m1, m2, 0b1100);
            vec4 p1 = blend4(m2, p2, 0b1100);
            m2 = blend4(m2, p2, 0b1110);
            p2 = blend4(p2, setzero4(), 0b1110);
            
            mul += 4;
            vec4 a0 = fmadd4(p0, load4(dif  ), loadu4(src  ));
            vec4 a1 = fmadd4(p1, load4(dif+4), loadu4(src+4));
            vec4 a2 = fmadd4(p2, broadcast1(dif+8), loadu4(src+8));
            
            storeu4(dst  , fnmadd4(m0, dd, a0));
            storeu4(dst+4, fnmadd4(m1, loadu4(dif+1), a1));
            storeu4(dst+8, fnmadd4(m2, loadu4(dif+5), a2));
            dif += 12; dst += 12; src += 12;
        } break;
        case 1: {
            // 3 vectors remaining
            vec4 m1 = broadcast1(mul  );
            vec4 m0 = blend4(mm, m1, 0b1000);
            vec4 m2 = broadcast1(mul+1);
            vec4 p0 = blend4(m1, m2, 0b1000);
            m1 = blend4(m1, m2, 0b1100);
            vec4 p1 = blend4(m2, setzero4(), 0b1100);
            
            mul += 3;
            vec4 a0 = fmadd4(p0, load4(dif  ), loadu4(src  ));
            vec4 a1 = fmadd4(p1, load4(dif+4), loadu4(src+4));
            vec4 a2 = broadcast1(src+8);
            
            storeu4(dst  , fnmadd4(m0, dd, a0));
            storeu4(dst+4, fnmadd4(m1, loadu4(dif+1), a1));
            store1(dst+8, fnmadd4(m2, broadcast1(dif+5), a2));
            //storelo(dst+8, fnmadd2(getlo(m2), loaddup2(dif+5), getlo(a2)));
            dif += 9; dst += 9; src += 9;
        } break;
        case 2: {
            // 2 vectors remaining
            vec4 m1 = broadcast1(mul);
            vec4 m0 = blend4(mm, m1, 0b1000);
            vec4 p0 = blend4(m1, setzero4(), 0b1000);
            
            mul += 2;
            vec4 a0 = fmadd4(p0, load4(dif), loadu4(src));
            vec4 a1 = loadu4(src+4);
            
            storeu4(dst  , fnmadd4(m0, dd, a0));
            store2(dst+4, fnmadd4(m1, broadcast2(dif+1), a1));
            //store2(dst+4, fnmadd2(getlo(m1), loadu2(dif+1), getlo(a1)));
            dif += 6; dst += 6; src += 6;
        } break;
        case 3: {
            // 1 vector remaining
            ++mul;
            store3(dst, fnmadd4(mm, dd, load3(src)));
            dif += 3; dst += 3; src += 3;
        } break;
        default:
            printf("unexpected case in projectForcesD3D_AVX!");
    }
    assert_true( mul == end+4 );
}

#endif

void testD(size_t cnt, void (*func)(size_t, const real*, const real*, const real*, real*), char const* str)
{
    real *x = nullptr, *y = nullptr, *z = nullptr;
    new_reals(x, y, z, 1.0);
    
    zero_real(ALOC, x);
    func(NBS, dir_, pos_, lag_, x);
    VecPrint::print(std::cout, std::min(DISP,NCO+2), x);

    TicToc::tic();
    for ( size_t ii=0; ii<cnt; ++ii )
    {
        func(NBS, dir_, x, y, z);
        // check the code with unaligned memory:
        func(NBS, dir_, y+2, z, x+2);
        func(NBS, dir_, z+4, x, y+4);
    }
    TicToc::toc(str);
    
    free_reals(x,y,z);
}


void testProjectionD(size_t cnt)
{
    std::cout << "testProjection DOWN " << DIM << "D\n";
    testD(cnt, projectForcesD_,    " D_   ");
    testD(cnt, projectForcesD_ADD, " D_ADD");
    testD(cnt, projectForcesD_FMA, " D_FMA");
    testD(cnt, projectForcesD_PTR, " D_PTR");
#if ( DIM == 2 ) && REAL_IS_DOUBLE && defined __SSE__
    testD(cnt, projectForcesD2D_SSE, " D_SSE");
#endif
#if ( DIM == 2 ) && REAL_IS_DOUBLE && defined __AVX__
    testD(cnt, projectForcesD2D_AVX, " D_AVX");
#endif
#if ( DIM == 3 ) && REAL_IS_DOUBLE && defined __AVX__
    testD(cnt, projectForcesD3D_AVX, " D_AVX");
#endif
}


//------------------------------------------------------------------------------
#pragma mark - Fiber::projectForces()

void projectForces(size_t nbs, const real* X, real* Y)
{
#if defined __AVX__ && ( DIM == 2 )
    projectForcesU2D_AVX(nbs, dir_, X, lag_);
#else
    projectForcesU_(nbs, dir_, X, lag_);
#endif
    
    // find Lagrange multipliers
    alsatian_xptts2(nbs, 1, diag_, upper_, lag_, NBS);

#if defined __AVX__ && ( DIM == 2 )
    projectForcesD2D_AVX(nbs, dir_, X, lag_, Y);
#else
    projectForcesD_(nbs, dir_, X, lag_, Y);
#endif
}

void projectDPTTS(size_t nbs, const real* X, real* Y)
{
    alsatian_xptts2(nbs, 1, diag_, upper_, lag_, NBS);
}

void scaleTangentially(size_t nbp, const real* X, const real* dir, real* Y)
{
    for ( size_t p = 0; p < nbp; ++p )
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

void scaleTangentiallyPTR(size_t nbp, const real* src, const real* dir, real* dst)
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

void projectTangent(size_t nbs, const real* X, real* Y)
{
    scaleTangentially(nbs+1, X, ani_, tmp_);
    projectForcesU_(nbs, dir_, tmp_, lag_);
    
    // find Lagrange multipliers
    alsatian_xptts2(nbs, 1, diag_, upper_, lag_, NBS);

    projectForcesD_(nbs, dir_, X, lag_, Y);
    scaleTangentially(nbs+1, Y, ani_, Y);
}

void projectScale(size_t nbs, const real* X, real* Y)
{
    scaleTangentially(nbs+1, X, ani_, Y);
    scaleTangentially(nbs+1, Y, ani_, Y);
}

void testProject(size_t cnt, void (*func)(size_t, const real*, real*), char const* str)
{
    real *x = nullptr, *y = nullptr, *z = nullptr;
    new_reals(x, y, z, 1.0);
    
    zero_real(ALOC, x);
    func(NBS, force_, x);
    VecPrint::print(std::cout, std::min(DISP,NCO+2), x);

    TicToc::tic();
    for ( size_t ii=0; ii<cnt; ++ii )
    {
        func(NBS, x, y);
        // check the code with unaligned memory:
        func(NBS, y+2, z);
        func(NBS, z, x);
    }
    TicToc::toc(str);
    free_reals(x,y,z);
}

//------------------------------------------------------------------------------
#pragma mark - Main

int main(int argc, char* argv[])
{
    RNG.seed();
    new_reals(force_, lag_, tmp_, 1.0);
    setRandom(NBS+1, force_, 1.0);
    setFilament(NBS+1, 1.0, 2.0);
    
    const size_t CNT = 1<<20;
    
    if ( 1 )
    {
        testProjectionU(CNT);
        testProjectionD(CNT);
    }
    if ( 0 )
    {
        setProjection();
        setAnisotropy();
        std::cout << "testProject " << DIM << "D\n";
        testProject(CNT, projectForces,  " projF");
        testProject(CNT, projectDPTTS,   " dptts");
        testProject(CNT, projectTangent, " projT");
        testProject(CNT, projectScale,   " scale");
    }
    
    free_reals(tmp_, lag_, force_);
    free_reals(pos_, dir_, ani_);
    free_reals(diag_, upper_, nullptr);

    return EXIT_SUCCESS;
}
