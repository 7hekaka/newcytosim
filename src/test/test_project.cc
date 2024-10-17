// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include <sys/time.h>
#include <limits>
#include <cfenv>

#define DIM 3

#include "assert_macro.h"
#include "real.h"
#include "timer.h"
#include "vector.h"
#include "random.h"
#include "vecprint.h"
#include "blas.h"
#include "lapack.h"
#include "cytoblas.h"
#include "simd.h"
#include "simd_float.h"

#include "../math/fp_exceptions.cc"
#include "../sim/mecafil_code.cc"


/// number of segments:
const index_t NSEG = 117;
/// number of vertex coordinates
const index_t NVAL = DIM * ( NSEG + 1 );
const index_t ALOC = NVAL + 8;
const index_t SUB = 64;

/// vectors for one filament
real *pos_=nullptr, *dir_=nullptr, *ani_=nullptr;
real *lag_=nullptr, *force_=nullptr, *tmp_=nullptr;
real *diag_=nullptr, *upper_=nullptr;

//

/// fill vector with signalling NANs
void nan_fill(index_t cnt, real * ptr)
{
    real n = std::numeric_limits<real>::signaling_NaN();
    for ( index_t u = 0; u < cnt; ++u )
        ptr[u] = n;
}


void new_nans(index_t cnt, real*& ptr)
{
    free_real(ptr);
    ptr = new_real(cnt);
    nan_fill(cnt, ptr);
    //printf("%p = new_nans(%lu)\n", ptr, cnt);
}

void new_nans(index_t cnt, real*& x, real*& y, real*& z)
{
    new_nans(cnt, x);
    new_nans(cnt, y);
    new_nans(cnt, z);
}

void randomize(index_t cnt, real*& x, real*& y, real*& z, real mag)
{
    for ( index_t i = 0; i < cnt; ++i )
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


//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

/// reference implementation
void projectForcesU_(index_t nbs, const real* dir, const real* src, real* mul)
{
    #pragma omp simd
    for ( index_t i = 0; i < nbs; ++i )
    {
        const real * X = src + DIM * i;
        const real * d = dir + DIM * i;
        mul[i] = d[0] * ( X[DIM  ] - X[0] )
               + d[1] * ( X[DIM+1] - X[1] )
#if ( DIM > 2 )
               + d[2] * ( X[DIM+2] - X[2] )
#endif
        ;
    }
}

inline void DPTTRF(index_t nbs)
{
    int info = 0;
    lapack::xpttrf(nbs, diag_, upper_, &info);
}

inline void DPTTS2(index_t nbs, real* Y)
{
    lapack::xptts2(nbs, 1, diag_, upper_, Y, 0);
}

//------------------------------------------------------------------------------

void setFilament(index_t nbs, real seg, real persistence_length, real mag)
{
    nbs = std::min(nbs, NSEG);
    index_t nbv = DIM * ( nbs + 1 );
    new_nans(nbv, pos_);
    new_nans(nbv, dir_);
    new_nans(nbv, lag_);
    new_nans(nbv, force_);

    real sigma = std::sqrt(2.0*seg/persistence_length);
    
    Vector pos(0,0,0);
    Vector dir = Vector::randU();
    
    pos.store(pos_);
    dir.store(dir_);
    for ( index_t p = 1 ; p < nbs; ++p )
    {
        pos += seg * dir;
        pos.store(pos_+DIM*p);
        dir.store(dir_+DIM*p);
        //rotate dir in a random direction:
        real a = sigma * RNG.gauss();
        dir = std::cos(a) * dir + dir.randOrthoU(std::sin(a));
    }
    pos.store(pos_+DIM*nbs);
    print_floating_point_exceptions("setFilament");

    for ( index_t i = 0 ; i < nbv; ++i )
        force_[i] = mag * RNG.sreal();
    
    print_floating_point_exceptions("setForce");
}


void setProjection(index_t nbs)
{
    new_nans(nbs, diag_);
    new_nans(nbs, upper_);
    
    index_t j = 0;
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
    upper_[j] = 0.0;
    print_floating_point_exceptions("setProjection");

    DPTTRF(nbs);
    print_floating_point_exceptions("DPTTRF");
    
    projectForcesU_(nbs, dir_, force_, lag_);
    DPTTS2(nbs, lag_);
    print_floating_point_exceptions("DPTTS2");

#if ( 0 )
    VecPrint::edges("D", nbs, diag_, 3);
    VecPrint::edges("U", nbs-1, upper_, 3);
    VecPrint::edges("lag", nbs-1, lag_, 3);
#endif
}


void setAnisotropy(index_t nbs)
{
    const index_t sup = DIM * nbs;
    new_nans(sup+DIM, ani_);
    new_nans(sup+DIM, tmp_);

    // for the extremities, the direction of the nearby segment is used.
    for ( index_t d = 0; d < DIM; ++d )
    {
        ani_[d]     = dir_[d];
        ani_[d+sup] = dir_[d+sup-DIM];
    }
    // for intermediate points, the directions of the flanking segments are averaged
    for ( index_t p = DIM ; p < sup; ++p )
        ani_[p] = 0.5 * ( dir_[p-DIM] + dir_[p] );
    
    print_floating_point_exceptions("setAnisotropy");
}


void setRandom(int np, real * vec, real mag)
{
    for ( index_t p = 0 ; p < ALOC; ++p )
        vec[p] = mag * RNG.sreal();
}

//------------------------------------------------------------------------------
#pragma mark - PROJECT UP

/**
 Perform first calculation needed by projectForces:
 */
inline void projectForcesU_TWO(index_t nbs, const real* dir, const real* src, real* mul)
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
        mul[0] = dir[0] * (x3 - x0) + dir[1] * (x4 - x1);
#elif ( DIM >= 3 )
        x5 = src[2];
        mul[0] = dir[0] * (x3 - x0) + dir[1] * (x4 - x1) + dir[2] * (x5 - x2);
        x2 = x5;
#endif
        ++mul;
        src += DIM;
        dir += DIM;
        x0 = x3;
        x1 = x4;
    }
    
    while ( mul < end )
    {
        x3 = src[0];
        x4 = src[1];
#if ( DIM == 2 )
        mul[0] = dir[0] * (x3 - x0) + dir[1] * (x4 - x1);
#elif ( DIM >= 3 )
        x5 = src[2];
        mul[0] = dir[0] * (x3 - x0) + dir[1] * (x4 - x1) + dir[2] * (x5 - x2);
#endif
        
#if ( DIM == 2 )
        x0 = src[2];
        x1 = src[3];
        mul[1] = dir[2] * (x0 - x3) + dir[3] * (x1 - x4);
#elif ( DIM >= 3 )
        x0 = src[3];
        x1 = src[4];
        x2 = src[5];
        mul[1] = dir[3] * (x0 - x3) + dir[4] * (x1 - x4) + dir[5] * (x2 - x5);
#endif
        
        mul += 2;
        src += 2*DIM;
        dir += 2*DIM;
    }
    assert_true( mul == end );
}

#if REAL_IS_DOUBLE && defined(__AVX__)

/**
 Attention: this does not check the boundaries and will write beyond the
 nbs-th point of tmp, which should be allocated accordingly.
 F. Nedelec, 11.01.2018
 */
void projectForcesU2D_AVY(index_t nbs, const real* dir, const real* src, real* mul)
{
    real *const end = mul - 3 + nbs;
    
    // calculate the terms 4 by 4
    while ( mul < end )
    {
        vec4 a = mul4(sub4(loadu4(src+2), loadu4(src  )), load4(dir  ));
        vec4 b = mul4(sub4(loadu4(src+6), loadu4(src+4)), load4(dir+4));
        dir += 8;
        src += 8;
        vec4 p = unpacklo2f128(a, b);
        vec4 q = unpackhi2f128(a, b);
        store4(mul, add4(unpacklo4(p, q), unpackhi4(p, q)));
        mul += 4;
    }
}

#endif


//------------------------------------------------------------------------------
#pragma mark - PROJECT DOWN

/**
 Perform second calculation needed by projectForces:
 Y <- X + Jt * tmp
 */
void projectForcesD_(const index_t nbs, const real* dir, const real* X, const real* mul, real* Y)
{
    for ( index_t d = 0; d < DIM; ++d )
        Y[d] = X[d] + dir[d] * mul[0];

    for ( index_t e = DIM*nbs; e < DIM*(nbs+1); ++e )
        Y[e] = X[e] - dir[e-DIM] * mul[nbs-1];

    for ( index_t j = 1; j < nbs; ++j )
    {
        const index_t kk = DIM * j;
        Y[kk  ] = X[kk  ] + dir[kk  ] * mul[j] - dir[kk-DIM  ] * mul[j-1];
        Y[kk+1] = X[kk+1] + dir[kk+1] * mul[j] - dir[kk-DIM+1] * mul[j-1];
#if ( DIM > 2 )
        Y[kk+2] = X[kk+2] + dir[kk+2] * mul[j] - dir[kk-DIM+2] * mul[j-1];
#endif
    }
}


/**
 1xMUL 2xADD
 */
void projectForcesD_ADD(index_t nbs, const real* dir, const real* X, const real* mul, real* Y)
{
    real a0 = X[0];
    real a1 = X[1];
#if ( DIM > 2 )
    real a2 = X[2];
#endif
    
    for ( index_t jj = 0; jj < nbs; ++jj )
    {
        real b0 = dir[DIM*jj  ] * mul[jj];
        real b1 = dir[DIM*jj+1] * mul[jj];
#if ( DIM > 2 )
        real b2 = dir[DIM*jj+2] * mul[jj];
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
    
    const index_t ee = DIM * nbs;
    Y[ee  ] = a0;
    Y[ee+1] = a1;
#if ( DIM > 2 )
    Y[ee+2] = a2;
#endif
}


/// proceeding downward, this could be fused with downward part of DPTTS2
void projectForcesD_REV(index_t nbs, const real* dir, const real* X, const real* mul, real* Y)
{
    const real * E = X + DIM * nbs;
    real a0 = E[0];
    real a1 = E[1];
#if ( DIM > 2 )
    real a2 = E[2];
#endif
    
    for ( index_t i = nbs; i-- > 0; )
    {
        real m0 = dir[DIM*i  ] * mul[i];
        real m1 = dir[DIM*i+1] * mul[i];
#if ( DIM > 2 )
        real m2 = dir[DIM*i+2] * mul[i];
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
void projectForcesD_RIV(index_t nbs, const real* dir, const real* X, const real* mul, real* Y)
{
    const real * E = X + DIM * nbs;
    real a0 = E[0];
    real a1 = E[1];
#if ( DIM > 2 )
    real a2 = E[2];
#endif

    mul += nbs-1;
    dir += DIM*(nbs-1);
    Y += DIM*nbs;
    
    while ( E > X )
    {
        E -= DIM;
        real m0 = dir[0] * mul[0];
        real m1 = dir[1] * mul[0];
#if ( DIM > 2 )
        real m2 = dir[2] * mul[0];
#endif
        --mul;
        dir -= DIM;
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
void projectForcesD_FMA(index_t nbs, const real* dir, const real* X, const real* mul, real* Y)
{
    real a0 = X[0];
    real a1 = X[1];
#if ( DIM > 2 )
    real a2 = X[2];
#endif
    
    for ( index_t jj = 0; jj < nbs; ++jj )
    {
        real m = mul[jj];
        real d0 = dir[DIM*jj  ];
        real d1 = dir[DIM*jj+1];
#if ( DIM > 2 )
        real d2 = dir[DIM*jj+2];
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
    
    const index_t ee = DIM * nbs;
    Y[ee  ] = a0;
    Y[ee+1] = a1;
#if ( DIM > 2 )
    Y[ee+2] = a2;
#endif
}

/*
 3D workflow:
 
     Y[0] = X[0] + dir[0] * mul[0] - dir[-3] * mul[-1];
     Y[1] = X[1] + dir[1] * mul[0] - dir[-2] * mul[-1];
     Y[2] = X[2] + dir[2] * mul[0] - dir[-1] * mul[-1];

     Y[3] = X[3] + dir[3] * mul[1] - dir[0] * mul[0];
     Y[4] = X[4] + dir[4] * mul[1] - dir[1] * mul[0];
     Y[5] = X[5] + dir[5] * mul[1] - dir[2] * mul[0];

     Y[6] = X[6] + dir[6] * mul[2] - dir[3] * mul[1];
     Y[7] = X[7] + dir[7] * mul[2] - dir[4] * mul[1];
     Y[8] = X[8] + dir[8] * mul[2] - dir[5] * mul[1];

     Y[9] = X[9] + dir[9] * mul[3] - dir[6] * mul[2];
     Y[A] = X[A] + dir[A] * mul[3] - dir[7] * mul[2];
     Y[B] = X[B] + dir[B] * mul[3] - dir[8] * mul[2];
 */

#if ( DIM == 3 ) && USE_SIMD
// SSE version using FMA
void projectForcesD3D_sse(index_t nbs, const double* dir,
                          const double* src, const double* mul, double* dst)
{
    double const*const end = mul + nbs - 1;
    vec2 d0 = load2(dir);
    vec2 d1 = load2(dir+2);
    vec2 d2;

    // first round involving 6 scalars = 3 SSE vectors
    if ( mul < end )
    {
        vec2 BC = loadu2(mul);
        vec2 BB = unpacklo2(BC, BC);
        vec2 CC = unpackhi2(BC, BC);
        vec2 AB = unpacklo2(setzero2(), BB);

        d2 = load2(dir+4);

        mul += 2;
        vec2 a0 = fmadd2(BB, d0, loadu2(src  ));
        vec2 a1 = fmadd2(BC, d1, loadu2(src+2));
        vec2 a2 = fmadd2(CC, d2, loadu2(src+4));

        storeu2(dst  , a0);
        storeu2(dst+2, fnmadd2(AB, catshift(d0, d0), a1));
        storeu2(dst+4, fnmadd2(BB, catshift(d0, d1), a2));
        dir += 6;
        dst += 6;
        src += 6;
    }
    else
    {
        // only 2 points overall and one multiplicator
        d1 = unpacklo2(d1, setzero2());

        vec2 BB = loaddup2(mul);
        vec2 AB = unpacklo2(setzero2(), BB);
        vec2 BC = unpacklo2(BB, setzero2());

        vec2 a0 = fmadd2(BB, d0, loadu2(src  ));
        vec2 a1 = fmadd2(BC, d1, loadu2(src+2));
        
        storeu2(dst  , a0);
        storeu2(dst+2, fnmadd2(AB, catshift(d0, d0), a1));
        storeu2(dst+4, fnmadd2(BB, catshift(d0, d1), loadu2(src+4)));
        return;
    }
#if 1
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
        vec2 dA = d1; // dir-4
        vec2 dB = d2; // dir-2
        vec2 dC = load2(dir  );
        vec2 dD = load2(dir+2);
        vec2 dE = load2(dir+4);
        d0 = load2(dir+6);
        d1 = load2(dir+8);
        d2 = load2(dir+10);

        vec2 a0 = fnmadd2(catshift(dA, dB), AA, loadu2(src  ));
        vec2 a1 = fnmadd2(catshift(dB, dC), AB, loadu2(src+2));
        vec2 a2 = fnmadd2(catshift(dC, dD), BB, loadu2(src+4));
        vec2 b0 = fnmadd2(catshift(dD, dE), CC, loadu2(src+6));
        vec2 b1 = fnmadd2(catshift(dE, d0), CD, loadu2(src+8));
        vec2 b2 = fnmadd2(catshift(d0, d1), DD, loadu2(src+10));

        storeu2(dst  , fmadd2(dC, BB, a0));
        storeu2(dst+2, fmadd2(dD, BC, a1));
        storeu2(dst+4, fmadd2(dE, CC, a2));
        storeu2(dst+ 6, fmadd2(d0, DD, b0));
        storeu2(dst+ 8, fmadd2(d1, DE, b1));
        storeu2(dst+10, fmadd2(d2, EE, b2));
        
        dir += 12;
        dst += 12;
        src += 12;
    }
#endif
    // bulk work processing 6 scalars = 2 vectors
    while ( mul < end )
    {
        vec2 AB = loadu2(mul-1);
        vec2 BC = loadu2(mul);
        vec2 AA = unpacklo2(AB, AB);
        vec2 BB = unpackhi2(AB, AB);
        vec2 CC = unpackhi2(BC, BC);

        mul += 2;
        vec2 dA = d1;
        vec2 dB = d2;
        d0 = load2(dir  );
        d1 = load2(dir+2);
        d2 = load2(dir+4);
        
        vec2 a0 = fnmadd2(AA, catshift(dA, dB), loadu2(src  ));
        vec2 a1 = fnmadd2(AB, catshift(dB, d0), loadu2(src+2));
        vec2 a2 = fnmadd2(BB, catshift(d0, d1), loadu2(src+4));

        storeu2(dst  , fmadd2(BB, d0, a0));
        storeu2(dst+2, fmadd2(BC, d1, a1));
        storeu2(dst+4, fmadd2(CC, d2, a2));
        
        dir += 6;
        dst += 6;
        src += 6;
    }
    // end points
    if ( nbs & 1 )
    {
        // 2 points remaining = 6 scalars
        vec2 AA = loaddup2(mul-1);
        vec2 BB = loaddup2(mul);
        vec2 AB = unpacklo2(AA, BB);

        vec2 dA = d1;
        vec2 dB = d2;
        d0 = load2(dir);
        d1 = unpacklo2(load2(dir+2), setzero2());

        vec2 a0 = fnmadd2(AA, catshift(dA, dB), loadu2(src  ));
        vec2 a1 = fnmadd2(AB, catshift(dB, d0), loadu2(src+2));
        vec2 a2 = fnmadd2(BB, catshift(d0, d1), loadu2(src+4));

        storeu2(dst  , fmadd2(BB, d0, a0));
        storeu2(dst+2, fmadd2(BB, d1, a1));
        storeu2(dst+4, a2);
    }
    else
    {
        // 1 point remains = 3 scalars
        vec2 AA = loaddup2(mul-1);

        storeu2(dst, fnmadd2(AA, catshift(d1, d2), loadu2(src)));
        store1(dst+2, fnmadd1(AA, catshift(d2, d2), load1(src+2)));
    }
}
#endif

//------------------------------------------------------------------------------
void addProjectionDiff_(const index_t nbs, const real* mul, const real* X, real* Y)
{
    for ( index_t i = 0; i < nbs; ++i )
    {
        for ( index_t d = 0; d < DIM; ++d )
        {
            const real w = mul[i] * ( X[DIM*i+DIM+d] - X[DIM*i+d] );
            Y[DIM*i+(    d)] += w;
            Y[DIM*i+(DIM+d)] -= w;
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - Comparison against reference implementation


void projectForces(index_t nbs, const real* X, real* Y)
{
    projectForcesU_(nbs, dir_, X, lag_);
    DPTTS2(nbs, lag_);
    projectForcesD_(nbs, dir_, X, lag_, Y);
}


#if REAL_IS_DOUBLE && USE_SIMD
void projectForces_SSE(index_t nbs, const double* X, real* Y)
{
#if ( DIM == 2 )
    projectForcesU2D_SSE(nbs, dir_, X, lag_);
#else
    projectForcesU3D_SSE(nbs, dir_, X, lag_);
#endif
    DPTTS2(nbs, lag_);
#if ( DIM == 2 )
    projectForcesD2D_SSE(nbs, dir_, X, lag_, Y);
#else
    projectForcesD3D_SSE(nbs, dir_, X, lag_, Y);
#endif
}
#elif defined(__SSE3__)
void projectForces_SSE(index_t nbs, const float* X, real* Y)
{
    projectForcesU3D_SSE(nbs, dir_, X, lag_);
    DPTTS2(nbs, lag_);
    projectForcesD3D_SSE(nbs, dir_, X, lag_, Y);
}
#endif


#if REAL_IS_DOUBLE && defined(__AVX__)
void projectForces_AVX(index_t nbs, const double* X, real* Y)
{
#if ( DIM == 2 )
    projectForcesU2D_AVX(nbs, dir_, X, lag_);
#else
    projectForcesU3D_AVX(nbs, dir_, X, lag_);
#endif
    DPTTS2(nbs, lag_);
#if ( DIM == 2 )
    projectForcesD2D_AVX(nbs, dir_, X, lag_, Y);
#else
    projectForcesD3D_AVX(nbs, dir_, X, lag_, Y);
#endif
}
#endif


template < void (*FUNC)(index_t, const real*, real*) >
void checkProject(index_t nbs, const char msg[])
{
    printf("%s projectForces %2u ", msg, nbs);
    real *x = nullptr, *y = nullptr, *z = nullptr;
    index_t nbv = DIM * ( nbs + 1 );
    new_nans(nbv+4, x, y, z);
    setFilament(nbs, 0.1, 20.0, 1.0);
    setProjection(nbs);
    std::feclearexcept(FE_ALL_EXCEPT);
    
    for ( index_t i = 0; i < nbv; ++i )
        x[i] = RNG.sreal();
    
    nan_fill(nbs, lag_);
    projectForces(nbs, x, y);
    print_floating_point_exceptions("SCA");
    
    nan_fill(nbs, lag_);
    FUNC(nbs, x, z);
    print_floating_point_exceptions(msg);
    
    real err = blas::difference(nbv, y, z);
    if ( not ( abs_real(err) < 64*REAL_EPSILON ) )
    {
        printf(" WRONG! %e:\n", err);
        VecPrint::edges("ref ", nbv, y);
        VecPrint::edges("    ", nbv, z);
    }
    else
        printf(" okay--> %e\n", err);
    free_reals(x,y,z);
}


//------------------------------------------------------------------------------
#pragma mark - Adaptors

void onlyDPTTS(index_t nbs, const real* X, real* Y)
{
    lapack::xptts2(nbs, 1, diag_, upper_, Y, NSEG);
}

void projectDiff_(index_t nbs, const real* X, real* Y)
{
    addProjectionDiff_(nbs, lag_, X, Y);
}

void projectDiff_R(index_t nbs, const real* X, real* Y)
{
    addProjectionDiff_R(nbs, lag_, X, Y);
}

void projectDiff_F(index_t nbs, const real* X, real* Y)
{
    addProjectionDiff_F(nbs, lag_, X, Y);
}

#if USE_SIMD
void projectDiff_SSE(index_t nbs, const real* X, real* Y)
{
#if ( DIM == 2 )
    addProjectionDiff2D_SSE(nbs, lag_, X, Y);
#elif ( DIM == 3 ) && REAL_IS_DOUBLE
    addProjectionDiff3D_SSE(nbs, lag_, X, Y);
#endif
}
#endif
#if ( DIM == 2 ) && REAL_IS_DOUBLE && defined(__AVX__)
void projectDiff_AVX(index_t nbs, const real* X, real* Y)
{
    addProjectionDiff_AVX(nbs, lag_, X, Y);
}
#endif


//------------------------------------------------------------------------------
#pragma mark - Test UP & Down

template <void (*FUNC)(index_t, const real*, const real*, real*)>
void testU(index_t cnt, char const* str)
{
    real *x = nullptr, *y = nullptr, *z = nullptr;
    new_nans(ALOC, x, y, z);
    randomize(NVAL, x, y, z, 1.0);
    nan_fill(NVAL, lag_);

    FUNC(NSEG, dir_, force_, lag_);
    VecPrint::edges(NSEG, lag_);
    print_floating_point_exceptions("");
    
    tick();
    for ( index_t i=0; i<cnt; ++i )
    {
        randomize(NVAL, x, y, z, 1.0);
        for ( index_t j=0; j<SUB; ++j )
        {
            FUNC(NSEG, dir_, y, z);
            FUNC(NSEG, dir_, x, y);
            FUNC(NSEG, dir_, z, y);
        }
    }
    printf("  %8s %5.0f", str, tock());
    print_floating_point_exceptions("");
    printf("\n");
    free_reals(x,y,z);
}


template < void (*FUNC)(index_t, const real*, const real*, const real*, real*) >
void testD(index_t cnt, char const* str)
{
    real *x = nullptr, *y = nullptr, *z = nullptr;
    new_nans(ALOC, x, y, z);
    randomize(NVAL, x, y, z, 1.0);
    nan_fill(NVAL, x);
    
    FUNC(NSEG, dir_, pos_, lag_, x);
    VecPrint::edges(NVAL, x);
    print_floating_point_exceptions("");
    
    tick();
    for ( index_t i=0; i<cnt; ++i )
    {
        randomize(NVAL, x, y, z, 1.0);
        for ( index_t j=0; j<SUB; ++j )
        {
            FUNC(NSEG, dir_, x, lag_, z);
            FUNC(NSEG, dir_, y, lag_, x);
            FUNC(NSEG, dir_, z, lag_, y);
        }
    }
    printf("  %8s %5.0f", str, tock());
    print_floating_point_exceptions("");
    printf("\n");
    free_reals(x,y,z);
}


template < void (*FUNC)(index_t, const real*, real*) >
void timeProject(index_t cnt, char const* str)
{
    real *x = nullptr, *y = nullptr, *z = nullptr;
    new_nans(ALOC, x, y, z);
    randomize(NVAL, x, y, z, 1.0);

    FUNC(NSEG, force_, x);
    VecPrint::edges(NVAL+2, x, 4);
    print_floating_point_exceptions("");

    tick();
    for ( index_t i=0; i<cnt; ++i )
    {
        randomize(NVAL, x, y, z, 1.0);
        for ( index_t j=0; j<SUB; ++j )
        {
            FUNC(NSEG, x, y);
            FUNC(NSEG, y, z);
            FUNC(NSEG, z, x);
        }
    }
    printf("  %8s %5.0f", str, tock());
    print_floating_point_exceptions("");
    printf("\n");
    free_reals(x,y,z);
}

//------------------------------------------------------------------------------
#pragma mark - Speed tests

void testProjectionU(index_t cnt)
{
    std::cout << DIM << "D testProjection UP -- " << NSEG << " segments\n";
    testU<projectForcesU_>(cnt,    "U_   ");
    testU<projectForcesU_PTR>(cnt, "U_PTR");
    testU<projectForcesU_TWO>(cnt, "U_TWO");
#if ( DIM == 2 ) && REAL_IS_DOUBLE && USE_SIMD
    testU<projectForcesU2D_SSE>(cnt, "U_SSE");
#endif
#if ( DIM == 2 ) && REAL_IS_DOUBLE && defined(__AVX__)
    testU<projectForcesU2D_AVX>(cnt, "U_AVX");
    testU<projectForcesU2D_AVY>(cnt, "U_AVY");
#endif
#if ( DIM == 3 ) && REAL_IS_DOUBLE && USE_SIMD
    testU<projectForcesU3D_SSE>(cnt, "U_SSE");
#endif
#if ( DIM == 3 ) && REAL_IS_DOUBLE && defined(__AVX__)
    testU<projectForcesU3D_AVX>(cnt, "U_AVX");
#endif
}

void testProjectionD(index_t cnt)
{
    std::cout << DIM << "D testProjection DOWN -- " << NSEG << " segments\n";
    projectForcesU_(NSEG, dir_, force_, lag_);
    testD<projectForcesD_>(cnt,    "D_   ");
    testD<projectForcesD_ADD>(cnt, "D_ADD");
    testD<projectForcesD_REV>(cnt, "D_REV");
    testD<projectForcesD_RIV>(cnt, "D_RIV");
    testD<projectForcesD_FMA>(cnt, "D_FMA");
    testD<projectForcesD_PTR>(cnt, "D_PTR");
#if REAL_IS_DOUBLE && ( DIM == 3 ) && USE_SIMD
    testD<projectForcesD3D_SSE>(cnt, "D_SSE");
#endif
#if REAL_IS_DOUBLE
#if ( DIM == 3 ) && USE_SIMD
    testD<projectForcesD3D_sse>(cnt, "D_sse");
#endif
#if ( DIM == 3 ) && defined(__AVX__)
    testD<projectForcesD3D_AVX>(cnt, "D_AVX");
#endif
#if ( DIM == 2 ) && USE_SIMD
    testD<projectForcesD2D_SSE>(cnt, "D_SSE");
#endif
#if ( DIM == 2 ) && defined(__AVX__)
    testD<projectForcesD2D_AVX>(cnt, "D_AVX");
#endif
#endif
}

void testProjection(index_t cnt)
{
    std::cout << DIM << "D projectForces -- " << NSEG << " segments\n";
    setProjection(NSEG);
    setAnisotropy(NSEG);
    timeProject<projectForces>(cnt,    "projF");
#if REAL_IS_DOUBLE
#  if USE_SIMD
    timeProject<projectForces_SSE>(cnt, "prSSE");
#  endif
#  if defined(__AVX__)
    timeProject<projectForces_AVX>(cnt, "prAVX");
#  endif
#endif
    std::cout << "misc:\n";
    timeProject<onlyDPTTS>(cnt, "*dptts");
}

void testProjectionDiff(index_t cnt)
{
    setFilament(NSEG, 0.2, 20.0, 1.0);
    setProjection(NSEG);
    std::cout << DIM << "D addProjectionDiff -- " << NSEG << " segments\n";
    timeProject<projectDiff_>(cnt, "pdiff");
    timeProject<projectDiff_R>(cnt, "pdifR");
    timeProject<projectDiff_F>(cnt, "pdifF");
#if REAL_IS_DOUBLE && USE_SIMD
    timeProject<projectDiff_SSE>(cnt, "pdSSE");
#endif
#if ( DIM == 2 ) && REAL_IS_DOUBLE && defined(__AVX__)
    timeProject<projectDiff_AVX>(cnt, "pdAVX");
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Main

int main(int argc, char* argv[])
{
    const index_t REP = 4096;
    RNG.seed();
    std::cout << "DIM=" << DIM << "  real " << sizeof(real) << " bytes   " << __VERSION__ << "\n";
    //enable_floating_point_exceptions();
    if ( 0 )
    {
#if REAL_IS_DOUBLE
#  if defined(__AVX__)
        for ( index_t nbs = std::min(NSEG,(index_t)11); nbs > 0; --nbs )
            checkProject<projectForces_AVX>(nbs, "AVX");
#  endif
#  if USE_SIMD
        for ( index_t nbs = std::min(NSEG,(index_t)11); nbs > 0; --nbs )
            checkProject<projectForces_SSE>(nbs, "SSE");
#  endif
#endif
    }
    if ( 1 )
    {
        setFilament(NSEG, 0.2, 20.0, 1.0);
        testProjectionU(REP);
        testProjectionD(REP);
    }
    if ( 1 )
    {
        setFilament(NSEG, 0.2, 20.0, 1.0);
        testProjection(REP);
    }
    if ( 1 )
    {
        setFilament(NSEG, 0.2, 20.0, 1.0);
        testProjectionDiff(REP);
    }

    free_reals(tmp_, lag_, force_, pos_);
    free_reals(dir_, ani_, diag_, upper_);
}
