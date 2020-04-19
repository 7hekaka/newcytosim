//  Cytosim was created by Francois Nedelec.
//  Copyright FJN 2020 Sainsbury Laboratory, Cambridge University


/**
 This will perform:

     Y = X + (TT') X

 Where T is the local direction of the fiber given by `dir`.
 This is used to multiply the tangential component of force `X` by a factor 2,
 without changing the orthogonal components
 */
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


/**
 Perform second calculation needed by projectForces:
 */
void projectForcesD__(size_t nbs, const real* dif,
                      const real* X, const real* mul, real* Y)
{
    real a0 = dif[0] * mul[0];
    real a1 = dif[1] * mul[0];
#if ( DIM > 2 )
    real a2 = dif[2] * mul[0];
#endif
    
    Y[0] = X[0] + a0;
    Y[1] = X[1] + a1;
#if ( DIM > 2 )
    Y[2] = X[2] + a2;
#endif
    
    for ( size_t jj = 1; jj < nbs; ++jj )
    {
        const size_t kk = DIM * jj;
        real b0 = dif[kk  ] * mul[jj];
        Y[kk  ] = X[kk  ] + b0 - a0;
        a0 = b0;
        
        real b1 = dif[kk+1] * mul[jj];
        Y[kk+1] = X[kk+1] + b1 - a1;
        a1 = b1;
        
#if ( DIM > 2 )
        real b2 = dif[kk+2] * mul[jj];
        Y[kk+2] = X[kk+2] + b2 - a2;
        a2 = b2;
#endif
    }
    
    const size_t ee = DIM * nbs;
    Y[ee  ] = X[ee  ] - a0;
    Y[ee+1] = X[ee+1] - a1;
#if ( DIM > 2 )
    Y[ee+2] = X[ee+2] - a2;
#endif
}


/**
 Perform first calculation needed by projectForces:
 */
inline void projectForcesU_PTR(size_t nbs, const real* dif, const real* src, real* mul)
{
    real x3, x0 = src[0];
    real x4, x1 = src[1];
#if ( DIM >= 3 )
    real x5, x2 = src[2];
#endif
    src += DIM;
    const real *const end = mul + nbs;

    //normally optimized version
    while ( mul < end )
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


//------------------------------------------------------------------------------
#pragma mark - 2D SIMD


#if defined(__SSE3__) && REAL_IS_DOUBLE

#include "simd.h"

/**
 Perform first calculation needed by projectForces:
 */
inline void projectForcesU2D_SSE(size_t nbs, const real* dif, const real* X, real* mul)
{
    const real* pD = dif;
    const real* pX = X;
    real const*const end = mul + nbs - 2;
    real* pT = mul;
    
    vec2 y, x = load2(pX);
    while ( pT <= end )
    {
        y = load2(pX+2);
        pX += 4;
        vec2 a = mul2(sub2(y, x), load2(pD));
        x = load2(pX);
        vec2 b = mul2(sub2(x, y), load2(pD+2));
        pD += 4;
        //storeu2(pT, hadd2(a, b));
        storeu2(pT, add2(unpacklo2(a, b), unpackhi2(a, b)));
        pT += 2;
    }
    
    if ( pT < end+2 )
    {
        y = load2(pX+2);
        vec2 a = mul2(sub2(y, x), load2(pD));
        //storelo(pT, hadd2(a, a));
        storelo(pT, add2(a, unpackhi2(a, a)));
    }
}

/**
 Perform second calculation needed by projectForces:
 */
inline void projectForcesD2D_SSE(size_t nbs, const real* dif,
                               const real* X, const real* mul, real* Y)
{
    real *pY = Y;
    const real* pX = X;
    const real* pD = dif;
    
    vec2 cc = load2(X);
    
    real const* pM = mul;
    real const*const end = mul + nbs;
    while ( pM < end )
    {
        pX += DIM;
        vec2 d = mul2(load2(pD), loaddup2(pM));
        ++pM;
        pD += DIM;
        store2(pY, add2(cc, d));
        pY += DIM;
        cc = sub2(load2(pX), d);
    }
    store2(pY, cc);
}

#endif

#if defined(__AVX__) && REAL_IS_DOUBLE

#include "simd.h"

/**
 Perform first calculation needed by projectForces

 F. Nedelec, 9.12.2016, 6.9.2018
 */
inline void projectForcesU2D_AVX(size_t nbs, const real* dif, const real* X, real* mul)
{
    const real* pD = dif;
    const real* pX = X;
    real const*const end = mul + nbs - 4;
    real* pT = mul;

    while ( pT <= end )
    {
        vec4 a = mul4(sub4(loadu4(pX+2), loadu4(pX  )), load4(pD  ));
        vec4 b = mul4(sub4(loadu4(pX+6), loadu4(pX+4)), load4(pD+4));
        pD += 8;
        pX += 8;
        //store4(pT, hadd4(permute2f128(a,b,0x20), permute2f128(a,b,0x31)));
        vec4 p = permute2f128(a,b,0x20);
        vec4 q = permute2f128(a,b,0x31);
        store4(pT, add4(unpacklo4(p, q), unpackhi4(p, q)));
        pT += 4;
    }
    
    while ( pT <= end+2 )
    {
        vec4 d = mul4(sub4(loadu4(pX+2), loadu4(pX)), load4(pD));
        pX += 4;
        pD += 4;
        vec2 h = gethi(d);
        storeu2(pT, add2(unpacklo2(getlo(d),h), unpackhi2(getlo(d),h)));
        pT += 2;
    }
    
    if ( pT < end+4 )
    {
        vec2 a = mul2(sub2(load2(pX+2), load2(pX)), load2(pD));
        //storelo(pT, hadd2(a, a));
        storelo(pT, add2(a, unpackhi2(a, a)));
    }
}


/**
 Perform second calculation needed by projectForces

 ATTENTION: memory X and Y are not necessarily aligned since they are chunck from
 an array containing contiguous coordinates
 F. Nedelec, 9.12.2016, 23.03.2018
 */
inline void projectForcesD2D_AVX(size_t nbs, const real* dif,
                                 const real* X, const real* mul, real* Y)
{
    real *pY = Y;
    const real* pX = X;
    const real* pD = dif;
    
    vec4 cc = setzero4();
    
    const bool odd = nbs & 1;
    real const* pM = mul;
    real const*const end = mul + nbs - odd;
    
    while ( pM < end )
    {
        vec4 t = broadcast2(pM);
        vec4 x = loadu4(pX);
        pM += 2;
        vec4 m = permute4(t, 0b1100);
        vec4 d = mul4(m, load4(pD));
        pD += 4;
        vec4 n = permute2f128(cc,d,0x21);
        cc = d;
        vec4 z = add4(x, sub4(d, n));
        pX += 4;
        storeu4(pY, z);
        pY += 4;
    }
    
    vec2 c = gethi(cc);
    
    if ( odd )
    {
        assert( pM + 1 == mul + nbs );
        vec2 m = loaddup2(pM);
        vec2 x = mul2(m, load2(pD));
        vec2 z = add2(load2(pX), sub2(x, c));
        storeu2(pY, z);
        c = x;
        pY += 2;
        pX += 2;
    }
    
    vec2 z = sub2(load2(pX), c);
    storeu2(pY, z);
    assert( pY == Y + DIM * nbs );
    assert( pX == X + DIM * nbs );
}

#endif

//------------------------------------------------------------------------------
#pragma mark - 3D SIMD

#if defined(__AVX__) && REAL_IS_DOUBLE

#include "simd.h"

/// FJN @ Strasbourg, 17 and 18.04.2020
void projectForcesU3D_AVX(size_t nbs, const real* dif, const real* src, real* mul)
{
    const real *const end = mul + nbs - 4;
    while ( mul <= end )
    {
        /*
         *mul = dif[0] * ( src[DIM  ] - src[0] )
              + dif[1] * ( src[DIM+1] - src[1] )
              + dif[2] * ( src[DIM+2] - src[2] );
         */
        vec4 s0 = mul4(load4(dif  ), sub4(loadu4(src+ 3), load4(src  )));
        vec4 s1 = mul4(load4(dif+4), sub4(loadu4(src+ 7), load4(src+4)));
        vec4 s2 = mul4(load4(dif+8), sub4(loadu4(src+11), load4(src+8)));

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
    while ( mul <= end+2 )
    {
        /*
         *mul = dif[0] * ( src[DIM  ] - src[0] )
              + dif[1] * ( src[DIM+1] - src[1] )
              + dif[2] * ( src[DIM+2] - src[2] );
         */
        vec4 s0 = mul4(load4(dif  ), sub4(loadu4(src+ 3), load4(src  )));
        vec4 s1 = mul4(load4(dif+4), sub4(loadu4(src+ 7), load4(src+4)));

        vec4 xy = blend4(s0, s1, 0b1100);
        vec4 zx = permute2f128(s0, s0, 0x21);
        
        vec4 mm = shuffle4(xy, s1, 0b0101);
        mm = add4(mm, blend4(zx, xy, 0b0101));
        store2(mul, add4(mm, blend4(zx, s1, 0b1010)));
        
        src += 6;
        dif += 6;
        mul += 2;
    }
    while ( mul < end+4 )
    {
        /*
         *mul = dif[0] * ( src[DIM  ] - src[0] )
              + dif[1] * ( src[DIM+1] - src[1] )
              + dif[2] * ( src[DIM+2] - src[2] );
         */
        vec4 x = mul4(load4(dif), sub4(loadu4(src+3), load4(src)));
        vec4 z = permute2f128(x, x, 0x21);
        vec4 y = permute4(x, 0b0011);
        
        store1(mul, add4(add4(x, y), z));
        
        src += 3;
        dif += 3;
        ++mul;
    }
}


/*
 Ugly piece of code to harvest AVX power...
 FJN @ Strasbourg, 18 and 19.04.2020
 */
void projectForcesD3D_AVX(size_t nbs, const real* dif, const real* src, const real* mul, real* dst)
{
    const real* const end = mul + nbs - 4;
    /*
     This follows the standard pattern defined below, except
     that the negative terms are not present on the first vector.
     This handles 12 scalars (4 vectors) in one round.
     */
    if ( mul <= end )
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
        vec4 a0 = fmadd4(p0, load4(dif  ), load4(src  ));
        vec4 a1 = fmadd4(p1, load4(dif+4), load4(src+4));
        vec4 a2 = fmadd4(p2, load4(dif+8), load4(src+8));

        store4(dst  , fnmadd4(m0, broadcast1(dif), a0));
        store4(dst+4, fnmadd4(m1, loadu4(dif+1), a1));
        store4(dst+8, fnmadd4(m2, loadu4(dif+5), a2));
        dif += 12;
        dst += 12;
        src += 12;
    }
    /*
     This is where the bulk of the work is done, handing 12 scalars per pass.
     The loop can be unrolled.
     */
    while ( mul <= end )
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
        vec4 a0 = fmadd4(p0, load4(dif  ), load4(src  ));
        vec4 a1 = fmadd4(p1, load4(dif+4), load4(src+4));
        vec4 a2 = fmadd4(p2, load4(dif+8), load4(src+8));

        store4(dst  , fnmadd4(m0, loadu4(dif-3), a0));
        store4(dst+4, fnmadd4(m1, loadu4(dif+1), a1));
        store4(dst+8, fnmadd4(m2, loadu4(dif+5), a2));
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
        case 1: {
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
            vec4 a0 = fmadd4(p0, load4(dif  ), load4(src  ));
            vec4 a1 = fmadd4(p1, load4(dif+4), load4(src+4));
            vec4 a2 = fmadd4(p2, broadcast1(dif+8), load4(src+8));
            
            store4(dst  , fnmadd4(m0, dd, a0));
            store4(dst+4, fnmadd4(m1, loadu4(dif+1), a1));
            store4(dst+8, fnmadd4(m2, loadu4(dif+5), a2));
            //dif += 12; dst += 12; src += 12;
        } break;
        case 2: {
            // 3 vectors remaining
            vec4 m1 = broadcast1(mul  );
            vec4 m0 = blend4(mm, m1, 0b1000);
            vec4 m2 = broadcast1(mul+1);
            vec4 p0 = blend4(m1, m2, 0b1000);
            m1 = blend4(m1, m2, 0b1100);
            vec4 p1 = blend4(m2, setzero4(), 0b1100);
            
            mul += 3;
            vec4 a0 = fmadd4(p0, load4(dif  ), load4(src  ));
            vec4 a1 = fmadd4(p1, load4(dif+4), load4(src+4));
            vec4 a2 = broadcast1(src+8);
            
            store4(dst  , fnmadd4(m0, dd, a0));
            store4(dst+4, fnmadd4(m1, loadu4(dif+1), a1));
            store1(dst+8, fnmadd4(m2, broadcast1(dif+5), a2));
            //storelo(dst+8, fnmadd2(getlo(m2), loaddup2(dif+5), getlo(a2)));
            //dif += 9; dst += 9; src += 9;
        } break;
        case 3: {
            // 2 vectors remaining
            vec4 m1 = broadcast1(mul);
            vec4 m0 = blend4(mm, m1, 0b1000);
            vec4 p0 = blend4(m1, setzero4(), 0b1000);
            
            mul += 2;
            vec4 a0 = fmadd4(p0, load4(dif), load4(src));
            vec4 a1 = load4(src+4);
            
            store4(dst  , fnmadd4(m0, dd, a0));
            storeu2(dst+4, fnmadd4(m1, broadcast2(dif+1), a1));
            //store2(dst+4, fnmadd2(getlo(m1), loadu2(dif+1), getlo(a1)));
            //dif += 6; dst += 6; src += 6;
        } break;
        case 4: {
            // 1 vector remaining
            ++mul;
            store3(dst, fnmadd4(mm, dd, load3(src)));
            //dif += 3; dst += 3; src += 3;
        } break;
        default:
            printf("unexpected case in projectForcesD3D_AVX!");
    }
    assert_true( mul == end+5 );
}

#endif
