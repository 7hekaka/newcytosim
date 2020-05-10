//  Cytosim was created by Francois Nedelec.
//  Copyright FJN 2020 Sainsbury Laboratory, Cambridge University

/**
 This will perform:

     Y = X + (TT') X

 Where T is the local direction of the fiber given by `dir`.
 This is used to multiply the tangential component of force `X` by a factor 2,
 without changing the orthogonal components
 */
void scaleTangentially(size_t nbp, const real* src, const real* dir, real* dst)
{
    for ( size_t p = 0; p < nbp; ++p )
    {
        real const* xxx = src + DIM * p;
        real const* ddd = dir + DIM * p;
        real      * yyy = dst + DIM * p;
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


//------------------------------------------------------------------------------
#pragma mark - Rigidity

/*
 In this version the loop is unrolled
 */
void add_rigidity2(const size_t nbt, const real* X, const real rigid, real* Y)
{
    assert_true( X != Y );
    real dx = X[DIM  ] - X[0];
    real dy = X[DIM+1] - X[1];
#if ( DIM > 2 )
    real dz = X[DIM+2] - X[2];
#endif
    
    real * yv = Y;
    const real*const end = X + nbt + DIM;
    
    const real* xv = X+DIM;
    while ( xv < end )
    {
        real d0 = xv[DIM] - xv[0];
        real f0 = rigid * ( d0 - dx );
        dx = d0;
        yv[0    ] -=    f0;
        yv[DIM  ] += f0+f0;
        yv[DIM*2] -=    f0;
        ++yv;
        ++xv;
        
        real d1 = xv[DIM] - xv[0];
        real f1 = rigid * ( d1 - dy );
        dy = d1;
        yv[0    ] -=    f1;
        yv[DIM  ] += f1+f1;
        yv[DIM*2] -=    f1;
        ++yv;
        ++xv;
        
#if ( DIM > 2 )
        real d2 = xv[DIM] - xv[0];
        real f2 = rigid * ( d2 - dz );
        dz = d2;
        Y[0    ] -=    f2;
        Y[DIM  ] += f2*f2;
        Y[DIM*2] -=    f2;
        ++yv;
        ++xv;
#endif
    }
}


/*
 In this version the loop is unrolled, pointers are used
 and further optimization are made by replacing
 ( a0 -2*a1 + a2 ) by (a2-a1)-(a1-a0).
 */
void add_rigidity3(const size_t nbt, const real* X, const real rigid, real* Y)
{
    assert_true( X != Y );
    const real * xn = X + DIM;
    
    real x0 = xn[0];
    real x1 = xn[1];
#if ( DIM >= 3 )
    real x2 = xn[2];
#endif
    
    real d0 = x0 - X[0];
    real d1 = x1 - X[1];
#if ( DIM >= 3 )
    real d2 = x2 - X[2];
#endif
    
    real df0 = 0, of0 = 0, odf0 = 0;
    real df1 = 0, of1 = 0, odf1 = 0;
#if ( DIM >= 3 )
    real df2 = 0, of2 = 0, odf2 = 0;
#endif
    
    xn += DIM;
    
    real * yp = Y;
    real *const end = Y + nbt;
    while ( yp < end )
    {
        real e0 = *xn - x0;
        x0 = *xn;
        ++xn;
        real f0 = rigid * ( e0 - d0 );
        d0      = e0;
        df0     = f0 - of0;
        of0     = f0;
        *yp    += odf0 - df0;
        odf0    = df0;
        ++yp;
        
        real e1 = *xn - x1;
        x1 = *xn;
        ++xn;
        real f1 = rigid * ( e1 - d1 );
        d1      = e1;
        df1     = f1 - of1;
        of1     = f1;
        *yp    += odf1 - df1;
        odf1    = df1;
        ++yp;
        
#if ( DIM >= 3 )
        real e2 = *xn - x2;
        x2 = *xn;
        ++xn;
        real f2 = rigid * ( e2 - d2 );
        d2      = e2;
        df2     = f2 - of2;
        of2     = f2;
        *yp    += odf2 - df2;
        odf2    = df2;
        ++yp;
#endif
    }
    
    yp[0]   += df0 + of0;
    yp[1]   += df1 + of1;
#if ( DIM >= 3 )
    yp[2]   += df2 + of2;
#endif
    
    yp += DIM;
    
    yp[0] -= of0;
    yp[1] -= of1;
#if ( DIM >= 3 )
    yp[2] -= of2;
#endif
}

#if ( DIM == 2 ) && REAL_IS_DOUBLE

#include "simd.h"

#ifdef __SSE3__
/**
 2D implemention using SSE 128bit vector instructions with double precision
 */
void add_rigiditySSE(const size_t nbt, const real* X, const real rigid, real* Y)
{
    vec2 R = set2(rigid);
    real *const end = Y + nbt;

    vec2 nn = load2(X+2);
    vec2 oo = mul2(R, sub2(nn, load2(X)));
    vec2 yy = load2(Y);
    vec2 zz = load2(Y+2);
    
    while ( Y < end )
    {
        vec2 mm = load2(X+4);
        X += 2;
        vec2 dd = mul2(R, sub2(mm, nn));
        vec2 ff = sub2(dd, oo);
        oo = dd;
        nn = mm;
        store2(Y, sub2(yy, ff));
        yy = add2(zz, add2(ff, ff));
        zz = sub2(load2(Y+4), ff);
        Y += 2;
    }
    store2(Y, yy);
    store2(Y+2, zz);
}
#endif


#ifdef __AVX__
/**
 2D implemention using AVX 256bit vector instructions with double precision
 FJN 15.09.2018 -- 17.09.2018
 
 Note that the vectors X and Y are not aligned to memory!
 */
void add_rigidityAVX(const size_t nbt, const real* X, const real rigid, real* Y)
{
    vec4 R = set4(rigid);
    vec4 two = set4(2.0);
    
    real *const end = Y + nbt - 8;
    
    vec4 xxx = loadu4(X);
    vec4 eee = setzero4();

    // process data 8 by 8:
    while ( Y < end )
    {
        vec4 nnn = loadu4(X+4);
        vec4 iii = permute2f128(xxx, nnn, 0x21);
        vec4 ddd = sub4(sub4(nnn, iii), sub4(iii, xxx));
        xxx = loadu4(X+8);
        X += 8;
        vec4 ppp = permute2f128(eee, ddd, 0x21);
        vec4 jjj = permute2f128(nnn, xxx, 0x21);
#ifdef __FMA__
        storeu4(Y, fmadd4(R, fmsub4(two, ppp, add4(eee, ddd)), loadu4(Y)));
#else
        storeu4(Y, add4(mul4(R, sub4(add4(ppp, ppp), add4(eee, ddd))), loadu4(Y)));
#endif
        eee = sub4(sub4(xxx, jjj), sub4(jjj, nnn));
        ppp = permute2f128(ddd, eee, 0x21);
#ifdef __FMA__
        storeu4(Y+4, fmadd4(R, fmsub4(two, ppp, add4(ddd, eee)), loadu4(Y+4)));
#else
        storeu4(Y+4, add4(mul4(R, sub4(add4(ppp, ppp), add4(eee, ddd))), loadu4(Y+4)));
#endif
        Y += 8;
    }

    // process data 4 by 4:
    if ( Y < end+4 )
    {
        vec4 nnn = loadu4(X+4);
        X += 4;
        vec4 iii = permute2f128(xxx, nnn, 0x21);
        vec4 ddd = sub4(sub4(nnn, iii), sub4(iii, xxx));
        xxx = nnn;
        vec4 ppp = permute2f128(eee, ddd, 0x21);
#ifdef __FMA__
        storeu4(Y, fmadd4(R, fmsub4(two, ppp, add4(eee, ddd)), loadu4(Y)));
#else
        storeu4(Y, add4(mul4(R, sub4(add4(ppp, ppp), add4(eee, ddd))), loadu4(Y)));
#endif
        eee = ddd;
        Y += 4;
    }

    // process data 2 by 2 using SSE instructions:
    vec2 nn = gethi(xxx);
    vec2 oo = sub2(nn, getlo(xxx));
    vec2 ee = gethi(eee);
    vec2 yy = fmsub2(getlo(two), ee, getlo(eee));
    while ( Y < end+8 )
    {
        vec2 mm = loadu2(X+4);
        X += 2;
        vec2 ff = sub2(mm, nn);
        vec2 dd = sub2(ff, oo);
        nn = mm;
        oo = ff;
        storeu2(Y, fmadd2(getlo(R), sub2(yy, dd), loadu2(Y)));
#ifdef __FMA__
        yy = fmsub2(getlo(two), dd, ee);
#else
        yy = sub2(add2(dd, dd), ee);
#endif
        ee = dd;
        Y += 2;
    }
    storeu2(Y  ,  fmadd2(getlo(R), yy, loadu2(Y  )));
    storeu2(Y+2, fnmadd2(getlo(R), ee, loadu2(Y+2)));
}
#endif
#endif

//------------------------------------------------------------------------------
#pragma mark - ProjectForces


/**
 Perform second calculation needed by projectForces:
 */
void projectForcesD__(size_t nbs, const real* dif,
                      const real* src, const real* mul, real* dst)
{
    real a0 = dif[0] * mul[0];
    real a1 = dif[1] * mul[0];
#if ( DIM > 2 )
    real a2 = dif[2] * mul[0];
#endif
    
    dst[0] = src[0] + a0;
    dst[1] = src[1] + a1;
#if ( DIM > 2 )
    dst[2] = src[2] + a2;
#endif
    
    for ( size_t jj = 1; jj < nbs; ++jj )
    {
        const size_t kk = DIM * jj;
        real b0 = dif[kk  ] * mul[jj];
        dst[kk  ] = src[kk  ] + b0 - a0;
        a0 = b0;
        
        real b1 = dif[kk+1] * mul[jj];
        dst[kk+1] = src[kk+1] + b1 - a1;
        a1 = b1;
        
#if ( DIM > 2 )
        real b2 = dif[kk+2] * mul[jj];
        dst[kk+2] = src[kk+2] + b2 - a2;
        a2 = b2;
#endif
    }
    
    const size_t ee = DIM * nbs;
    dst[ee  ] = src[ee  ] - a0;
    dst[ee+1] = src[ee+1] - a1;
#if ( DIM > 2 )
    dst[ee+2] = src[ee+2] - a2;
#endif
}


/**
 Perform first calculation needed by projectForces:
 */
void projectForcesU_PTR(size_t nbs, const real* dif, const real* src, real* mul)
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
void projectForcesU2D_SSE(size_t nbs, const real* dif, const real* src, real* mul)
{
    real const*const end = mul - 1 + nbs;

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
        //store1(mul, hadd2(a, a));
        store1(mul, add2(a, unpackhi2(a, a)));
    }
}

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

#if defined(__AVX__) && REAL_IS_DOUBLE

#include "simd.h"

/**
 Perform first calculation needed by projectForces

 F. Nedelec, 9.12.2016, 6.9.2018
 */
void projectForcesU2D_AVX(size_t nbs, const real* dif, const real* src, real* mul)
{
    real const*const end = mul - 3 + nbs;

    while ( mul < end )
    {
        vec4 a = mul4(sub4(loadu4(src+2), loadu4(src  )), load4(dif  ));
        vec4 b = mul4(sub4(loadu4(src+6), loadu4(src+4)), load4(dif+4));
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
        store1(mul, add2(a, permute2(a,1)));
    }
}

/**
 Perform second calculation needed by projectForces

 ATTENTION: memory X and Y are not necessarily aligned since they are chunck from
 an array containing contiguous coordinates
 F. Nedelec, 9.12.2016, 23.03.2018
 */
void projectForcesD2D_AVX(size_t nbs, const real* dif,
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

//------------------------------------------------------------------------------
#pragma mark - 3D SIMD

#if defined(__AVX__) && REAL_IS_DOUBLE

#include "simd.h"

/// FJN @ Strasbourg, 17 and 18.04.2020
void projectForcesU3D_AVX(size_t nbs, const real* dif, const real* src, real* mul)
{
    const real *const end = mul - 3 + nbs;
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

        vec2 mm = shuffle2(xy, s1, 0b01);
        zx = add2(blend2(zx, s1, 0b10), blend2(zx, xy, 0b01));
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




/*
 Ugly piece of code to harvest AVX power...
 FJN @ Strasbourg, 18 and 19.04.2020
 */
void projectForcesD3D_AVX(size_t nbs, const real* dif, const real* src, const real* mul, real* dst)
{
    const real* const end = mul - 3 + nbs;
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
        a0 = fnmadd4(m0, loadu4(dif-3), a0);
        a1 = fnmadd4(m1, loadu4(dif+1), a1);
        a2 = fnmadd4(m2, loadu4(dif+5), a2);
        storeu4(dst  , a0);
        storeu4(dst+4, a1);
        storeu4(dst+8, a2);
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
        dd = blend4(mm, broadcast1(dif), 0b1000);
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
            
            mul += 3;
            vec4 a0 = fmadd4(p0, load4(dif  ), loadu4(src  ));
            vec4 a1 = fmadd4(p1, load4(dif+4), loadu4(src+4));
            vec4 a2 = fmadd4(p2, broadcast1(dif+8), loadu4(src+8));
            a0 = fnmadd4(m0, dd, a0);
            a1 = fnmadd4(m1, loadu4(dif+1), a1);
            a2 = fnmadd4(m2, loadu4(dif+5), a2);
            storeu4(dst  , a0);
            storeu4(dst+4, a1);
            storeu4(dst+8, a2);
            //dif += 12; dst += 12; src += 12;
        } break;
        case 1: {
            // 3 vectors remaining
            vec4 m1 = broadcast1(mul  );
            vec4 m0 = blend4(mm, m1, 0b1000);
            vec4 m2 = broadcast1(mul+1);
            vec4 p0 = blend4(m1, m2, 0b1000);
            m1 = blend4(m1, m2, 0b1100);
            vec4 p1 = blend4(m2, setzero4(), 0b1100);
            
            mul += 2;
            vec4 a0 = fmadd4(p0, load4(dif  ), loadu4(src  ));
            vec4 a1 = fmadd4(p1, load2Z(dif+4), loadu4(src+4));
            a0 = fnmadd4(m0, dd, a0);
            a1 = fnmadd4(m1, loadu4(dif+1), a1);
            vec2 a2 = fnmadd2(getlo(m2), load1(dif+5), load1(src+8));
            storeu4(dst  , a0);
            storeu4(dst+4, a1);
             store1(dst+8, a2);
            //dif += 9; dst += 9; src += 9;
        } break;
        case 2: {
            // 2 vectors remaining
            vec4 m1 = broadcast1(mul);
            vec4 m0 = blend4(mm, m1, 0b1000);
            vec4 p0 = blend4(m1, setzero4(), 0b1000);
            
            mul += 1;
            vec4 a0 = fmadd4(p0, load3(dif), loadu4(src));
            a0 = fnmadd4(m0, dd, a0);
            vec2 a2 = fnmadd2(getlo(m1), loadu2(dif+1), loadu2(src+4));
            storeu4(dst  , a0);
            storeu2(dst+4, a2);
            //dif += 6; dst += 6; src += 6;
        } break;
        case 3: {
            // 1 vector remaining
            //store3(dst, fnmadd4(mm, blend4(dd, setzero4(), 0b1000), load3(src)));
            storeu2(dst, fnmadd2(getlo(mm), getlo(dd), loadu2(src)));
            store1(dst+2, fnmadd1(getlo(mm), gethi(dd), load1(src+2)));
            //dif += 3; dst += 3; src += 3;
        } break;
        default:
            printf("unexpected case in projectForcesD3D_AVX!");
    }
    assert_true( mul == end+3 );
}

#endif

//------------------------------------------------------------------------------
#pragma mark - ProjectionDiff


///expanded implementation:
void add_projectiondiffR(const size_t nbs, const real* mul, const real* X, real* Y)
{
    // this loop cannot be unrolled as there is an OUTPUT dependency in Y
    for ( size_t jj = 0; jj < nbs; ++jj )
    {
        const real x = mul[jj] * ( X[DIM*jj+DIM  ] - X[DIM*jj  ] );
        Y[DIM*jj      ] += x;
        Y[DIM*jj+DIM  ] -= x;
#if ( DIM > 1 )
        const real y = mul[jj] * ( X[DIM*jj+DIM+1] - X[DIM*jj+1] );
        Y[DIM*jj    +1] += y;
        Y[DIM*jj+DIM+1] -= y;
#endif
#if ( DIM > 2 )
        const real z = mul[jj] * ( X[DIM*jj+DIM+2] - X[DIM*jj+2] );
        Y[DIM*jj    +2] += z;
        Y[DIM*jj+DIM+2] -= z;
#endif
    }
}


///scalar implementation
void add_projectiondiffF(const size_t nbs, const real* mul, const real* X, real* Y)
{
    real px0 = X[0];
    real px1 = X[1];
    real pw0 = 0;
    real pw1 = 0;
#if ( DIM >= 3 )
    real px2 = X[2];
    real pw2 = 0;
#endif
    
    for ( size_t jj = 0; jj < nbs; ++jj )
    {
        const real m = mul[jj];
        real x0 = X[DIM*jj+DIM  ];
        real x1 = X[DIM*jj+DIM+1];
#if ( DIM >= 3 )
        real x2 = X[DIM*jj+DIM+2];
#endif
        real w0 = m * ( x0 - px0 );
        real w1 = m * ( x1 - px1 );
#if ( DIM >= 3 )
        real w2 = m * ( x2 - px2 );
#endif
        px0 = x0;
        px1 = x1;
#if ( DIM >= 3 )
        px2 = x2;
#endif
        Y[DIM*jj  ] += w0 - pw0;
        Y[DIM*jj+1] += w1 - pw1;
#if ( DIM >= 3 )
        Y[DIM*jj+2] += w2 - pw2;
#endif
        pw0 = w0;
        pw1 = w1;
#if ( DIM >= 3 )
        pw2 = w2;
#endif
    }
    Y[DIM*nbs  ] -= pw0;
    Y[DIM*nbs+1] -= pw1;
#if ( DIM >= 3 )
    Y[DIM*nbs+2] -= pw2;
#endif
}

#if ( DIM == 2 ) && defined(__SSE3__) && REAL_IS_DOUBLE

#include "simd.h"

void add_projectiondiffSSE(const size_t nbs, const real* mul, const real* X, real* Y)
{
    vec2 px = load2(X);
    vec2 pw = setzero2();
    
    for ( size_t jj = 0; jj < nbs; ++jj )
    {
        vec2 m = loaddup2(mul+jj);
        vec2 x = load2(X+DIM*jj+DIM);
        vec2 y = load2(Y+DIM*jj);
        vec2 w = mul2(m, sub2(x, px));
        px = x;
        store2(Y+DIM*jj, add2(y, sub2(w, pw)));
        pw = w;
    }
    store2(Y+DIM*nbs, sub2(load2(Y+DIM*nbs), pw));
}

#endif

#if ( DIM == 2 ) && defined(__AVX__) && REAL_IS_DOUBLE

#include "simd.h"

void add_projectiondiffAVX(const size_t nbs, const real* mul, const real* X, real* Y)
{
    real * pY = Y;
    real const* pX = X;
    real const* pM = mul;
    
    if ( nbs & 1 )
    {
        vec2 m = loaddup2(pM);
        ++pM;
        vec2 s = mul2(sub2(load2(pX+DIM), load2(pX)), m);
        pX += DIM;
        storeu2(pY    , add2(load2(pY    ), s));
        storeu2(pY+DIM, sub2(load2(pY+DIM), s));
        pY += DIM;
    }
    
    real const*const end = mul + nbs;
    while ( pM < end )
    {
        vec4 a = broadcast2(pM);
        vec4 m = permute4(a, 0b1100);

        pM += DIM;
        vec4 s = mul4(m, sub4(loadu4(pX+2), loadu4(pX)));
        pX += 2*DIM;
        
        // this will not be fast, since the two vector are not independent:
        storeu4(pY  , add4(loadu4(pY  ), s));
        storeu4(pY+2, sub4(loadu4(pY+2), s));
        pY += 2*DIM;
    }
    assert_true(pM==end);
}

#endif
