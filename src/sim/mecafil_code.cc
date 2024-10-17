//  Cytosim was created by Francois Nedelec.
//  Copyright FJN 2020 Sainsbury Laboratory, Cambridge University

#include "simd.h"
#include "simd_float.h"


//------------------------------------------------------------------------------
#pragma mark - Rigidity

/*
 In this version the loop is unrolled
 */
void add_rigidity2(const int nbp, const real* X, const real rigid, real* Y)
{
    assert_true( X != Y );
    real dx = X[DIM  ] - X[0];
    real dy = X[DIM+1] - X[1];
#if ( DIM > 2 )
    real dz = X[DIM+2] - X[2];
#endif
    
    real * yv = Y;
    const real two = 2.0;
    const real*const end = X + DIM * ( nbp - 1 );
    
    const real* xv = X+DIM;
    while ( xv < end )
    {
        real d0 = xv[DIM] - xv[0];
        real f0 = rigid * ( d0 - dx );
        dx = d0;
        yv[0    ] -= f0;
        yv[DIM  ] += f0 * two;
        yv[DIM*2] -= f0;
        ++yv;
        ++xv;
        
        real d1 = xv[DIM] - xv[0];
        real f1 = rigid * ( d1 - dy );
        dy = d1;
        yv[0    ] -= f1;
        yv[DIM  ] += f1 * two;
        yv[DIM*2] -= f1;
        ++yv;
        ++xv;
        
#if ( DIM > 2 )
        real d2 = xv[DIM] - xv[0];
        real f2 = rigid * ( d2 - dz );
        dz = d2;
        Y[0    ] -= f2;
        Y[DIM  ] += f2 * two;
        Y[DIM*2] -= f2;
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
void add_rigidity3(const int nbp, const real* X, const real rigid, real* Y)
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
    real *const end = Y + DIM * ( nbp - 2 );
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

#if USE_SIMD
/**
 2D implemention using SSE 128bit vector instructions with double precision
 */
void add_rigidity2D_SSE(const int nbp, const double* X, const double rigid, double* Y)
{
    vec2 R = set2(rigid);
    double *const end = Y + DIM * ( nbp - 2 );

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
void add_rigidity2D_AVX(const int nbp, const double* X, const double rigid, double* Y)
{
    vec4 R = set4(rigid);
    vec4 two = set4(2.0);
    
    double *const end = Y + DIM * ( nbp - 2 ) - 8;
    
    vec4 xxx = loadu4(X);
    vec4 eee = setzero4();

    // process data 8 by 8:
    while ( Y < end )
    {
        vec4 nnn = loadu4(X+4);
        vec4 iii = catshift2(xxx, nnn);
        vec4 ddd = sub4(sub4(nnn, iii), sub4(iii, xxx));
        xxx = loadu4(X+8);
        X += 8;
        vec4 ppp = catshift2(eee, ddd);
        vec4 jjj = catshift2(nnn, xxx);
#if defined(__FMA__)
        storeu4(Y, fnmadd4(R, fnmadd4(two, ppp, add4(eee, ddd)), loadu4(Y)));
#else
        storeu4(Y, add4(mul4(R, sub4(add4(ppp, ppp), add4(eee, ddd))), loadu4(Y)));
#endif
        eee = sub4(sub4(xxx, jjj), sub4(jjj, nnn));
        ppp = catshift2(ddd, eee);
#if defined(__FMA__)
        storeu4(Y+4, fnmadd4(R, fnmadd4(two, ppp, add4(ddd, eee)), loadu4(Y+4)));
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
        vec4 iii = catshift2(xxx, nnn);
        vec4 ddd = sub4(sub4(nnn, iii), sub4(iii, xxx));
        xxx = nnn;
        vec4 ppp = catshift2(eee, ddd);
#if defined(__FMA__)
        storeu4(Y, fnmadd4(R, fnmadd4(two, ppp, add4(eee, ddd)), loadu4(Y)));
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
    vec2 yy = fnmadd2(getlo(two), ee, getlo(eee));
    while ( Y < end+8 )
    {
        vec2 mm = loadu2(X+4);
        X += 2;
        vec2 ff = sub2(mm, nn);
        vec2 dd = sub2(ff, oo);
        nn = mm;
        oo = ff;
        storeu2(Y, fmadd2(getlo(R), add2(ee, yy), loadu2(Y)));
#if defined(__FMA__)
        yy = fnmadd2(getlo(two), dd, ee);
#else
        yy = sub2(ee, add2(dd, dd));
#endif
        ee = dd;
        Y += 2;
    }
    storeu2(Y  , fnmadd2(getlo(R), yy, loadu2(Y  )));
    storeu2(Y+2, fnmadd2(getlo(R), ee, loadu2(Y+2)));
}
#endif
#endif

//------------------------------------------------------------------------------
#pragma mark - ProjectForces


/**
 Performs second calculation needed by projectForces.
 This version only calculates the difference of `src` once
 */
void projectForcesD__(size_t nbs, const real* dir,
                      const real* src, const real* mul, real* dst)
{
    real a0 = dir[0] * mul[0];
    real a1 = dir[1] * mul[0];
#if ( DIM > 2 )
    real a2 = dir[2] * mul[0];
#endif
    
    dst[0] = src[0] + a0;
    dst[1] = src[1] + a1;
#if ( DIM > 2 )
    dst[2] = src[2] + a2;
#endif
    
    for ( size_t jj = 1; jj < nbs; ++jj )
    {
        const size_t kk = DIM * jj;
        real b0 = dir[kk  ] * mul[jj];
        dst[kk  ] = src[kk  ] + b0 - a0;
        a0 = b0;
        
        real b1 = dir[kk+1] * mul[jj];
        dst[kk+1] = src[kk+1] + b1 - a1;
        a1 = b1;
        
#if ( DIM > 2 )
        real b2 = dir[kk+2] * mul[jj];
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
void projectForcesU_PTR(index_t nbs, const real* dir, const real* src, real* mul)
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
}


/**
 Perform second calculation needed by projectForces:
 */
void projectForcesD_PTR(index_t nbs, const real* dir,
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
        real b0 = dir[0] * mul[0];
        real b1 = dir[1] * mul[0];
#if ( DIM > 2 )
        real b2 = dir[2] * mul[0];
#endif
        dir += DIM;
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


#if REAL_IS_DOUBLE && USE_SIMD

/**
 Perform first calculation needed by projectForces:
 */
void projectForcesU2D_SSE(index_t nbs, const double* dir, const double* src, double* mul)
{
    real const*const end = mul - 1 + nbs;

    vec2 y, x = load2(src);
    while ( mul < end )
    {
        y = load2(src+2);
        src += 4;
        vec2 a = mul2(sub2(y, x), load2(dir));
        x = load2(src);
        vec2 b = mul2(sub2(x, y), load2(dir+2));
        dir += 4;
        storeu2(mul, add2(unpacklo2(a, b), unpackhi2(a, b)));
        mul += 2;
    }
    
    if ( mul < end+1 )
    {
        y = load2(src+2);
        vec2 a = mul2(sub2(y, x), load2(dir));
        store1(mul, add1(a, unpackhi2(a, a)));
    }
}

/**
 Perform second calculation needed by projectForces:
 */
void projectForcesD2D_SSE(index_t nbs, const double* dir,
                          const double* src, const double* mul, double* dst)
{
    vec2 cc = load2(src);
    
    real const*const end = mul + nbs;
    while ( mul < end )
    {
        src += DIM;
        vec2 d = mul2(load2(dir), loaddup2(mul));
        ++mul;
        dir += DIM;
        store2(dst, add2(cc, d));
        dst += DIM;
        cc = sub2(load2(src), d);
    }
    store2(dst, cc);
}

#endif

#if defined(__AVX__) && REAL_IS_DOUBLE

/**
 Perform first calculation needed by projectForces

 F. Nedelec, 9.12.2016, 6.9.2018
 */
void projectForcesU2D_AVX(index_t nbs, const double* dir, const double* src, double* mul)
{
    real const*const end = mul - 3 + nbs;

    while ( mul < end )
    {
        vec4 a = mul4(sub4(loadu4(src+2), loadu4(src  )), load4(dir  ));
        vec4 b = mul4(sub4(loadu4(src+6), loadu4(src+4)), load4(dir+4));
        dir += 8;
        src += 8;
        vec4 p = unpacklo2f128(a,b);
        vec4 q = unpackhi2f128(a,b);
        storeu4(mul, add4(unpacklo4(p, q), unpackhi4(p, q)));
        mul += 4;
    }

    while ( mul < end+2 )
    {
        //mul[jj] = dir[0] * ( X[DIM] - X[0] ) + dir[1] * ( X[DIM+1] - X[1] )
        vec4 d = mul4(sub4(loadu4(src+2), loadu4(src)), load4(dir));
        src += 4;
        dir += 4;
        vec2 h = gethi(d);
        storeu2(mul, add2(unpacklo2(getlo(d),h), unpackhi2(getlo(d),h)));
        mul += 2;
    }
    
    if ( mul < end+3 )
    {
        vec2 a = mul2(sub2(load2(src+2), load2(src)), load2(dir));
        store1(mul, add1(a, permute2(a)));
    }
}

/**
 Perform second calculation needed by projectForces

 ATTENTION: memory X and Y are not necessarily aligned since they are chunck from
 an array containing contiguous coordinates
 F. Nedelec, 9.12.2016, 23.03.2018
 */
void projectForcesD2D_AVX(index_t nbs, const double* dir,
                          const double* src, const double* mul, double* dst)
{
    vec4 cc = setzero4();
    
    const bool odd = ( nbs & 1 );
    real const*const end = mul - odd + nbs;
    
    while ( mul < end )
    {
        vec4 t = broadcast2(mul);
        vec4 x = loadu4(src);
        mul += 2;
        vec4 m = duplohi4(t);
        vec4 d = mul4(m, load4(dir));
        dir += 4;
        vec4 n = catshift2(cc,d);
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
        vec2 x = mul2(m, load2(dir));
        dir += 2;
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

#if REAL_IS_DOUBLE && USE_SIMD

void projectForcesU3D_SSE(index_t nbs, const double* dir, const double* src, double* mul)
{
    const double *const end = mul + nbs - 1;
#if 1
    // unrolled bulk of the calculation, processing 4x3 scalars
    while ( mul < end-2 )
    {
        /*
         *mul = dir[0] * ( src[DIM  ] - src[0] )
              + dir[1] * ( src[DIM+1] - src[1] )
              + dir[2] * ( src[DIM+2] - src[2] );
         */
        vec2 s0 = mul2(load2(dir  ), sub2(loadu2(src+3), loadu2(src  )));
        vec2 s1 = mul2(load2(dir+2), sub2(loadu2(src+5), loadu2(src+2)));
        vec2 s2 = mul2(load2(dir+4), sub2(loadu2(src+7), loadu2(src+4)));
        vec2 ss = unpacklo2(s0, s2);
        vec2 tt = unpackhi2(s0, s2);

        vec2 t0 = mul2(load2(dir+ 6), sub2(loadu2(src+ 9), loadu2(src+ 6)));
        vec2 t1 = mul2(load2(dir+ 8), sub2(loadu2(src+11), loadu2(src+ 8)));
        vec2 t2 = mul2(load2(dir+10), sub2(loadu2(src+13), loadu2(src+10)));
        vec2 xx = unpacklo2(t0, t2);
        vec2 yy = unpackhi2(t0, t2);

        store2(mul, add2(tt, add2(s1, ss)));
        store2(mul+2, add2(yy, add2(t1, xx)));

        src += 12;
        dir += 12;
        mul += 4;
    }
#endif
    // bulk of the calculation, processing 2x3 scalars
    while ( mul < end )
    {
        /*
         *mul = dir[0] * ( src[DIM  ] - src[0] )
              + dir[1] * ( src[DIM+1] - src[1] )
              + dir[2] * ( src[DIM+2] - src[2] );
         */
        vec2 s0 = mul2(load2(dir  ), sub2(loadu2(src+3), loadu2(src  )));
        vec2 s1 = mul2(load2(dir+2), sub2(loadu2(src+5), loadu2(src+2)));
        vec2 s2 = mul2(load2(dir+4), sub2(loadu2(src+7), loadu2(src+4)));

        vec2 ss = unpacklo2(s0, s2);
        vec2 tt = unpackhi2(s0, s2);

        store2(mul, add2(tt, add2(s1, ss)));
        
        src += 6;
        dir += 6;
        mul += 2;
    }
    // process remaining vector, 3 scalars
    if ( mul <= end )
    {
        /*
         *mul = dir[0] * ( src[DIM  ] - src[0] )
         + dir[1] * ( src[DIM+1] - src[1] )
         + dir[2] * ( src[DIM+2] - src[2] );
         */
        vec2 s0 = mul2(loadu2(dir ), sub2(loadu2(src+3), loadu2(src)));
        vec2 s1 = mul1(load1(dir+2), sub1(load1(src+5), load1(src+2)));
        vec2 ss = unpackhi2(s0, setzero2());
        
        store1(mul, add1(ss, add1(s0, s1)));
        //src += 3;
        //dir += 3;
        //mul += 1;
   }
}


/*
 3D workflow:
     Y[0] = dir[0] * mul[0] + "X[0] - dir[-3] * mul[-1]";
     Y[1] = dir[1] * mul[0] + "X[1] - dir[-2] * mul[-1]";

     Y[2] = dir[2] * mul[0] + "X[2] - dir[-1] * mul[-1]";
     Y[3] = dir[3] * mul[1] + X[3] - dir[0] * mul[0];

     Y[4] = dir[4] * mul[1] + X[4] - dir[1] * mul[0];
     Y[5] = dir[5] * mul[1] + X[5] - dir[2] * mul[0];
 */
void projectForcesD3D_SSE(index_t nbs, const double* dir,
                          const double* src, const double* mul, double* dst)
{
    double const*const end = mul + nbs - 1;
    vec2 p0 = loadu2(src);
    vec2 p1 = loadu2(src+2);
    
    if ( mul >= end )
    {
        // special case with only 2 points overall and one multiplicator
        vec2 p2 = loadu2(src+4);
        vec2 m0 = loaddup2(mul);
        vec2 x0 = mul2(load2(dir), m0);
        vec2 x1 = mul1(load1(dir+2), m0); // only lower value used
        vec2 a0 = add2(p0, x0);
        vec2 a1 = add1(p1, x1); // only lower value used
        vec2 s1 = sub2(p1, catshift(x0, x0)); // only upper value used
        vec2 s2 = sub2(p2, catshift(x0, x1));
        storeu2(dst  , a0);
        storeu2(dst+2, blend11(a1, s1));
        storeu2(dst+4, s2);
        return;
    }
#if defined(__AVX__)
    vec4 p01 = concatenate22(p0, p1);
    while ( mul < end-3 )
    {
        // processing 12 scalars = 6 SSE vectors
        vec4 p23 = loadu4(src+4);
        vec4 p45 = loadu4(src+8);
        vec4 m0 = loadu4(mul);
        vec4 m1 = duplo2f128(m0);
        vec4 m3 = duphi2f128(m0);
        m0 = duplo4(m1);
        m1 = duphi4(m1);
        vec4 m2 = duplo4(m3);
        m3 = duphi4(m3);
        vec4 x01 = mul4(load4(dir  ), blend31(m0, m1));
        vec4 x23 = mul4(load4(dir+4), blend22(m1, m2));
        vec4 x45 = mul4(load4(dir+8), blend13(m2, m3));
        storeu4(dst  , add4(x01, sub4(p01, catshift1(setzero4(), x01))));
        storeu4(dst+4, add4(x23, sub4(p23, catshift1(x01, x23))));
        storeu4(dst+8, add4(x45, sub4(p45, catshift1(x23, x45))));
        dir += 12;
        dst += 12;
        src += 12;
        mul += 4;
        p01 = sub4(loadu4(src), catshift1(x45, setzero4())); // may load crap at src+2
    }
    p0 = getlo(p01);
    p1 = gethi(p01);
#elif 1
    while ( mul < end-3 )
    {
        // processing 12 scalars = 6 SSE vectors
        vec2 p2 = loadu2(src+4);
        vec2 p3 = loadu2(src+6);
        vec2 p4 = loadu2(src+8);
        vec2 p5 = loadu2(src+10);
        vec2 m0 = loaddup2(mul);
        vec2 m1 = loaddup2(mul+1);
        vec2 m2 = loaddup2(mul+2);
        vec2 m3 = loaddup2(mul+3);
        vec2 x0 = mul2(load2(dir  ), m0);
        vec2 x1 = mul2(load2(dir+2), catshift(m0, m1));
        vec2 x2 = mul2(load2(dir+4), m1);
        vec2 x3 = mul2(load2(dir+6), m2);
        vec2 x4 = mul2(load2(dir+8), catshift(m2, m3));
        vec2 x5 = mul2(load2(dir+10), m3);
        storeu2(dst  , add2(x0, p0));
        storeu2(dst+2, add2(x1, sub2(p1, catshift(setzero2(), x0))));
        storeu2(dst+4, add2(x2, sub2(p2, catshift(x0, x1))));
        storeu2(dst+6, add2(x3, sub2(p3, catshift(x1, x2))));
        storeu2(dst+8, add2(x4, sub2(p4, catshift(x2, x3))));
        storeu2(dst+10, add2(x5, sub2(p5, catshift(x3, x4))));
        dir += 12;
        dst += 12;
        src += 12;
        mul += 4;
        p0 = sub2(loadu2(src), catshift(x4, x5));
        p1 = sub2(loadu2(src+2), catshift(x5, setzero2()));  // may load crap at src+2
    }
#endif
    while ( mul < end )
    {
        // processing 6 scalars = 3 SSE vectors
        vec2 p2 = loadu2(src+4);
        vec2 m0 = loaddup2(mul);
        vec2 m1 = loaddup2(mul+1);
        vec2 x0 = mul2(load2(dir  ), m0);
        vec2 x1 = mul2(load2(dir+2), catshift(m0, m1));
        vec2 x2 = mul2(load2(dir+4), m1);
        storeu2(dst  , add2(x0, p0));
        storeu2(dst+2, add2(x1, sub2(p1, catshift(setzero2(), x0))));
        storeu2(dst+4, add2(x2, sub2(p2, catshift(x0, x1))));
        dir += 6;
        dst += 6;
        src += 6;
        mul += 2;
        p0 = sub2(loadu2(src), catshift(x1, x2));
        p1 = sub2(loadu2(src+2), catshift(x2, setzero2()));  // may load crap at src+2
    }
    if ( nbs & 1 )
    {
        // two points remaining, 6 scalars
        vec2 p2 = loadu2(src+4);
        vec2 m0 = loaddup2(mul);
        vec2 x0 = mul2(load2(dir), m0);
        vec2 x1 = mul1(load1(dir+2), m0); //only lower value used
        storeu2(dst  , add2(x0, p0));
        storeu2(dst+2, add2(x1, sub2(p1, catshift(setzero2(), x0))));
        storeu2(dst+4, sub2(p2, catshift(x0, x1)));
    }
    else
    {
        // last point, 3 scalars
        storeu2(dst, p0);
        store1(dst+2, p1);
    }
}
#endif


#if defined(__AVX__) && REAL_IS_DOUBLE

/// FJN @ Strasbourg, 17 and 18.04.2020
void projectForcesU3D_AVX(index_t nbs, const double* dir, const double* src, double* mul)
{
    const double *const end = mul - 3 + nbs;
    while ( mul < end )
    {
        /*
         *mul = dir[0] * ( src[DIM  ] - src[0] )
              + dir[1] * ( src[DIM+1] - src[1] )
              + dir[2] * ( src[DIM+2] - src[2] );
         */
        vec4 s0 = mul4(load4(dir  ), sub4(loadu4(src+ 3), loadu4(src  )));
        vec4 s1 = mul4(load4(dir+4), sub4(loadu4(src+ 7), loadu4(src+4)));
        vec4 s2 = mul4(load4(dir+8), sub4(loadu4(src+11), loadu4(src+8)));

        vec4 zx = blend22(s2, s0);
        vec4 xy = blend22(s0, s1);
        zx = swap2f128(zx);
        vec4 yz = blend22(s1, s2);
        
        vec4 mm = shuffle4(xy, yz, 0b0101);
        zx = add4(blend0101(zx, yz), blend0101(xy, zx));
        storeu4(mul, add4(mm, zx));
        
        src += 12;
        dir += 12;
        mul += 4;
    }
    while ( mul < end+2 )
    {
        vec4 s0 = mul4(load4(dir  ), sub4(loadu4(src+3), loadu4(src)));
        vec2 s1 = mul2(load2(dir+4), sub2(loadu2(src+7), loadu2(src+4)));

        vec2 xy = getlo(s0);
        vec2 zx = gethi(s0);

        vec2 mm = catshift(xy, s1);
        zx = add2(blend11(zx, s1), blend11(xy, zx));
        storeu2(mul, add2(mm, zx));

        src += 6;
        dir += 6;
        mul += 2;
    }
    while ( mul < end+3 )
    {
        vec2 x = mul2(load2(dir), sub2(loadu2(src+3), loadu2(src)));
        vec2 z = mul1(load2(dir+2), sub1(load1(src+5), load1(src+2)));
        vec2 y = permute2(x);
        
        store1(mul, add1(add1(x, y), z));
        
        src += 3;
        dir += 3;
        ++mul;
    }
    assert_true(mul==end+3);
}


/*
 Ugly piece of code to harvest AVX power...
 FJN @ Strasbourg, 18 and 19.04.2020, improved 22.08.2021
 */
void projectForcesD3D_AVX(index_t nbs, const double* dir,
                          const double* src, const double* mul, double* dst)
{
    // we can assume nbs > 0
    const double* const end = mul - 3 + nbs;
    /*
     This is where the bulk of the work is done, handing 12 scalars per pass.
     The loop is executed if nbs >= 4
     */
    vec4 d2 = setzero4(), m3 = setzero4();
    while ( mul < end )
    {
        vec4 m1 = broadcast1(mul  );
        vec4 m0 = blend31(m3, m1);  // <- m3 = broadcast1(mul-1)
        vec4 m2 = broadcast1(mul+1);
        vec4 p0 = blend31(m1, m2);
        vec4 p2 = broadcast1(mul+2);
        m1 = blend22(m1, m2);
        vec4 p1 = blend22(m2, p2);
        m2 = blend13(m2, p2);
        m3 = broadcast1(mul+3);
        p2 = blend13(p2, m3);
        
        mul += 4;
        vec4 dA = d2;
        vec4 d0 = load4(dir  );
        vec4 d1 = load4(dir+4);
        d2 = load4(dir+8);
        
        vec4 a0 = fmadd4(p0, d0, loadu4(src  ));
        vec4 a1 = fmadd4(p1, d1, loadu4(src+4));
        vec4 a2 = fmadd4(p2, d2, loadu4(src+8));
        
        storeu4(dst  , fnmadd4(m0, catshift1(dA, d0), a0));
        storeu4(dst+4, fnmadd4(m1, catshift1(d0, d1), a1));
        storeu4(dst+8, fnmadd4(m2, catshift1(d1, d2), a2));
        
        dir += 12;
        dst += 12;
        src += 12;
    }
    /*
     We need to consider here multiple cases depending on how many vectors are
     left, since the positive terms are not present on the last vector.
     In addition the above code was not executed if there was 4 or less vectors,
     and the first vector is still a special case.
     */
    switch ( mul - end )
    {
        case 0: {
            assert_true( mul == end );
            // 4 vectors remaining
            vec4 m1 = broadcast1(mul);
            vec4 m0 = blend31(m3, m1);  // m3 from previous round
            vec4 m2 = broadcast1(mul+1);
            vec4 p0 = blend31(m1, m2);
            vec4 p2 = broadcast1(mul+2);
            m1 = blend22(m1, m2);
            vec4 p1 = blend22(m2, p2);
            m2 = blend13(m2, p2);
            p2 = blend13(p2, setzero4()); // { P 0 0 0 }
            vec4 dA = d2;
            vec4 d0 = load4(dir);
            vec4 d1 = load4(dir+4);
            d2 = blend13(cast4(load2(dir+8)), setzero4()); // loading crap, only [0] used
            vec4 a0 = fmadd4(p0, d0, loadu4(src  ));
            vec4 a1 = fmadd4(p1, d1, loadu4(src+4));
            vec4 a2 = fmadd4(p2, d2, loadu4(src+8));
            storeu4(dst  , fnmadd4(m0, catshift1(dA, d0), a0));
            storeu4(dst+4, fnmadd4(m1, catshift1(d0, d1), a1));
            storeu4(dst+8, fnmadd4(m2, catshift1(d1, d2), a2));
            //dir += 12; dst += 12; src += 12;
        } break;
        case 1: {
            assert_true( mul == end+1 );
            // 3 vectors remaining
            vec4 m1 = broadcast1(mul);
            vec4 m0 = blend31(m3, m1);
            vec4 m2 = broadcast1(mul+1);
            vec4 p0 = blend31(m1, m2);
            m1 = blend22(m1, m2);
            vec4 p1 = blend22(m2, setzero4()); // { M M 0 0 }
            vec4 d0 = load4(dir);
            vec2 d1 = load2(dir+4);
            vec4 a0 = fmadd4(p0, d0, loadu4(src));
            vec4 a1 = fmadd4(p1, concatenate22(d1, setzero2()), loadu4(src+4));
            storeu4(dst  , fnmadd4(m0, catshift1(d2, d0), a0));
            storeu4(dst+4, fnmadd4(m1, catshift1(d0, cast4(d1)), a1));
             store1(dst+8, fnmadd2(getlo(m2), unpackhi2(d1,d1), load1(src+8)));
            //dir += 9; dst += 9; src += 9;
        } break;
        case 2: {
            assert_true( mul == end+2 );
            // 2 vectors remaining
            vec4 m1 = broadcast1(mul);
            vec4 m0 = blend31(m3, m1);
            vec4 p0 = blend31(m1, setzero4()); // { P 0 0 0 }
            vec4 d0 = clear4th(load4(dir));  // loading crap; load3(dir);
            vec4 a0 = fmadd4(p0, d0, loadu4(src));
            storeu4(dst  , fnmadd4(m0, catshift1(d2, d0), a0));
            storeu2(dst+4, fnmadd2(getlo(m1), getlo(catshift1(d0, d0)), loadu2(src+4)));
            //dir += 6; dst += 6; src += 6;
        } break;
        case 3: {
            assert_true( mul == end+3 );
            // 1 vector remaining
            vec4 aa = catshift1(d2, setzero4());
            //store3(dst, fnmadd4(mm, aa, load3(src)));
            storeu2(dst, fnmadd2(getlo(m3), getlo(aa), loadu2(src)));
            store1(dst+2, fnmadd1(getlo(m3), gethi(aa), load1(src+2)));
            //dir += 3; dst += 3; src += 3;
        } break;
        default:
            puts("ERROR: unexpected case in projectForcesD3D_AVX!");
    }
}
#endif

#if defined(__SSE3__)
/*
 Ugly piece of code to harvest SSE power for single precision
 FJN @ Strasbourg, 23.08.2021
 */
void projectForcesD3D_SSE(index_t nbs, const float* dir,
                          const float* src, const float* mul, float* dst)
{
    // we can assume nbs > 0
    const float* const end = mul - 3 + nbs;
    /*
     This is where the bulk of the work is done, handing 12 scalars per pass.
     The loop is executed if nbs >= 4
     */
    vec4f d2 = setzero4f(), m3 = setzero4f();
    while ( mul < end )
    {
        vec4f m1 = broadcast1f(mul  );
        vec4f m0 = blend31f(m3, m1);  // <- m3 = broadcast1(mul-1)
        vec4f m2 = broadcast1f(mul+1);
        vec4f p0 = blend31f(m1, m2);
        vec4f p2 = broadcast1f(mul+2);
        m1 = blend22f(m1, m2);
        vec4f p1 = blend22f(m2, p2);
        m2 = blend13f(m2, p2);
        m3 = broadcast1f(mul+3);
        p2 = blend13f(p2, m3);
        
        mul += 4;
        vec4f dA = d2;
        vec4f d0 = load4f(dir  );
        vec4f d1 = load4f(dir+4);
        d2 = load4f(dir+8);
        
        vec4f a0 = fmadd4f(p0, d0, loadu4f(src  ));
        vec4f a1 = fmadd4f(p1, d1, loadu4f(src+4));
        vec4f a2 = fmadd4f(p2, d2, loadu4f(src+8));
        
        storeu4f(dst  , fnmadd4f(m0, catshift1f(dA, d0), a0));
        storeu4f(dst+4, fnmadd4f(m1, catshift1f(d0, d1), a1));
        storeu4f(dst+8, fnmadd4f(m2, catshift1f(d1, d2), a2));
        
        dir += 12;
        dst += 12;
        src += 12;
    }
    /*
     We need to consider here multiple cases depending on how many vectors are
     left, since the positive terms are not present on the last vector.
     In addition the above code was not executed if there was 4 or less vectors,
     and the first vector is still a special case.
     */
    switch ( mul - end )
    {
        case 0: {
            assert_true( mul == end );
            // 4 vectors remaining
            vec4f m1 = broadcast1f(mul  );
            vec4f m0 = blend31f(m3, m1);  // m3 from previous round
            vec4f m2 = broadcast1f(mul+1);
            vec4f p0 = blend31f(m1, m2);
            vec4f p2 = broadcast1f(mul+2);
            m1 = blend22f(m1, m2);
            vec4f p1 = blend22f(m2, p2);
            m2 = blend13f(m2, p2);
            p2 = blend13f(p2, setzero4f()); // { P 0 0 0 }
            vec4f dA = d2;
            vec4f d0 = load4f(dir  );
            vec4f d1 = load4f(dir+4);
            d2 = load1f(dir+8);
            vec4f a0 = fmadd4f(p0, d0, loadu4f(src  ));
            vec4f a1 = fmadd4f(p1, d1, loadu4f(src+4));
            vec4f a2 = fmadd4f(p2, d2, loadu4f(src+8));
            storeu4f(dst  , fnmadd4f(m0, catshift1f(dA, d0), a0));
            storeu4f(dst+4, fnmadd4f(m1, catshift1f(d0, d1), a1));
            storeu4f(dst+8, fnmadd4f(m2, catshift1f(d1, d2), a2));
            //dir += 12; dst += 12; src += 12;
        } break;
        case 1: {
            assert_true( mul == end+1 );
            // 3 vectors remaining
            vec4f m1 = broadcast1f(mul  );
            vec4f m0 = blend31f(m3, m1);
            vec4f m2 = broadcast1f(mul+1);
            vec4f p0 = blend31f(m1, m2);
            m1 = blend22f(m1, m2);
            vec4f p1 = blend22f(m2, setzero4f()); // { M M 0 0 }
            vec4f d0 = load4f(dir  );
            vec4f d1 = load2f(dir+4);
            vec4f a0 = fmadd4f(p0, d0, loadu4f(src  ));
            vec4f a1 = fmadd4f(p1, blend22f(d1, setzero4f()), loadu4f(src+4));
            storeu4f(dst  , fnmadd4f(m0, catshift1f(d2, d0), a0));
            storeu4f(dst+4, fnmadd4f(m1, catshift1f(d0, d1), a1));
             store1f(dst+8, fnmadd2f(m2, catshift1f(d1, d1), load1f(src+8)));
            //dir += 9; dst += 9; src += 9;
        } break;
        case 2: {
            assert_true( mul == end+2 );
            // 2 vectors remaining
            vec4f m1 = broadcast1f(mul);
            vec4f m0 = blend31f(m3, m1);
            vec4f p0 = blend31f(m1, setzero4f()); // { P 0 0 0 }
            vec4f d0 = clear4th(load4f(dir));  // loading crap; load3(dir);
            vec4f a0 = fmadd4f(p0, d0, loadu4f(src));
            storeu4f(dst, fnmadd4f(m0, catshift1f(d2, d0), a0));
            store2f(dst+4, fnmadd2f(m1, catshift1f(d0, d0), load2f(src+4)));
            //dir += 6; dst += 6; src += 6;
        } break;
        case 3: {
            assert_true( mul == end+3 );
            // 1 vector remaining
            vec4f aa = catshift1f(d2, setzero4f());
            //store3(dst, fnmadd4(mm, aa, load3(src)));
            store2f(dst, fnmadd2f(m3, aa, load2f(src)));
            store1f(dst+2, fnmadd2f(m3, gethi2f(aa), load1f(src+2)));
            //dir += 3; dst += 3; src += 3;
        } break;
        default:
            puts("ERROR: unexpected case in projectForcesD3D_AVX!");
    }
}
#endif

//------------------------------------------------------------------------------
#pragma mark - Correction terms to the Projection


///expanded implementation:
void addProjectionDiff_R(const index_t nbs, const real* mul, const real* X, real* Y)
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


/// scalar implementation
void addProjectionDiff_F(const index_t nbs, const real* mul, const real* X, real* Y)
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

#if ( DIM == 2 ) && REAL_IS_DOUBLE && USE_SIMD

void addProjectionDiff2D_SSE(const index_t nbs, const double* mul, const double* src, double* dst)
{
    vec2 p = load2(src);
    vec2 n = load2(dst);
    
    const double* end = mul + nbs;
    src += DIM;
    
    while ( mul < end )
    {
        vec2 x = load2(src);
        src += DIM;
        vec2 w = mul2(loaddup2(mul), sub2(x, p));
        p = x;
        store2(dst, add2(n, w));
        dst += DIM;
        n = sub2(load2(dst), w);
        ++mul;
    }
    store2(dst, n);
}

#endif

#if ( DIM == 2 ) && REAL_IS_DOUBLE && defined(__AVX__)

void addProjectionDiff_AVX(const index_t nbs, const double* mul, const double* X, double* Y)
{
    double * pY = Y;
    double const* pX = X;
    double const* pM = mul;
    
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
    
    double const*const end = mul + nbs;
    while ( pM < end )
    {
        vec4 a = broadcast2(pM);
        vec4 m = duplohi4(a);

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

#if ( DIM == 3 ) && REAL_IS_DOUBLE && USE_SIMD

/**
 Idea of the calculation:
    const real w = mul[i] * ( X[d+DIM*i+DIM] - X[d+DIM*i] );
    Y[d+DIM*i    ] += w;
    Y[d+DIM*i+DIM] -= w;
 */
void addProjectionDiff3D_SSE(const index_t nbs, const double* mul, const double* src, double* dst)
{
    // there should be at least 2 3D vectors, and these loads are safe:
    vec2 p0 = loadu2(src);
    vec2 p1 = load1(src+2); // upper value not used
    vec2 n0 = loadu2(dst);
    vec2 n1 = loadu2(dst+2);

    const double* end = mul + nbs - 1;
    src += DIM;
    
#if 0
    p1 = catshift(p0, p1);
    p0 = unpacklo2(p0, p0); // lower value will not be used
    while ( mul < end )
    {
        // process 2 3D vectors = 6 scalars
        vec2 x0 = loadu2(src);
        vec2 x1 = loadu2(src+2);
        vec2 x2 = loadu2(src+4);
        src += DIM * 2;
        vec2 m = loadu2(mul);
        vec2 m0 = unpacklo2(m, m);
        vec2 m1 = unpackhi2(m, m);
        vec2 w0 = mul2(sub2(x0, catshift(p0, p1)), m0);
        vec2 w1 = mul2(sub2(x1, catshift(p1, x0)), blend11(m0, m1));
        vec2 w2 = mul2(sub2(x2, catshift(x0, x1)), m1);
        p0 = x1;
        p1 = x2;
        x1 = sub2(n1, catshift(setzero2(), w0));
        x2 = sub2(loadu2(dst+4), catshift(w0, w1));
        storeu2(dst  , add2(w0, n0));
        storeu2(dst+2, add2(w1, x1));
        storeu2(dst+4, add2(w2, x2));
        dst += DIM * 2;
        n0 = sub2(loadu2(dst), catshift(w1, w2));
        n1 = sub2(loadu2(dst+2), catshift(w2, setzero2()));
        mul += 2;
    }
    p0 = catshift(p0, p1);
    p1 = unpackhi2(p1, p1); // upper value not used
#endif
    
    while ( mul <= end )
    {
        vec2 x0 = loadu2(src);
        vec2 x1 = load1(src+2);
        src += DIM;
        vec2 m = loaddup2(mul);
        vec2 w0 = mul2(m, sub2(x0, p0));
        vec2 w1 = mul1(m, sub1(x1, p1));
        p0 = x0;
        p1 = x1;
        storeu2(dst, add2(n0, w0));
        store1(dst+2, add1(n1, w1));
        dst += DIM;
        n0 = sub2(loadu2(dst), w0);
        n1 = sub1(load1(dst+2), w1);
        ++mul;
    }
    storeu2(dst, n0);
    store1(dst+2, n1);
}

#endif
