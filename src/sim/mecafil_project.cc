// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "exceptions.h"
#include "vecprint.h"

/**
 This is a C-translation of the LAPACK reference implementation of dpttrf()
 (that is Thomas' algorithm to factorize a tridiagonal matrix)
*/
void lapack_xpttrf(int N, real* D, real* E, int* INFO)
{
    for ( int i = 0; i < N-1; ++i )
    {
        if ( D[i] < 0 )
        {
            *INFO = i;
            return;
        }
        real e = E[i];
        E[i] = e / D[i];
        D[i+1] = D[i+1] - e * E[i];
    }
}

/**
 This is a C-translation of the LAPACK reference implementation of dptts2()
 
 *     Solve A * X = B using the factorization A = L*D*L**T,
 *     overwriting each right hand side vector with its solution.
 
     DO I = 2, N
         B( I ) = B( I ) - B( I-1 ) * E( I-1 )
     CONTINUE
 
     B( N ) = B( N ) / D( N )
 
     DO I = N - 1, 1, -1
         B( I ) = B( I ) / D( I ) - B( I+1 ) * E( I )
     CONTINUE
 */
void lapack_xptts2(int N, int NRHS, const real* D, const real* E, real * B, int LDB)
{
    assert_true( NRHS == 1 ); // in this case, LDB is not used

    for ( int i = 1; i < N; ++i )
        B[i] = B[i] - B[i-1] * E[i-1];
    
    B[N-1] = B[N-1] / D[N-1];
    
    for ( int i = N-2; i >= 0; --i )
        B[i] = B[i] / D[i] - B[i+1] * E[i];
}


/**
 Italian factorization that uses different operations
 From
     Numerical Mathematics
     Springer (2000) ISBN 0-387-98959-5
     A. Quarteroni, R. Sacco and F. Saleri
     Page 93
 
 PSEUDO code with indices starting at 1:
     b=[0;b];
     c=[c;0];
     gamma(1) = 1.0 / a(1);
     for i = 2:N
         gamma(i) = 1.0 / ( a(i) - b(i) * gamma(i-1) * c(i-1) );
 */
void italian_xpttrf(int size, real* D, real* E, int* INFO)
{
    D[0] = 1.0 / D[0];
    
    for ( int n = 1; n < size; ++n )
        D[n] = 1.0 / ( D[n] - E[n-1] * E[n-1] * D[n-1] );
}


/**
 Italian version without divisions
 
 PSEUDO code with indices starting at 1:
     y(1) = gamma(1) * f(1);
     for i = 2:N
         y(i) = gamma(i) * ( f(i) - b(i) * y(i-1) );
     x(N) = y(N);
     for i = N-1:-1:1
         x(i) = y(i) - gamma(i) * c(i) * x(i+1);
 */
void italian_xptts2(int size, int nrhs, real const* D, real const* E, real * B, int LDB)
{
    assert_true(nrhs == 1); // in this case, LDB is not used
 
    B[0] = D[0] * B[0];
    
    // upward recursion on B[]
    for ( int n = 1; n < size; ++n )
        B[n] = D[n] * ( B[n] - B[n-1] * E[n-1] );
    
    // downward recursion on B[]
    for ( int n = size-2; n >= 0; --n )
        B[n] = B[n] - D[n] * E[n] * B[n+1];
}


/**
 Based on the 'Italian' version, precalculating constant products:
 
     DEL[n] <-  D[n] * E[n-1]
     E[n]   <-  D[n] * E[n]
 */
void alsatian_xpttrf(int size, real* D, real* E, int* INFO)
{
    real * DEL = D + size;
    D[0] = 1.0 / D[0];
    
    real x = D[0];
    for ( int n = 1; n < size; ++n )
    {
        //D[n] = 1.0 / ( D[n] - E[n-1] * E[n-1] * D[n-1] );
        x = 1.0 / ( D[n] - E[n-1] * E[n-1] * x );
        D[n] = x;
        // precalculate product that is constant:
        DEL[n] = x * E[n-1];
    }
    
    // precalculate two products that are constant:
    for ( int n = 0; n < size-1; ++n )
        E[n] = D[n] * E[n];
}


/**
 Based on the 'Italian' version, using precalculated constant terms
 */
void alsatian_xptts2(int size, int nrhs, real const* D, real const* E, real * B, int LDB)
{
    assert_true(nrhs == 1); // in this case, LDB is not used
    real const* DEL = D + size;

    B[0] = D[0] * B[0];
    
    // upward recursion on B[]
    for ( int n = 1; n < size; ++n )
    {
        // B[n] = D[n] * ( B[n] - B[n-1] * E[n-1] );
        B[n] = D[n] * B[n] - DEL[n] * B[n-1];
    }
    
    // downward recursion on B[]
    for ( int n = size-2; n >= 0; --n )
    {
        // B[n] = B[n] - D[n] * E[n] * B[n+1];
        B[n] = B[n] - E[n] * B[n+1];
    }
}


/**
Enable the optimized version of setSpeedsFromForces, which only works in 2D.
It is built on the custom factorization 'alsatian_xpttrf'
should be defined as ZERO
*/
#define NEW_CUSTOM_XPTTRF 0 //( DIM == 2 )

/*
 Selection of LAPACK routines
 The LAPACK implementation is the safest choice
 The Alsatian version may be faster as it avoids divisions
 */

#if NEW_CUSTOM_XPTTRF
#  define DPTTRF alsatian_xpttrf
#  define DPTTS2 alsatian_xptts2
#elif ( 1 )
#  define DPTTRF alsatian_xpttrf
#  define DPTTS2 alsatian_xptts2
#elif ( 0 )
#  define DPTTRF italian_xpttrf
#  define DPTTS2 italian_xptts2
#elif ( 0 )
#  define DPTTRF lapack_xpttrf
#  define DPTTS2 lapack_xptts2
#else
#  define DPTTRF lapack::xpttrf
#  define DPTTS2 lapack::xptts2
#endif

//------------------------------------------------------------------------------
#pragma mark -


void Mecafil::buildProjection()
{
    //reset all variables for the projections:
    rfAllocated  = 0;
    mtJJt        = nullptr;
    mtJJtiJforce = nullptr;
}


void Mecafil::allocateProjection(const size_t nbp)
{
    if ( rfAllocated < nbp )
    {
        //std::clog << reference() << "allocateProjection(" << nbp << ")\n";
        if ( mtJJt )
            free_real(mtJJt);
        
        // make a multiple of chunk to align memory:
        rfAllocated  = chunk_real(nbp);
        assert_true(rfAllocated > 0);
        
        real * mem = new_real(4*rfAllocated);
        //zero_real(4*rfAllocated, mem);
        
        mtJJt        = mem;
        mtJJtU       = mem + rfAllocated * 2;
        mtJJtiJforce = mem + rfAllocated * 3;
    }
}


void Mecafil::destroyProjection()
{
    //std::clog << reference() << "destroyProjection\n";
    if ( mtJJt ) free_real(mtJJt);
    mtJJt        = nullptr;
    mtJJtU       = nullptr;
    mtJJtiJforce = nullptr;
}


#if NEW_ANISOTROPIC_FIBER_DRAG

/// version for drag coefficients doubled in directions orthogonal to fiber
void Mecafil::makeProjection()
{
    assert_true( nbPoints() >= 2 );
    assert_true( rfAllocated >= nbPoints() );

    //set the diagonal and off-diagonal of J*J'
    const unsigned nbu = nbPoints() - 2;
    const real*const dif = rfDiff;
   real b = 1;

    for ( unsigned jj = 0; jj < nbu; ++jj )
    {
        const real* X = dif + DIM * jj;
#if ( DIM == 2 )
        real xn = X[0]*X[2] + X[1]*X[3];
#else
        real xn = X[0]*X[3] + X[1]*X[4] + X[2]*X[5];
#endif
        real a = 0.25 * ( 1 + xn ) * ( 1 + xn );
        mtJJt[jj]  = 2.0 + a + b;
        mtJJtU[jj] = -xn - a;
        b = a;
    }
    
    mtJJt[nbu] = 3.0 + b;
    
    int info = 0;
    DPTTRF(nbu+1, mtJJt, mtJJtU, &info);

    if ( 0 )
    {
        std::clog << "D "; VecPrint::print(std::clog, nbu+1, mtJJt, 3);
        std::clog << "E "; VecPrint::print(std::clog, nbu, mtJJtU, 3);
        //std::clog << "X="; VecPrint::print(std::clog, DIM*(nbu+2), pPos);
    }

    if ( info )
    {
        std::clog << "Mecafil::makeProjection failed (info = " << info << ")\n";
        throw Exception("could not build Fiber's projection matrix");
    }
}

#else

/// standard version with isotropic drag coefficient
void Mecafil::makeProjection()
{
    assert_true( nbPoints() >= 2 );
    assert_true( rfAllocated >= nbPoints() );

    //set the diagonal and off-diagonal of J*J'
    const unsigned nbu = nbPoints() - 2;
    const real*const dif = rfDiff;

    for ( unsigned jj = 0; jj < nbu; ++jj )
    {
        const real* X = dif + DIM * jj;
#if ( DIM == 2 )
        real xn = X[0]*X[2] + X[1]*X[3];
#else
        real xn = X[0]*X[3] + X[1]*X[4] + X[2]*X[5];
#endif
        
#if ( DIM == 2 )
        mtJJt[jj] = 2.0 * ( X[0]*X[0] + X[1]*X[1] );
#else
        mtJJt[jj] = 2.0 * ( X[0]*X[0] + X[1]*X[1] + X[2]*X[2] );
#endif
        // the diagonal term should be nearly equal to 2, since dif[] vectors are normalized
        //mtJJt[jj]  = 2.0;
        mtJJtU[jj] = -xn;
    }
    
    const real* X = dif + DIM*nbu;
#if ( DIM == 2 )
    mtJJt[nbu] = 2.0 * ( X[0]*X[0] + X[1]*X[1] );
#else
    mtJJt[nbu] = 2.0 * ( X[0]*X[0] + X[1]*X[1] + X[2]*X[2] );
#endif
    // the diagonal term should be nearly equal to 2, since dif[] vectors are normalized
    //mtJJt[nbu] = 2.0;

    int info = 0;
    DPTTRF(nbu+1, mtJJt, mtJJtU, &info);

    if ( 0 )
    {
        std::clog << "D "; VecPrint::print(std::clog, nbu+1, mtJJt, 3);
        std::clog << "E "; VecPrint::print(std::clog, nbu, mtJJtU, 3);
        //std::clog << "X="; VecPrint::print(std::clog, DIM*(nbu+2), pPos);
    }

    if ( info )
    {
        std::clog << "Mecafil::makeProjection failed (info = " << info << ")\n";
        throw Exception("could not build Fiber's projection matrix");
    }
}

#endif
//------------------------------------------------------------------------------
#pragma mark -

/**
 Perform first calculation needed by projectForces:
 tmp <- J * X
 */
void projectForcesU_(unsigned nbs, const real* dif, const real* vec, real* mul)
{
    #pragma vector unaligned
    for ( unsigned jj = 0; jj < nbs; ++jj )
    {
        const real * X = vec + DIM * jj;
        const real * d = dif + DIM * jj;
        mul[jj] = d[0] * ( X[DIM  ] - X[0] )
                + d[1] * ( X[DIM+1] - X[1] )
#if ( DIM > 2 )
                + d[2] * ( X[DIM+2] - X[2] )
#endif
        ;
    }
}

/**
 Perform second calculation needed by projectForces:
 Y <- s * ( X + Jt * tmp )
 */
void projectForcesD_(unsigned nbs, const real* dif, const real alpha, const real* X, const real* mul, real* Y)
{
    for ( unsigned d = 0, e = DIM*nbs; d < DIM; ++d, ++e )
    {
        Y[d] = alpha * ( X[d] + dif[d    ] * mul[    0] );
        Y[e] = alpha * ( X[e] - dif[e-DIM] * mul[nbs-1] );
    }
    
    for ( unsigned jj = 1; jj < nbs; ++jj )
    {
        const unsigned kk = DIM*jj;
        Y[kk  ] = alpha * ( X[kk  ] + dif[kk  ] * mul[jj] - dif[kk-DIM  ] * mul[jj-1] );
        Y[kk+1] = alpha * ( X[kk+1] + dif[kk+1] * mul[jj] - dif[kk-DIM+1] * mul[jj-1] );
#if ( DIM > 2 )
        Y[kk+2] = alpha * ( X[kk+2] + dif[kk+2] * mul[jj] - dif[kk-DIM+2] * mul[jj-1] );
#endif
    }
}

/**
 Perform second calculation needed by projectForces:
 */
void projectForcesD__(unsigned nbs, const real* dif, const real alpha,
                      const real* X, const real* mul, real* Y)
{
    real a0 = dif[0] * mul[0];
    real a1 = dif[1] * mul[0];
#if ( DIM > 2 )
    real a2 = dif[2] * mul[0];
#endif
    
    Y[0] = alpha * ( X[0] + a0 );
    Y[1] = alpha * ( X[1] + a1 );
#if ( DIM > 2 )
    Y[2] = alpha * ( X[2] + a2 );
#endif
    
    for ( unsigned jj = 1; jj < nbs; ++jj )
    {
        const unsigned kk = DIM * jj;
        real b0 = dif[kk  ] * mul[jj];
        Y[kk  ] = alpha * ( X[kk  ] + b0 - a0 );
        a0 = b0;
        
        real b1 = dif[kk+1] * mul[jj];
        Y[kk+1] = alpha * ( X[kk+1] + b1 - a1 );
        a1 = b1;
        
#if ( DIM > 2 )
        real b2 = dif[kk+2] * mul[jj];
        Y[kk+2] = alpha * ( X[kk+2] + b2 - a2 );
        a2 = b2;
#endif
    }
    
    const unsigned ee = DIM * nbs;
    Y[ee  ] = alpha * ( X[ee  ] - a0 );
    Y[ee+1] = alpha * ( X[ee+1] - a1 );
#if ( DIM > 2 )
    Y[ee+2] = alpha * ( X[ee+2] - a2 );
#endif
}


void projectForcesD___(unsigned nbs, const real* dif, const real alpha,
                       const real* X, const real* mul, real* Y)
{
    // set Y, using values in X and lag
    //projectForcesD(nbs, rfDiff, alpha, X, mul, Y);
    long nn = DIM * nbs;
    long pp = nn - DIM;
    
    real d0 = dif[pp  ] * mul[nbs-1];
    real d1 = dif[pp+1] * mul[nbs-1];
    
    Y[nn  ] = alpha * ( X[nn  ] - d0 );
    Y[nn+1] = alpha * ( X[nn+1] - d1 );
    
    real p0 = X[pp  ] - d0;
    real p1 = X[pp+1] - d1;
    
    for ( long n = nbs-1; n > 0; --n )
    {
        nn = DIM * n;
        pp = nn - DIM;
        
        d0 = dif[pp  ] * mul[n-1];
        d1 = dif[pp+1] * mul[n-1];

        Y[nn  ] = alpha * ( p0 - d0 );
        Y[nn+1] = alpha * ( p1 - d1 );

        p0 = X[pp  ] + d0;
        p1 = X[pp+1] + d1;
    }

    Y[0] = alpha * p0;
    Y[1] = alpha * p1;
}

#if ( DIM == 2 ) && defined(__SSE3__) && REAL_IS_DOUBLE

#include "simd.h"

/**
 Perform first calculation needed by projectForces:
 */
inline void projectForcesU_SSE(unsigned nbs, const real* dif, const real* X, real* mul)
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
        storeu2(pT, hadd2(a, b));
        pT += 2;
    }
    
    if ( pT < end+2 )
    {
        y = load2(pX+2);
        vec2 a = mul2(sub2(y, x), load2(pD));
        storelo(pT, hadd2(a, a));
    }
}

/**
 Perform second calculation needed by projectForces:
 */
inline void projectForcesD_SSE(unsigned nbs, const real* dif, const real alpha,
                               const real* X, const real* mul, real* Y)
{
    real *pY = Y;
    const real* pX = X;
    const real* pD = dif;
    
    vec2 cc = load2(X);
    vec2 ss = set2(alpha);
    
    real const* pM = mul;
    real const*const end = mul + nbs;
    while ( pM < end )
    {
        pX += DIM;
        vec2 d = mul2(load2(pD), loaddup2(pM));
        ++pM;
        pD += DIM;
        store2(pY, mul2(ss, add2(cc, d)));
        pY += DIM;
        cc = sub2(load2(pX), d);
    }
    store2(pY, mul2(ss, cc));
}

#endif

#if ( DIM == 2 ) && defined(__AVX__) && REAL_IS_DOUBLE

#include "simd.h"

/**
 Perform first calculation needed by projectForces:
 tmp <- J * X
 F. Nedelec, 9.12.2016, 6.9.2018
 */
inline void projectForcesU_AVX(unsigned nbs, const real* dif, const real* X, real* mul)
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
        vec4 p = permute2f128(a,b,0x20), q = permute2f128(a,b,0x31);
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
        storelo(pT, hadd2(a, a));
    }
}


/**
 Perform second calculation needed by projectForces:
 Y <- s * ( X + Jt * tmp )
 ATTENTION: memory X and Y are not necessarily aligned since they are chunck from
 an array containing contiguous coordinates
 F. Nedelec, 9.12.2016, 23.03.2018
 */
inline void projectForcesD_AVX(unsigned nbs, const real* dif, const real alpha,
                               const real* X, const real* mul, real* Y)
{
    real *pY = Y;
    const real* pX = X;
    const real* pD = dif;
    
    vec4 cc = setzero4();
    vec4 ss = set4(alpha);
    
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
        storeu4(pY, mul4(ss, z));
        pY += 4;
    }
    
    vec2 c = gethi(cc);
    
    if ( odd )
    {
        assert( pM + 1 == mul + nbs );
        vec2 m = loaddup2(pM);
        vec2 x = mul2(m, load2(pD));
        vec2 z = add2(load2(pX), sub2(x, c));
        storeu2(pY, mul2(set2(alpha), z));
        c = x;
        pY += 2;
        pX += 2;
    }
    
    vec2 z = sub2(load2(pX), c);
    storeu2(pY, mul2(set2(alpha), z));
    assert( pY == Y + DIM * nbs );
    assert( pX == X + DIM * nbs );
}

#endif


/**
 Perform first calculation needed by projectForces:
 */
inline void projectForcesU_PTR(unsigned nbs, const real* dif, const real* X, real* mul)
{
    const real * pX = X + DIM;
    const real * pM = dif;
    real x3, x0 = X[0];
    real x4, x1 = X[1];
#if ( DIM >= 3 )
    real x5, x2 = X[2];
#endif
    real *const end = mul + nbs;
    
    //normally optimized version
    for ( real* pT = mul; pT < end; ++pT )
    {
        x3 = pX[0];
        x4 = pX[1];
#if ( DIM == 2 )
        pT[0] = pM[0] * (x3 - x0) + pM[1] * (x4 - x1);
#elif ( DIM >= 3 )
        x5 = pX[2];
        pT[0] = pM[0] * (x3 - x0) + pM[1] * (x4 - x1) + pM[2] * (x5 - x2);
        x2 = x5;
#endif
        pX += DIM;
        pM += DIM;
        x0 = x3;
        x1 = x4;
    }
}


/**
 Perform first calculation needed by projectForces:
 */
inline void projectForcesU_PTR2(unsigned nbs, const real* dif, const real* X, real* mul)
{
    const real * pX = X + DIM;
    const real * pM = dif;
    real x3, x0 = X[0];
    real x4, x1 = X[1];
#if ( DIM >= 3 )
    real x5, x2 = X[2];
#endif
    real *const end = mul + nbs;
    
    //further optimization with manual loop-unrolling
    real* pT = mul;
    if ( nbs & 1 )
    {
        x3 = pX[0];
        x4 = pX[1];
#if ( DIM == 2 )
        pT[0] = pM[0] * (x3 - x0) + pM[1] * (x4 - x1);
#elif ( DIM >= 3 )
        x5 = pX[2];
        pT[0] = pM[0] * (x3 - x0) + pM[1] * (x4 - x1) + pM[2] * (x5 - x2);
        x2 = x5;
#endif
        ++pT;
        pX += DIM;
        pM += DIM;
        x0 = x3;
        x1 = x4;
    }
    
    while ( pT < end )
    {
        x3 = pX[0];
        x4 = pX[1];
#if ( DIM == 2 )
        pT[0] = pM[0] * (x3 - x0) + pM[1] * (x4 - x1);
#elif ( DIM >= 3 )
        x5 = pX[2];
        pT[0] = pM[0] * (x3 - x0) + pM[1] * (x4 - x1) + pM[2] * (x5 - x2);
#endif
        
#if ( DIM == 2 )
        x0 = pX[2];
        x1 = pX[3];
        pT[1] = pM[2] * (x0 - x3) + pM[3] * (x1 - x4);
#elif ( DIM >= 3 )
        x0 = pX[3];
        x1 = pX[4];
        x2 = pX[5];
        pT[1] = pM[3] * (x0 - x3) + pM[4] * (x1 - x4) + pM[5] * (x2 - x5);
#endif
        
        pT += 2;
        pX += 2*DIM;
        pM += 2*DIM;
    }
    assert_true( pT == end );
}

/**
 Perform second calculation needed by projectForces:
 */
void projectForcesD_PTR(unsigned nbs, const real* dif, const real alpha,
                        const real* X, const real* mul, real* Y)
{
    // Y <- X + Jt * tmp :
    real x0 = X[0];
    real x1 = X[1];
#if ( DIM > 2 )
    real x2 = X[2];
#endif
    
    const real* pX = X+DIM;
    const real* pM = dif;
    real *pY = Y;
    real const*const end = mul+nbs;
    for ( real const* pT = mul; pT < end; ++pT )
    {
        real y0 = *pT * pM[0];
        real y1 = *pT * pM[1];
#if ( DIM > 2 )
        real y2 = *pT * pM[2];
#endif
        pM  += DIM;
        pY[0]  = alpha * ( x0 + y0 );
        pY[1]  = alpha * ( x1 + y1 );
#if ( DIM > 2 )
        pY[2]  = alpha * ( x2 + y2 );
#endif
        pY  += DIM;
        x0     = pX[0] - y0;
        x1     = pX[1] - y1;
#if ( DIM > 2 )
        x2     = pX[2] - y2;
#endif
        pX  += DIM;
    }
    pY[0] = alpha * x0;
    pY[1] = alpha * x1;
#if ( DIM > 2 )
    pY[2] = alpha * x2;
#endif
}


//------------------------------------------------------------------------------
#pragma mark -


#if ( DIM == 2 ) && REAL_IS_DOUBLE
#  if defined(__AVX__)
#    warning "Using AVX implementation"
#    define projectForcesU projectForcesU_AVX
#    define projectForcesD projectForcesD_AVX
#  elif defined(__SSE3__)
#    warning "Using SSE3 implementation"
#    define projectForcesU projectForcesU_SSE
#    define projectForcesD projectForcesD_SSE
#  else
#    define projectForcesU projectForcesU_
#    define projectForcesD projectForcesD_
#  endif
#else
#  define projectForcesU projectForcesU_
#  define projectForcesD projectForcesD_
#endif

#if ( NEW_CUSTOM_XPTTRF == 0 )

void Mecafil::setSpeedsFromForces(const real* X, const real alpha, real* Y) const
{
    const unsigned nbs = nbSegments();
    //printf("X  "); VecPrint::print(std::clog, DIM*nbPoints(), X);

    // calculate `lag` without modifying `X`
#if NEW_ANISOTROPIC_FIBER_DRAG
    scaleTangentially(nPoints, X, rfDir, rfVTP);
    projectForcesU(nbs, rfDiff, rfVTP, rfLLG);
#else
    projectForcesU(nbs, rfDiff, X, rfLLG);
#endif

    // lag <- inv( J * Jt ) * lag to find the Lagrange multipliers
    DPTTS2(nbs, 1, mtJJt, mtJJtU, rfLLG, nbs);
    
    // set Y, using values in X and lag
    projectForcesD(nbs, rfDiff, alpha/rfDragPoint, X, rfLLG, Y);

#if NEW_ANISOTROPIC_FIBER_DRAG
    scaleTangentially(nPoints, Y, rfDir, Y);
#endif
    //printf("Y  "); VecPrint::print(std::clog, DIM*nbPoints(), Y);
}

#else

/**
This optimized setSpeedFromForces() works only in 2D,
and only with the 'alsatian' factorization:
*/
void Mecafil::setSpeedsFromForces(const real* X, const real alpha, real* Y) const
{
    const long nbs = nbSegments();
    vec2 LL;
    const real* pF = rfDiff;
    real* pM = rfLLG;
#if NEW_ANISOTROPIC_FIBER_DRAG
    scaleTangentially(nPoints, X, rfDir, rfVTP);
    real const* pX = rfVTP;
#else
    real const* pX = X;
#endif

    {
    //printf("X  "); VecPrint::print(std::clog, DIM*nbPoints(), X);
    real const* D = mtJJt;
    real const* DEL = mtJJt + nbs;

    // calculate `lag` without modifying `X`
    //projectForcesU(nbs, rfDiff, X, mul);
    real L0, L1 = 0.0;
    vec2 LL = setzero2();
    
    real const*const end = rfLLG + nbs;
    
    long n = 0;
    while ( pM <= end )
    {
        vec4 d = mul4(sub4(loadu4(pX+2), loadu4(pX)), load4(pF));
        pX += 4;
        pF += 4;
        vec2 h = gethi(d);
        vec2 a = add2(unpacklo2(getlo(d),h), unpackhi2(getlo(d),h));
        //mul[jj] = dif[0] * ( X[DIM] - X[0] ) + dif[1] * ( X[DIM+1] - X[1] )
        a = mul2(load2(D), a);
        D += 2;
#if 1
        //L0 = D[n  ] * a[0] - DEL[n  ] * L1;
        //L1 = D[n+1] * a[1] - DEL[n+1] * L0;
        L0 = a[0] - DEL[0] * L1;
        L1 = a[1] - DEL[1] * L0;
        pM[0] = L0;
        pM[1] = L1;
#else
        h = load2(DEL);
        // assuming LL contains old value in higher register
        vec2 N = fnmadd1(unpackhi2(LL, LL), h, a);
        LL = fnmadd2(shuffle2(LL, N, 1), h, a);
        store2(pM, LL);
#endif
        DEL += 2;
        pM += 2;
    }
    }


#if ( 1 )
    real const* E = mtJJtU;
    real * mul = rfLLG;
    // downward recursion
    for ( long n = nbs-2; n >= 0; --n )
        mul[n] = mul[n] - E[n] * mul[n+1];

    projectForcesD(nbs, rfDiff, alpha/rfDragPoint, X, mul, Y);
#else
    // set Y, using values in X and lag
    const vec2 aa = set2(alpha/rfDragPoint);

    long e = DIM * nbs - DIM;
    pM = rfLLG + nbs - 1;
    pX = X + e;
    pF = rfDiff + e;
    real * pY = Y + DIM * nbs;

    LL = loaddup2(pM);
    //real d0 = rfDiff[p  ] * mul[nbs-1];
    //real d1 = rfDiff[p+1] * mul[nbs-1];
    vec2 x = mul2(load2(pF), LL);
    pF -= 2;
    pM -= 1;
    
    //Y[n  ] = alpha * ( X[n  ] - d0 );
    //Y[n+1] = alpha * ( X[n+1] - d1 );
    store2(pY, mul2(aa, sub2(load2(pX+DIM), x)));
    pY -= 2;
    
    //real p0 = X[p  ] - d0;
    //real p1 = X[p+1] - d1;
    x = sub2(load2(pX), x);
    pX -= 2;

    real const* pE = mtJJtU + nbs - 2;

    // nbs-1 iterations going downward:
    //for ( long n = nbs-2; n >= 0; --n )
    while ( pX >= X )
    {
    //IACA_START
        // mul[n] = mul[n] - E[n] * mul[n+1];
        // assuming LL contains the old value in the two registers
        //assert_true(pM==rfLLG+n);
        LL = fnmadd2(loaddup2(pE), LL, loaddup2(pM));
        storelo(pM, LL);
        pM -= 1;
        pE -= 1;

        //d0 = rfDiff[p  ] * L;
        //d1 = rfDiff[p+1] * L;
        //assert_true(pF==rfDiff+DIM*n);
        vec2 dd = mul2(load2(pF), LL);
        pF -= 2;
        
        //Y[p+DIM  ] = alpha * ( p0 - d0 );
        //Y[p+DIM+1] = alpha * ( p1 - d1 );
        //assert_true(pY==Y+DIM*n+DIM);
        store2(pY, mul2(aa, sub2(x, dd)));
        pY -= 2;
        
        //p0 = X[p  ] + d0;
        //p1 = X[p+1] + d1;
        //assert_true(pX==X+DIM*n);
        x = add2(load2(pX), dd);
        pX -= 2;
    }
    //IACA_END

    //Y[0] = alpha * p0;
    //Y[1] = alpha * p1;
    store2(Y, mul2(aa, x));
#endif

#if NEW_ANISOTROPIC_FIBER_DRAG
    scaleTangentially(nPoints, Y, rfDir, Y);
#endif
}
#endif


void Mecafil::computeTensions(const real* force)
{
    const unsigned nbs = nbSegments();
    
#if NEW_ANISOTROPIC_FIBER_DRAG
    
    scaleTangentially(nPoints, force, rfDir, rfVTP);
    projectForcesU(nbs, rfDiff, rfVTP, rfLag);
    
#else

    projectForcesU(nbs, rfDiff, force, rfLag);
    
#endif
    
    // tmp <- inv( J * Jt ) * tmp to find the multipliers
    DPTTS2(nbs, 1, mtJJt, mtJJtU, rfLag, nbs);
}

void Mecafil::storeTensions(const real*)
{
    copy_real(nPoints, rfLLG, rfLag);
}

//------------------------------------------------------------------------------
#pragma mark - Projection DIFF


void Mecafil::makeProjectionDiff(const real* force)
{
    const unsigned nbs = nbSegments();
    assert_true( nbs > 0 );
    
#if ( 0 )
    // verify that we have the correct Lagrange multipliers:
    copy_real(nbs, rfLag, rfLLG);
    computeTensions(force);
    real n = blas::max_diff(nbs, rfLLG, rfLag);
    if ( n > 1e-6 )
    {
        std::clog << "Error= \n" << n << "\n";
        std::clog << "Lagrange: "; VecPrint::print(std::clog, std::min(20u,nbs), rfLLG);
        std::clog << "Multipl.: "; VecPrint::print(std::clog, std::min(20u,nbs), rfLag);
        std::clog << "\n";
    }
#endif
    
    //----- we remove compressive forces ( negative Lagrange-multipliers )
    useProjectionDiff = false;
    for ( unsigned jj = 0; jj < nbs; ++jj )
    {
        if ( rfLag[jj] > 0 )
        {
            useProjectionDiff = true;
            break;
        }
    }
    
    if ( useProjectionDiff )
    {
        const real th = 0.0;
        const real sc = 1.0 / segmentation();
        #pragma vector unaligned
        for ( unsigned jj = 0; jj < nbs; ++jj )
            mtJJtiJforce[jj] = std::max(th, rfLag[jj] * sc);
        
        //std::clog << "projectionDiff: " << blas::nrm2(nbs, mtJJtiJforce) << std::endl;
        //std::clog << "projectionDiff:"; VecPrint::print(std::clog, std::min(20u,nbs), mtJJtiJforce);
    }
}


//------------------------------------------------------------------------------

///straightforward implementation:
inline void add_projectiondiff(const unsigned nbs, const real* mul, const real* X, real* Y)
{
    for ( unsigned jj = 0; jj < nbs; ++jj )
    {
        for ( unsigned d = 0; d < DIM; ++d )
        {
            const real w = mul[jj] * ( X[DIM*jj+DIM+d] - X[DIM*jj+d] );
            Y[DIM*jj    +d] += w;
            Y[DIM*jj+DIM+d] -= w;
        }
    }
}

///expanded implementation:
inline void add_projectiondiffR(const unsigned nbs, const real* mul, const real* X, real* Y)
{
    // this loop cannot be unrolled as there is an OUTPUT dependency in Y
    for ( unsigned jj = 0; jj < nbs; ++jj )
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
inline void add_projectiondiffF(const unsigned nbs, const real* mul, const real* X, real* Y)
{
    real px0 = X[0];
    real px1 = X[1];
    real pw0 = 0;
    real pw1 = 0;
#if ( DIM >= 3 )
    real px2 = X[2];
    real pw2 = 0;
#endif
    
    for ( unsigned jj = 0; jj < nbs; ++jj )
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

inline void add_projectiondiffSSE(const unsigned nbs, const real* mul, const real* X, real* Y)
{
    vec2 px = load2(X);
    vec2 pw = setzero2();
    
    for ( unsigned jj = 0; jj < nbs; ++jj )
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

inline void add_projectiondiffAVX(const unsigned nbs, const real* mul, const real* X, real* Y)
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


void Mecafil::addProjectionDiff(const real* X, real* Y) const
{
    assert_true(useProjectionDiff);
#if ( 0 )
    // debug code to compare with default implementation
    unsigned nbp = nbPoints()*DIM;
    real * vec = new_real(nbp);
    copy_real(nbp, Y, vec);
    add_projectiondiff(nbSegments(), mtJJtiJforce, X, vec);
#endif

#if ( DIM == 2 ) && defined(__SSE3__) && REAL_IS_DOUBLE
    add_projectiondiffSSE(nbSegments(), mtJJtiJforce, X, Y);
    //add_projectiondiff(nbSegments(), mtJJtiJforce, X, Y);
    //add_projectiondiffAVX(nbSegments(), mtJJtiJforce, X, Y);
#else
    add_projectiondiffF(nbSegments(), mtJJtiJforce, X, Y);
#endif
    
    
#if ( 0 )
    // debug code to compare with default implementation
    real n = blas::max_diff(nbp, Y, vec);
    if ( n > 1e-6 )
    {
        std::clog << "proj_diff error " << n << " (" << nbp << ")\n";
        VecPrint::print(std::clog, std::min(20u,nbp), vec);
        VecPrint::print(std::clog, std::min(20u,nbp), Y);
    }
    free_real(vec);
#endif
}


