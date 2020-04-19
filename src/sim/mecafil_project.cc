// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "mecafil_code.cc"
#include "exceptions.h"
#include "vecprint.h"
#include "dpttrf.h"


/*
 Selection of LAPACK routines
 The LAPACK implementation is the safest choice
 The Alsatian version is faster as it avoids divisions
 */

#if ( 1 )
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


/*
Selection of projectForces() routines, depending on architecture
*/

#if ( DIM == 2 ) && REAL_IS_DOUBLE
#  if defined(__AVX__)
#    define projectForcesU projectForcesU2D_AVX
#    define projectForcesD projectForcesD2D_AVX
#  elif defined(__SSE3__)
#    warning "Using SSE3 Fiber::projectForces"
#    define projectForcesU projectForcesU2D_SSE
#    define projectForcesD projectForcesD2D_SSE
#  else
#    warning "Using scalar Fiber::projectForces"
#    define projectForcesU projectForcesU_
#    define projectForcesD projectForcesD_
#  endif
#elif ( DIM == 3 ) && REAL_IS_DOUBLE
#  if defined(__AVX__)
#    define projectForcesU projectForcesU3D_AVX
#    define projectForcesD projectForcesD3D_AVX
#  endif
#else
#  warning "Using scalar Fiber::projectForces"
#  define projectForcesU projectForcesU_
#  define projectForcesD projectForcesD_
#endif


//------------------------------------------------------------------------------
#pragma mark -


void Mecafil::buildProjection()
{
    //reset all variables for the projections:
    iJJt        = nullptr;
    iJJtiJforce = nullptr;
}


void Mecafil::allocateProjection(const size_t ms)
{
    //std::clog << reference() << "allocateProjection(" << nbp << ")\n";
    free_real(iJJt);
    real * mem = new_real(3*ms);
    //zero_real(3*ms, mem);
    iJJt        = mem;
    iJJtU       = mem + ms;
    iJJtiJforce = mem + ms * 2;
}


void Mecafil::destroyProjection()
{
    //std::clog << reference() << "destroyProjection\n";
    free_real(iJJt);
    iJJt        = nullptr;
    iJJtU       = nullptr;
    iJJtiJforce = nullptr;
}


#if NEW_ANISOTROPIC_FIBER_DRAG

/// version for drag coefficients doubled in directions orthogonal to fiber
void Mecafil::makeProjection()
{
    assert_true( nbPoints() >= 2 );

    //set the diagonal and off-diagonal of J*J'
    const size_t nbu = nbPoints() - 2;
    real b = 1;

    for ( size_t jj = 0; jj < nbu; ++jj )
    {
        const real* X = iDir + DIM * jj;
#if ( DIM == 2 )
        real xn = X[0]*X[2] + X[1]*X[3];
#else
        real xn = X[0]*X[3] + X[1]*X[4] + X[2]*X[5];
#endif
        real a = 0.25 * ( 1 + xn ) * ( 1 + xn );
        iJJt[jj]  = 2.0 + a + b;
        iJJtU[jj] = -xn - a;
        b = a;
    }
    
    iJJt[nbu] = 3.0 + b;
    
    int info = 0;
    DPTTRF(nbu+1, iJJt, iJJtU, &info);

    if ( 0 )
    {
        std::clog << "\nD "; VecPrint::print(std::clog, nbu+1, iJJt, 3);
        std::clog << "\nU "; VecPrint::print(std::clog, nbu, iJJtU, 3);
        //std::clog << "\nX="; VecPrint::print(std::clog, DIM*(nbu+2), pPos);
    }

    if ( info )
    {
        std::clog << "Mecafil::makeProjection failed (" << info << ")\n";
        throw Exception("could not build Fiber's projection matrix");
    }
}

#else

/// standard version with isotropic drag coefficient
void Mecafil::makeProjection()
{
    assert_true( nbPoints() >= 2 );

    //set the diagonal and off-diagonal of J*J'
    const size_t nbu = nbPoints() - 2;

    for ( size_t jj = 0; jj < nbu; ++jj )
    {
        const real* X = iDir + DIM * jj;
#if ( DIM == 2 )
        real xn = X[0]*X[2] + X[1]*X[3];
#else
        real xn = X[0]*X[3] + X[1]*X[4] + X[2]*X[5];
#endif
        
        // this term should be 2, since iDir[] vectors are normalized:
#if ( DIM == 2 )
        iJJt[jj] = 2.0 * ( X[0]*X[0] + X[1]*X[1] );
#else
        iJJt[jj] = 2.0 * ( X[0]*X[0] + X[1]*X[1] + X[2]*X[2] );
#endif
        // iJJt[jj]  = 2.0;
        
        iJJtU[jj] = -xn;
    }
    
    const real* X = iDir + DIM*nbu;
    // this term should be 2, since iDir[] vectors are normalized
#if ( DIM == 2 )
    iJJt[nbu] = 2.0 * ( X[0]*X[0] + X[1]*X[1] );
#else
    iJJt[nbu] = 2.0 * ( X[0]*X[0] + X[1]*X[1] + X[2]*X[2] );
#endif
    //iJJt[nbu] = 2.0;

    int info = 0;
    DPTTRF(nbu+1, iJJt, iJJtU, &info);

    if ( 0 )
    {
        std::clog << "\nD "; VecPrint::print(std::clog, nbu+1, iJJt, 3);
        std::clog << "\nU "; VecPrint::print(std::clog, nbu, iJJtU, 3);
        //std::clog << "\nX="; VecPrint::print(std::clog, DIM*(nbu+2), pPos);
    }

    if ( info )
    {
        std::clog << "Mecafil::makeProjection failed (" << info << ")\n";
        throw Exception("could not build Fiber's projection matrix");
    }
}

#endif
//------------------------------------------------------------------------------
#pragma mark - Reference (scalar) code

/**
 Perform first calculation needed by projectForces:
     mul[] <- dot(dif[], X[+DIM]-X[])
 which is:
     mul[i] <- dot(dif[i*DIM], X[i*DIM+DIM]-X[i*DIM])
     for i in [ 0, nbs-1 ]
 
 with 'nbs' = number of segments, and
      dif[] of size nbs*DIM
      vec[] of size (nbs+1)*DIM
      mul[] of size nbs
 */
void projectForcesU_(size_t nbs, const real* dif, const real* src, real* mul)
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
 Perform second calculation needed by projectForces:
     Y <- X +/- dif * mul
 which is:
     Y[i] <- X[i] + dif[i] * mul[i] - dif[i-1] * mul[i-1]
     for i in DIM * [ 0, nbs-1 ]

 with 'nbs' = number of segments, and
      dif[] of size nbs*DIM
      X[] and Y[] of size (nbs+1)*DIM
      mul[] of size nbs
 */
void projectForcesD_(size_t nbs, const real* dif, const real* X, const real* mul, real* Y)
{
    for ( size_t s = 0, e = DIM*nbs; s < DIM; ++s, ++e )
    {
        Y[s] = X[s] + dif[s    ] * mul[    0];
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


//------------------------------------------------------------------------------
#pragma mark -

/*
 Note that this works fine even if ( X == Y )
 */
void Mecafil::projectForces(const real* X, real* Y) const
{
    const size_t nbs = nbSegments();
    //printf("X  "); VecPrint::print(std::clog, DIM*nbPoints(), X);

    // calculate `iLLG` without modifying `X`
#if NEW_ANISOTROPIC_FIBER_DRAG
    scaleTangentially(nPoints, X, iAni, iVTP);
    projectForcesU(nbs, iDir, iVTP, iLLG);
#else
    projectForcesU(nbs, iDir, X, iLLG);
#endif
    
    // iLLG <- inv( J * Jt ) * iLLG to find the Lagrange multipliers
    DPTTS2(nbs, 1, iJJt, iJJtU, iLLG, nbs);

    // set Y, using values in X and iLLG
    projectForcesD(nbs, iDir, X, iLLG, Y);

#if NEW_ANISOTROPIC_FIBER_DRAG
    scaleTangentially(nPoints, Y, iAni, Y);
#endif
    //printf("Y  "); VecPrint::print(std::clog, DIM*nbPoints(), Y);
}


void Mecafil::computeTensions(const real* force)
{
    const size_t nbs = nbSegments();
    
#if NEW_ANISOTROPIC_FIBER_DRAG
    
    scaleTangentially(nPoints, force, iAni, iVTP);
    projectForcesU(nbs, iDir, iVTP, iLag);
    
#else

    projectForcesU(nbs, iDir, force, iLag);
    
#endif
    
    // tmp <- inv( J * Jt ) * tmp to find the multipliers
    DPTTS2(nbs, 1, iJJt, iJJtU, iLag, nbs);
}


void Mecafil::storeTensions(const real*)
{
    copy_real(nPoints, iLLG, iLag);
}


void Mecafil::printProjection(std::ostream& os) const
{
    const size_t nbv = DIM * nbPoints();
    real * res = new_real(nbv*nbv);
    real * src = new_real(nbv);
    real * dst = new_real(nbv);
    zero_real(nbv, src);
    zero_real(nbv, dst);
    for ( size_t i = 0; i < nbv; ++i )
    {
        src[i] = 1.0;
        projectForces(src, dst);
        copy_real(nbv, dst, res+nbv*i);
        src[i] = 0.0;
    }
    free_real(dst);
    free_real(src);
    os << "Mecafil:Projection  " << reference() << '\n';
    VecPrint::print(os, nbv, nbv, res, nbv);
    free_real(res);
}



//------------------------------------------------------------------------------
#pragma mark - Projection DIFF


void Mecafil::makeProjectionDiff(const real* force)
{
    const size_t nbs = nbSegments();
    assert_true( nbs > 0 );
    
#if ( 0 )
    // verify that we have the correct Lagrange multipliers:
    copy_real(nbs, iLag, iLLG);
    computeTensions(force);
    real n = blas::max_diff(nbs, iLLG, iLag);
    if ( n > 1e-6 )
    {
        std::clog << "Error= \n" << n << "\n";
        std::clog << "Lagrange: "; VecPrint::print(std::clog, std::min(20u,nbs), iLLG);
        std::clog << "Multipl.: "; VecPrint::print(std::clog, std::min(20u,nbs), iLag);
        std::clog << "\n";
    }
#endif
    
    //----- we remove compressive forces ( negative Lagrange-multipliers )
    useProjectionDiff = false;
    for ( size_t jj = 0; jj < nbs; ++jj )
    {
        if ( iLag[jj] > 0 )
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
        for ( size_t jj = 0; jj < nbs; ++jj )
            iJJtiJforce[jj] = std::max(th, iLag[jj] * sc);
        
        //std::clog << "projectionDiff: " << blas::nrm2(nbs, iJJtiJforce) << std::endl;
        //std::clog << "projectionDiff:"; VecPrint::print(std::clog, std::min(20u,nbs), iJJtiJforce);
    }
}


//------------------------------------------------------------------------------

///straightforward implementation:
inline void add_projectiondiff(const size_t nbs, const real* mul, const real* X, real* Y)
{
    for ( size_t jj = 0; jj < nbs; ++jj )
    {
        for ( size_t d = 0; d < DIM; ++d )
        {
            const real w = mul[jj] * ( X[DIM*jj+DIM+d] - X[DIM*jj+d] );
            Y[DIM*jj    +d] += w;
            Y[DIM*jj+DIM+d] -= w;
        }
    }
}

///expanded implementation:
inline void add_projectiondiffR(const size_t nbs, const real* mul, const real* X, real* Y)
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
inline void add_projectiondiffF(const size_t nbs, const real* mul, const real* X, real* Y)
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

inline void add_projectiondiffSSE(const size_t nbs, const real* mul, const real* X, real* Y)
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

inline void add_projectiondiffAVX(const size_t nbs, const real* mul, const real* X, real* Y)
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
    size_t nbp = nbPoints()*DIM;
    real * vec = new_real(nbp);
    copy_real(nbp, Y, vec);
    add_projectiondiff(nbSegments(), iJJtiJforce, X, vec);
#endif

#if ( DIM == 2 ) && defined(__SSE3__) && REAL_IS_DOUBLE
    add_projectiondiffSSE(nbSegments(), iJJtiJforce, X, Y);
    //add_projectiondiff(nbSegments(), iJJtiJforce, X, Y);
    //add_projectiondiffAVX(nbSegments(), iJJtiJforce, X, Y);
#else
    add_projectiondiffF(nbSegments(), iJJtiJforce, X, Y);
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


