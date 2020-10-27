// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "mecafil_code.cc"
#include "exceptions.h"
#include "vecprint.h"
#include "xpttrf.h"


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
Selection of projectForces() routines optimized for some architectures
*/

#if ( DIM == 3 ) && REAL_IS_DOUBLE && defined(__AVX__)
#  define projectForcesU projectForcesU3D_AVX
#  define projectForcesD projectForcesD3D_AVX
#elif ( DIM == 2 ) && REAL_IS_DOUBLE && defined(__AVX__)
#  define projectForcesU projectForcesU2D_AVX
#  define projectForcesD projectForcesD2D_AVX
#elif ( DIM == 2 ) && REAL_IS_DOUBLE && defined(__SSE3__)
#  define projectForcesU projectForcesU2D_SSE
#  define projectForcesD projectForcesD2D_SSE
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
        
        // this term should be 2.0, since iDir[] vectors are normalized:
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
     mul[] <- dot(dif[], src[+DIM] - src[])
 which is:
     mul[i] <- dot(dif[i*DIM], src[i*DIM+DIM] - src[i*DIM])
     for i in [ 0, nbs-1 ]
 
 with 'nbs' = number of segments, and
      dif[] of size nbs*DIM
      vec[] of size (nbs+1)*DIM
      mul[] of size nbs
 
 Note that this should work even if 'mul==src'
 */
void projectForcesU_(size_t nbs, const real* dir, const real* src, real* mul)
{
    const real *const end = mul + nbs;

    while ( mul < end )
    {
        *mul = dir[0] * ( src[DIM  ] - src[0] )
             + dir[1] * ( src[DIM+1] - src[1] )
#if ( DIM > 2 )
             + dir[2] * ( src[DIM+2] - src[2] )
#endif
        ;
        src += DIM;
        dir += DIM;
        ++mul;
    }
}

/**
 Perform second calculation needed by projectForces:
     dst <- src +/- dif * mul
 which is:
     dst[i] <- src[i] + dif[i] * mul[i] - dif[i-1] * mul[i-1]
     for i in DIM * [ 0, nbs-1 ]

 with 'nbs' = number of segments, and
      dif[] of size nbs*DIM
      src[] and dst[] of size (nbs+1)*DIM
      mul[] of size nbs
 
 Note that this should work even if 'dst==src'
 */
void projectForcesD_(size_t nbs, const real* dir, const real* src, const real* mul, real* dst)
{
    for ( size_t s = 0, e = DIM*nbs; s < DIM; ++s, ++e )
    {
        dst[s] = src[s] + dir[s    ] * mul[    0];
        dst[e] = src[e] - dir[e-DIM] * mul[nbs-1];
    }
    
    for ( size_t jj = 1; jj < nbs; ++jj )
    {
        const size_t kk = DIM*jj;
        const real M = mul[jj], P = mul[jj-1];
        dst[kk  ] = src[kk  ] + dir[kk  ] * M - dir[kk-DIM  ] * P;
        dst[kk+1] = src[kk+1] + dir[kk+1] * M - dir[kk-DIM+1] * P;
#if ( DIM > 2 )
        dst[kk+2] = src[kk+2] + dir[kk+2] * M - dir[kk-DIM+2] * P;
#endif
    }
}


//------------------------------------------------------------------------------

/*
 Note that this works fine even if ( X == Y )
 */
void Mecafil::projectForces(const real* X, real* Y) const
{
#if NEW_SKIP_PROJECTION
    if ( skipProjection )
        return copy_real(DIM*nPoints, X, Y);
#endif
    
    const size_t nbs = nbSegments();
    //printf("X  "); VecPrint::print(std::clog, DIM*nbPoints(), X);

    // calculate `iLLG` without modifying `X`
#if NEW_ANISOTROPIC_FIBER_DRAG
    // iLLG is used as a temporary space to store nPoints * DIM scalars:
    scaleTangentially(nPoints, X, iAni, iLLG);
    projectForcesU(nbs, iDir, iLLG, iLLG);
#else
    projectForcesU(nbs, iDir, X, iLLG);
#endif
    
    // Lagrange multipliers <- inv( J * Jt ) * iLLG
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
    // iLLG is used as a temporary space to store nPoints * DIM scalars:
    scaleTangentially(nPoints, force, iAni, iLLG);
    projectForcesU(nbs, iDir, iLLG, iLag);
#else
    projectForcesU(nbs, iDir, force, iLag);
#endif
    
    // tmp <- inv( J * Jt ) * tmp to find the multipliers
    DPTTS2(nbs, 1, iJJt, iJJtU, iLag, nbs);
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
//#include "cytoblas.h"

void Mecafil::makeProjectionDiff(const real* force)
{
    const size_t nbs = nbSegments();
    assert_true( nbs > 0 );
    
#if NEW_SKIP_PROJECTION
    if ( skipProjection )
    {
        useProjectionDiff = false;
        return;
    }
#endif

#if 0
    // Check here that iLLG[] contains the correct Lagrange multipliers
    // compute Lagrange multipliers corresponding to 'force' in iLag:
    computeTensions(force);
    real n = blas::max_diff(nbs, iLLG, iLag);
    if ( n > 1e-6 )
    {
        fprintf(stderr, "\n|iLag - iLLG| = %e", n);
        fprintf(stderr, "\niLag "); VecPrint::print(std::clog, std::min(20LU,nbs), iLag);
        fprintf(stderr, "\niLLG "); VecPrint::print(std::clog, std::min(20LU,nbs), iLLG);
    }
#endif
    
    // use Lagrange multipliers computed from the last projectForces() in iLLG

    // remove compressive forces ( negative Lagrange-multipliers )
    useProjectionDiff = false;
    for ( size_t jj = 0; jj < nbs; ++jj )
    {
        if ( iLLG[jj] > 0 )
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
            iJJtiJforce[jj] = std::max(th, iLLG[jj] * sc);
        
        //std::clog << "projectionDiff: " << blas::nrm2(nbs, iJJtiJforce) << '\n';
        //std::clog << "projectionDiff:"; VecPrint::print(std::clog, std::min(20u,nbs), iJJtiJforce);
    }
}


/// Reference (scalar) code
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
    //add_projectiondiffAVX(nbSegments(), iJJtiJforce, X, Y);
#else
    add_projectiondiffF(nbSegments(), iJJtiJforce, X, Y);
    //add_projectiondiff(nbSegments(), iJJtiJforce, X, Y);
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


