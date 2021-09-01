// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020

#include "mecafil_code.cc"
#include "exceptions.h"
#include "xpttrf.h"

// required for debugging:
#include "cytoblas.h"
#include "vecprint.h"

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
#elif ( DIM == 3 ) && REAL_IS_DOUBLE && defined(__SSE3__)
#  define projectForcesU projectForcesU3D_SSE
#  define projectForcesD projectForcesD3D_SSE
#elif ( DIM == 3 ) && defined(__SSE3__)
#  define projectForcesU projectForcesU_
#  define projectForcesD projectForcesD3D_SSE
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


void Mecafil::initProjection()
{
    //reset all variables for the projections:
    iJJt   = nullptr;
#if ADD_PROJECTION_DIFF
    iJJtJF = nullptr;
#endif
}


void Mecafil::allocateProjection(const size_t ms)
{
    //std::clog << reference() << "allocateProjection(" << nbp << ")\n";
    free_real(iJJt);
#if ADD_PROJECTION_DIFF
    real * mem = new_real(3*ms);
    //zero_real(3*ms, mem);
    iJJt   = mem;
    iJJtU  = mem + ms;
    iJJtJF = mem + ms * 2;
#else
    real * mem = new_real(2*ms);
    //zero_real(2*ms, mem);
    iJJt   = mem;
    iJJtU  = mem + ms;
#endif
}


void Mecafil::destroyProjection()
{
    //std::clog << reference() << "destroyProjection\n";
    free_real(iJJt);
    iJJt   = nullptr;
    iJJtU  = nullptr;
#if ADD_PROJECTION_DIFF
    iJJtJF = nullptr;
#endif
}


#if NEW_ANISOTROPIC_FIBER_DRAG

/** This version assumes that drag coefficients are double in the directions
 orthogonal to the fiber axis, ie. it works with anisotropic fiber drag  */
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
        VecPrint::print("D", nbu+1, iJJt, 2);
        VecPrint::print("U", nbu, iJJtU, 2);
        //VecPrint::print("X", DIM*(nbu+2), pPos, 2);
    }

    if ( info )
    {
        std::clog << "Mecafil::makeProjection failed (" << info << ")\n";
        throw Exception("could not build Fiber's projection matrix");
    }
}

#else

/** This is the standard version assuming isotropic drag coefficients */
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
        iJJt[jj] = 2 * ( X[0]*X[0] + X[1]*X[1] );
#else
        iJJt[jj] = 2 * ( X[0]*X[0] + X[1]*X[1] + X[2]*X[2] );
#endif
        // iJJt[jj]  = 2.0;
        
        iJJtU[jj] = -xn;
    }
    
    const real* X = iDir + DIM*nbu;
    // this term should be 2, since iDir[] vectors are normalized
#if ( DIM == 2 )
    iJJt[nbu] = 2 * ( X[0]*X[0] + X[1]*X[1] );
#else
    iJJt[nbu] = 2 * ( X[0]*X[0] + X[1]*X[1] + X[2]*X[2] );
#endif
    //iJJt[nbu] = 2.0;

    int info = 0;
    DPTTRF(nbu+1, iJJt, iJJtU, &info);

    if ( 1 )
    {
        VecPrint::print("D", nbu+1, iJJt, 2);
        VecPrint::print("U", nbu, iJJtU, 2);
        //VecPrint::print("X", DIM*(nbu+2), pPos, 2);
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
void projectForcesD_(const size_t nbs, const real* dir, const real* src, const real* mul, real* dst)
{
    for ( size_t s = 0; s < DIM; ++s )
        dst[s] = src[s] + dir[s] * mul[0];
    
    for ( size_t e = DIM*nbs; e < DIM*(nbs+1); ++e )
        dst[e] = src[e] - dir[e-DIM] * mul[nbs-1];
    
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
#if UNCONSTRAINED_LENGTH
    if ( unconstrainLength )
        return copy_real(DIM*nPoints, X, Y);
#endif
    
    const size_t nbs = nbSegments();
    //VecPrint::print("X", DIM*nbPoints(), X);

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
    //VecPrint::print("Y", DIM*nbPoints(), Y);
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


/** This extracts the matrix underlying the 'Mecafil::projectForces()' */
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
    os << "Mecafil:Projection  " << reference() << " (" << nbPoints() << ")\n";
    VecPrint::print(os, nbv, nbv, res, nbv);
    free_real(res);
}



//------------------------------------------------------------------------------
#pragma mark - Correction terms to the Projection

#if ADD_PROJECTION_DIFF

// add debug code to compare with reference implementation
#define CHECK_PROJECTION_DIFF 0


void Mecafil::makeProjectionDiff(const real* force)
{
    useProjectionDiff = false;
#if UNCONSTRAINED_LENGTH
    if ( unconstrainLength )
        return;
#endif
    const size_t nbs = nbSegments();
    assert_true( nbs > 0 );
    
#if CHECK_PROJECTION_DIFF
    // Check here that iLLG[] contains the correct Lagrange multipliers
    // compute Lagrange multipliers corresponding to 'force' in iLag:
    computeTensions(force);
    real e = blas::difference(nbs, iLLG, iLag);
    if ( e > 1e-6 )
    {
        std::clog << "\n|iLag - iLLG| = " << e << "\n";
        VecPrint::print("\niLag ", std::min(16LU,nbs), iLag);
        VecPrint::print("\niLLG ", std::min(16LU,nbs), iLLG);
    }
#endif
    
    const real threshold = 0.0;

    // use Lagrange multipliers computed from the last projectForces() in iLLG
    // check for extensile ( positive ) multipliers
    for ( size_t i = 0; i < nbs; ++i )
    {
        if ( iLLG[i] > threshold )
        {
            useProjectionDiff = true;
            break;
        }
    }
    
    // remove compressive ( negative ) multipliers
    if ( useProjectionDiff )
    {
        const real alpha = 1.0 / segmentation();
        #pragma vector unaligned
        for ( size_t jj = 0; jj < nbs; ++jj )
            iJJtJF[jj] = std::max(threshold, alpha * iLLG[jj]);
        
        //std::clog << "projectionDiff: " << blas::nrm2(nbs, iJJtJF) << '\n';
        //VecPrint::print("projectionDiff:", std::min(20u,nbs), iJJtJF);
    }
}


/// Reference (scalar) code
inline void addProjectionDiff_(const size_t nbs, const real* mul, const real* X, real* Y)
{
    for ( size_t i = 0; i < nbs; ++i )
    {
        for ( size_t d = 0; d < DIM; ++d )
        {
            const real w = mul[i] * ( X[DIM*i+DIM+d] - X[DIM*i+d] );
            Y[DIM*i+(    d)] += w;
            Y[DIM*i+(DIM+d)] -= w;
        }
    }
}


void Mecafil::addProjectionDiff(const real* X, real* Y) const
{
#if CHECK_PROJECTION_DIFF
    size_t nbp = nbPoints()*DIM;
    real * vec = new_real(nbp);
    copy_real(nbp, Y, vec);
    addProjectionDiff_(nbSegments(), iJJtJF, X, vec);
#endif

#if ( DIM == 2 ) && defined(__SSE3__) && REAL_IS_DOUBLE
    addProjectionDiff_SSE(nbSegments(), iJJtJF, X, Y);
    //addProjectionDiff_AVX(nbSegments(), iJJtJF, X, Y);
#else
    addProjectionDiff_F(nbSegments(), iJJtJF, X, Y);
    //addProjectionDiff_(nbSegments(), iJJtJF, X, Y);
#endif
    
#if CHECK_PROJECTION_DIFF
    real e = blas::difference(nbp, Y, vec);
    if ( e > 1e-6 )
    {
        std::clog << "\naddProjectionDiff(" << nbp << ") error " << e;
        VecPrint::print("\nref ", std::min(16UL,nbp), vec);
        VecPrint::print("\nopt ", std::min(16UL,nbp), Y);
    }
    free_real(vec);
#endif
}

#endif
