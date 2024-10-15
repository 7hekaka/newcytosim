// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020

#include "mecafil_code.cc"
#include "exceptions.h"

// required for debugging:
#include "cytoblas.h"
#include "vecprint.h"

/*
 Selection of LAPACK routines, the safest choice
 */

#define DPTTRF lapack::xpttrf
#define DPTTS2 lapack::xptts2


/*
Selection of projectForces() routines optimized for some architectures
*/

#if ( DIM == 3 ) && REAL_IS_DOUBLE && defined(__AVX__)
#  define projectForcesU projectForcesU3D_AVX
#  define projectForcesD projectForcesD3D_AVX
#elif ( DIM == 3 ) && REAL_IS_DOUBLE && USE_SIMD
#  define projectForcesU projectForcesU3D_SSE
#  define projectForcesD projectForcesD3D_SSE
#elif ( DIM == 3 ) && defined(__SSE3__)
#  define projectForcesU projectForcesU_
#  define projectForcesD projectForcesD3D_SSE
#elif ( DIM == 2 ) && REAL_IS_DOUBLE && defined(__AVX__)
#  define projectForcesU projectForcesU2D_AVX
#  define projectForcesD projectForcesD2D_AVX
#elif ( DIM == 2 ) && REAL_IS_DOUBLE && USE_SIMD
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


void Mecafil::allocateProjection(const size_t ms, real* mem)
{
    //std::clog << reference() << "allocateProjection(" << nbp << ")\n";
#if ADD_PROJECTION_DIFF
    //zero_real(3*ms, mem);
    iJJt   = mem;
    iJJtU  = mem + ms;
    iJJtJF = mem + ms * 2;
#else
    //zero_real(2*ms, mem);
    iJJt   = mem;
    iJJtU  = mem + ms;
#endif
}


void Mecafil::destroyProjection()
{
    //std::clog << reference() << "destroyProjection\n";
    iJJt   = nullptr;
    iJJtU  = nullptr;
#if ADD_PROJECTION_DIFF
    iJJtJF = nullptr;
#endif
}


/** This is the standard version assuming isotropic drag coefficients */
void Mecafil::makeProjection()
{
    assert_true( nbPoints() >= 2 );

    //set the diagonal and off-diagonal of J*J'
    const index_t nbu = nbPoints() - 2;

    for ( size_t jj = 0; jj < nbu; ++jj )
    {
        const real* X = iDir + DIM * jj;
#if ( DIM == 2 )
        iJJtU[jj] = -( X[0]*X[2] + X[1]*X[3] );
#else
        iJJtU[jj] = -( X[0]*X[3] + X[1]*X[4] + X[2]*X[5] );
#endif
        
        // the diagonal term should be 2.0, since iDir[] vectors are normalized:
#if ( DIM == 2 )
        iJJt[jj] = 2 * ( X[0]*X[0] + X[1]*X[1] );
#else
        iJJt[jj] = 2 * ( X[0]*X[0] + X[1]*X[1] + X[2]*X[2] );
#endif
        // iJJt[jj]  = 2.0;
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

    if ( 0 )
    {
        VecPrint::print(" D", nbu+1, iJJt, 2);
        VecPrint::print("_U", nbu, iJJtU, 2);
        //VecPrint::print("X", DIM*(nbu+2), pPos, 2);
    }

    if ( info )
    {
        std::clog << "Mecafil::makeProjection failed (" << info << ")\n";
        throw Exception("could not build Fiber's projection matrix");
    }
}


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
     for i in DIM * [ 0, nbs-1 ]
         dst[i] <- src[i] + dif[i] * mul[i] - dif[i-1] * mul[i-1]

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
    
    for ( size_t j = 1; j < nbs; ++j )
    {
        const size_t kk = DIM * j;
        const real M = mul[j], P = mul[j-1];
        dst[kk  ] = src[kk  ] + dir[kk  ] * M - dir[kk-DIM  ] * P;
        dst[kk+1] = src[kk+1] + dir[kk+1] * M - dir[kk-DIM+1] * P;
#if ( DIM > 2 )
        dst[kk+2] = src[kk+2] + dir[kk+2] * M - dir[kk-DIM+2] * P;
#endif
    }
}


//------------------------------------------------------------------------------

/*
 Y <- components of X that are compatible with the length constaints
 Note that this should work correctly even if ( X == Y ), which is always the case
 */
void Mecafil::projectForces(const real* X, real* Y) const
{
    const index_t nbs = nbSegments();
    //VecPrint::print("X", DIM*nbPoints(), X);
    
    // calculate `iLLG` without modifying `X`
    projectForcesU(nbs, iDir, X, iLLG);
    
    // Calculate Lagrange multipliers: iLLG <- inv( J * Jt ) * iLLG
    DPTTS2(nbs, 1, iJJt, iJJtU, iLLG, nbs);

    // set Y, using values in X and multipliers in iLLG
    projectForcesD(nbs, iDir, X, iLLG, Y);

    //VecPrint::print("Y", DIM*nbPoints(), Y);
}


/**
 This sets `iLag` corresponding to the given forces
 */
void Mecafil::computeTensions(const real* force)
{
    const index_t nbs = nbSegments();
    
    projectForcesU(nbs, iDir, force, iLag);
    
    // determine the multipliers: iLag <- inv( J * Jt ) * iLag
    DPTTS2(nbs, 1, iJJt, iJJtU, iLag, nbs);
    
    //fprintf(stderr, "\nmul "); VecPrint::print(stderr, nbs, iLag, 6);
}


/** This extracts the matrix underlying the 'Mecafil::projectForces()' */
void Mecafil::printProjection(FILE * file) const
{
    const index_t nbv = DIM * nbPoints();
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
    fprintf(file, " %s Projection ( %u )\n", reference().c_str(), nbPoints());
    VecPrint::full(file, nbv, nbv, res, nbv);
    free_real(res);
}



//------------------------------------------------------------------------------
#pragma mark - Correction terms to the Projection

#if ADD_PROJECTION_DIFF

// add debug code to compare with reference implementation
#define CHECK_PROJECTION_DIFF 0

/** This assumes that the Lagrange multipliers in 'iLLG' can be used */
void Mecafil::setProjectionDiff(const real threshold)
{
    const index_t nbs = nbSegments();

    // use Lagrange multipliers computed from the last projectForces() in iLLG
    // check for extensile ( positive ) multipliers
    for ( index_t i = 0; i < nbs; ++i )
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
        #pragma omp simd
        for ( size_t jj = 0; jj < nbs; ++jj )
            iJJtJF[jj] = std::max(threshold, alpha * iLLG[jj]);
        iJJtJF[nbs] = 0;
        
        //std::clog << "projectionDiff: " << blas::nrm2(nbs, iJJtJF) << '\n';
        //VecPrint::print("projectionDiff:", std::min(20u,nbs), iJJtJF);
    }
    else // this is not necessary:
        zero_real(nPoints, iJJtJF);
}


/// Reference (scalar) code
/**
 This looks similar to projectForcesD_() except with dir[i] = X[i+DIM]-X[i]
 */
inline void addProjectionDiff_(const size_t nbs, const real* mul, const real* X, real* Y)
{
    for ( size_t i = 0; i < nbs; ++i )
    {
        real const* xx = X + DIM*i;
        real * yy = Y + DIM*i;
        for ( size_t d = 0; d < DIM; ++d )
        {
            const real w = mul[i] * ( xx[DIM+d] - xx[d] );
            yy[    d] += w;
            yy[DIM+d] -= w;
        }
    }
}


/// Add projection-diff matrix
void Mecafil::addProjectionDiff(real* mat) const
{
    index_t nbs = nbSegments();
    index_t bks = DIM * nPoints;
#if CHECK_PROJECTION_DIFF
    real * dup = new_real(bks*bks);
    real * tmp = new_real(bks*bks);
    copy_real(bks*bks, mat, dup);
    zero_real(bks*bks, tmp);
    // use addProjectionDiff() to add all column vectors:
    for ( index_t i = 0; i < bks; ++i )
    {
        tmp[i] = 1;
        addProjectionDiff(tmp, dup+bks*i);
        tmp[i] = 0;
    }
    free_real(tmp);
#endif
    for ( index_t i = 0; i < nbs; ++i )
    {
        real w = iJJtJF[i];
        if ( w != 0 ) for ( unsigned d = 0; d < DIM; ++d )
        {
            real * yy = mat + (1+bks)*(DIM*i+d);
            yy[0  ] -= w;
            yy[DIM] += w;
            real * zz = yy + DIM * bks;
            zz[0  ] += w;
            zz[DIM] -= w;
        }
    }
#if CHECK_PROJECTION_DIFF
    real e = blas::difference(bks*bks, mat, dup);
    if ( e > 1e-6 )
    {
        VecPrint::print("iJJtJF", nPoints-1, iJJtJF);
        VecPrint::full("diffP", bks, bks, dup, bks);
        VecPrint::full("---dP", bks, bks, mat, bks);
    }
    free_real(dup);
#endif
}


void Mecafil::addProjectionDiff(const real* X, real* Y) const
{
#if CHECK_PROJECTION_DIFF
    index_t nbp = DIM * nbPoints();
    real * vec = new_real(nbp);
    copy_real(nbp, Y, vec);
    addProjectionDiff_(nbSegments(), iJJtJF, X, vec);
#endif

#if ( DIM == 2 ) && REAL_IS_DOUBLE && USE_SIMD
    addProjectionDiff2D_SSE(nbSegments(), iJJtJF, X, Y);
    //addProjectionDiff_AVX(nbSegments(), iJJtJF, X, Y);
#elif ( DIM == 3 ) && REAL_IS_DOUBLE && USE_SIMD
    addProjectionDiff3D_SSE(nbSegments(), iJJtJF, X, Y);
#else
    addProjectionDiff_F(nbSegments(), iJJtJF, X, Y);
    //addProjectionDiff_(nbSegments(), iJJtJF, X, Y);
#endif
    
#if CHECK_PROJECTION_DIFF
    real e = blas::difference(nbp, Y, vec);
    if ( e > 1e-6 )
    {
        std::clog << "Mecafil:addProjectionDiff(" << nbp << ") error " << e << "\n";
        VecPrint::edges("ref ", nbp, vec, 3);
        VecPrint::edges("sse ", nbp, Y, 3);
    }
    free_real(vec);
#endif
}

#endif

#if ADD_PROJECTION_DIFF == 7

/** This is the debug pathway */
void Mecafil::makeProjectionDiff(const real* force)
{
    if ( force )
    {
        // compute tensions in 'iLag' using the given forces
        computeTensions(force);
        return;
    }
    
    // Check that `iLLG` contains the same values as `iLag`
    real e = blas::difference(nbSegments(), iLLG, iLag);
    if ( e > 1e-6 )
    {
        std::clog << "Mecafil: |iLag - iLLG| = " << e << "\n";
        VecPrint::print("iLag ", nbSegments(), iLag);
        VecPrint::print("iLLG ", nbSegments(), iLLG);
    }

    setProjectionDiff(0);
}

#elif ADD_PROJECTION_DIFF

/** This is the normal pathway without verifications */
void Mecafil::makeProjectionDiff(const real*)
{
    setProjectionDiff(0);
}

#endif
