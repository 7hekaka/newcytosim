// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "assert_macro.h"
#include "chain.h"
#include "iowrapper.h"
#include "messages.h"
#include "mecapoint.h"
#include "interpolation.h"
#include "fiber_site.h"
#include "exceptions.h"
#include "glossary.h"
#include "blas.h"
#include "lapack.h"
#include "cytoblas.h"
#include "modulo.h"
#include "vector3.h"
#include "simul.h"
#include "vecprint.h"


// defined in "mecafil_code.cc"
extern void projectForcesD_(size_t, const real*, const real*, const real*, real*);


extern Modulo const* modulo;

/**
 This returns N+1, where N is the integer that minimizes
     abs_real( length / N - segmentation ),
 */
size_t Chain::bestNumberOfPoints(const real ratio)
{
    real n = 1 + std::floor(ratio);
    return (size_t)( n + ( ratio*(n-0.5) >= n*(n-1) ));
}


real Chain::contourLength(const real pts[], size_t n_pts)
{
    real len = 0;
    Vector a(pts), b;
    for ( size_t n = 1; n < n_pts; ++n )
    {
        b.load(pts+DIM*n);
        len += (b-a).norm();
        a = b;
    }
    return len;
}


/** number of time steps between each attempt to remove/add a point
 also number of states over which curvature information is averaged */
#define RECUT_PERIOD 8

Chain::Chain()
{
    fnCut          = 0;
    fnSegmentation = 0;
#if NEW_SKIP_PROJECTION
    skipProjection = false;
#endif
#if CURVATURE_DEPENDENT_SEGMENTATION
    autoCutVal  = 0;
    autoCutCnt  = 0;
    // set different 'seeds' to desynchronize the adjustSegmentation()
    autoCutIdx  = RNG.pint(RECUT_PERIOD);
#endif
    fnAbscissaM = 0;
    fnAbscissaP = 0;
    fnBirthTime = 0;
    needUpdate  = false;
    //iDirValid   = false;
#if FIBER_HAS_NORMAL
    fnNormal.set(0, 0, 1);
#endif
}


//------------------------------------------------------------------------------
#pragma mark -

/**
This does not change length or segmentation
*/
void Chain::setStraight(Vector const& pos, Vector const& dir)
{
    assert_true( pos.valid() );
    assert_true( dir.norm() > 0.1 );
    // 'dir' is normalized for safety:
    Vector dpts = dir * ( fnCut / dir.norm() );
    //
    for ( size_t p = 0 ; p < nPoints; ++p )
        setPoint(p, pos + p * dpts);
}


void Chain::setStraight(Vector const& pos, Vector const& dir, real len)
{
    assert_true( fnSegmentation > REAL_EPSILON );

    if ( len <= 0 )
        throw InvalidParameter("fiber:length must be > 0");

    size_t np = bestNumberOfPoints(len/fnSegmentation);
    
    setNbPoints(np);
    setSegmentation(len/(np-1));
    fnAbscissaP = fnAbscissaM + len;
    setStraight(pos, dir);
    updateFiber();
}


void Chain::placeEnd(const FiberEnd ref)
{
    switch( ref )
    {
        case NO_END:
            flipChainPolarity();
            break;
        
        case MINUS_END:
            translate(posMiddle()-posEndM());
            break;
            
        case PLUS_END:
            translate(posMiddle()-posEndM());
            flipChainPolarity();
            break;
            
        case CENTER:
            break;
            
        default:
            ABORT_NOW("invalid argument to Chain::placeEnd()");
    }
}


/**
 This will set the Fiber with `n_pts` points unless `n_pts == 0`, in which case
 the number of points will be set automatically from `fnSegmentation`.
 `pts[]` should provide `DIM * n_pts` coordinates.

 The given set of points do not need to be equally distributed.
 The MINUS_END and PLUS_END will be set to the first and last points in `pts[]`,
 and intermediate points will be interpolated at regular intervals on `pts[]`.
 
 The length of the resulting fiber will match the sum of given segment lengths,
 and the segments will be approximately equal to each other.
 Thus reshape() is called eventually to equalize the segments.
 */
void Chain::setShape(const real pts[], size_t n_pts, size_t np)
{
    assert_true(n_pts > 1);
    Vector a(pts), b;
    
    //calculate the total length
    real len = contourLength(pts, n_pts);
    
    if ( np == 0 )
    {
        assert_true( fnSegmentation > REAL_EPSILON );
        np = bestNumberOfPoints(len/fnSegmentation);
    }
    setNbPoints(np);
    setSegmentation(len/(np-1));
    fnAbscissaP = fnAbscissaM + len;
    
    a.load(pts);
    b.load(pts+DIM);
    setPoint(0, a);
    
    len = (b-a).norm();
    real h = 0;
    size_t p = 1;
    --np;
    
    for ( size_t n = 1; n < np; ++n )
    {
        h += fnCut;

        while ( h > len )
        {
            h -= len;
            a = b;
            ++p;
            assert_true(p<n_pts);
            b.load(pts+DIM*p);
            len = (b-a).norm();
        }
        
        setPoint(n, a+(h/len)*(b-a));
    }
    b.load(pts+DIM*n_pts-DIM);
    setPoint(np, b);
    //updateFiber();
    reshape();
}

/**
 The filament is set as a random walk with given persistence length

 This return a filament in a random direction, with the center of gravity at zero
 and the average orientation aligned with (1, 0, 0)
 */
void Chain::setEquilibrated(real len, real persistence_length)
{
    size_t np = bestNumberOfPoints(len/fnSegmentation);
    assert_true( np > 1 );
    
    setNbPoints(np);
    setSegmentation(len/(np-1));
    fnAbscissaP = fnAbscissaM + len;
    
    real sigma = sqrt(2*fnCut/persistence_length);
    
    Vector pos(0,0,0);
    Vector dir(1,0,0);
    setPoint(0, pos);
    
    for ( size_t p = 1 ; p < np; ++p )
    {
        pos += fnCut * dir;
        setPoint(p, pos);
        //rotate dir in a random direction:
        real a = sigma * RNG.gauss();
        dir = cos(a) * dir + dir.randOrthoU(sin(a));
    }
    
    // cancel out mean orientation and position:
    translate(-0.5*pos);
    if ( pos.normSqr() > 0.01 * fnCut )
    {
        Rotation rot = Rotation::rotationToVector(pos).transposed();
        rotate(rot);
    }
    updateFiber();
}


/**
 This adjusts the current `normal` or makes a new one if necessary
 (that is only used for display following the 'faked' actin or microtubule style)
 */
Vector3 Chain::adjustedNormal(Vector3 const& d) const
{
#if FIBER_HAS_NORMAL
    if ( fnNormal.normSqr() < 0.8 || dot(fnNormal, d) > 0.5 )
        fnNormal = d.orthogonal(1.0);
    else
        fnNormal = d.orthogonal(fnNormal, 1.0);
    return fnNormal;
#else
    return Vector3(1,0,0);
#endif
}


double Chain::age() const
{
    return simul().time() - fnBirthTime;
}


//===================================================================
#pragma mark -

/*
 This deals with Fiber having one segment only,
 for which the procedure is trivial
 */
void Chain::reshape_two(const real* src, real* dst, real cut)
{
    real X = src[  DIM] - src[0];
#if ( DIM == 1 )
    real s = 0.5 - 0.5 * (cut/abs_real(X));
#elif ( DIM == 2 )
    real Y = src[1+DIM] - src[1];
    real n = sqrt( X * X + Y * Y );
    real s = 0.5 - 0.5 * (cut/n);
#else
    real Y = src[1+DIM] - src[1];
    real Z = src[2+DIM] - src[2];
    real n = sqrt( X * X + Y * Y + Z * Z );
    real s = 0.5 - 0.5 * (cut/n);
#endif
    
    dst[0    ] = src[0    ] + s * X;
    dst[  DIM] = src[  DIM] - s * X;
#if ( DIM > 1 )
    dst[1    ] = src[1    ] + s * Y;
    dst[1+DIM] = src[1+DIM] - s * Y;
#endif
#if ( DIM > 2 )
    dst[2    ] = src[2    ] + s * Z;
    dst[2+DIM] = src[2+DIM] - s * Z;
#endif
}


/**
 Shorten segments to restore their length to 'cut'.
 We use a multidimensional Newton's method, to find iteratively the scalar
 coefficients that define the amount of displacement of each point.
 
     X[i] = vector of position
 
 We note 'dif' the differences between consecutive points:  dif[i] = X[i+1] - X[i]
 Given one scalar per segment: A[i], the point is displaced as:
 
     Y[i] = X[i] + A[i] * dif[i] - A[i-1] * dif[i-1]
 
 except for the first and last points, for which there is only one term:
 
     Y[0] = X[0] + A[  0] * dif[  0]
     Y[L] = X[L] - A[L-1] * dif[L-1]
 
 We want 'A[]' to restore the length of segments:
 
     ( Y[i+1] - Y[i] )^2 = cut^2
 
 i.e. 'A[]' should fulfill a set of equalities F[i] = 0, with:
 
     F[i] = ( Y[i+1] - Y[i] )^2 - cut^2
 
 Note that:
 
     Y[i+1] - Y[i] = A[i+1] * dif[i+1] + (1-2*A[i]) * dif[i] + A[i-1] * dif[i-1]
 
 Method: use all zeros as first guess for 'sca', and apply a multidimensional
 Newton's method to iteratively refine the guess.
 
 In practice, we calculate `A_next` from `A` using the relationship:
 
     J(A) * ( A_next - A ) = -F(A)
 
 Where J is the Jacobian matrix: J[i,j] = dF[i] / dA[j]
 
 For this problem, J is square and tri-diagonal but not symmetric,
 and must be recalculated at each iteration. A factor 2 can be factorized:

     A_next = A - 1/2 inv(K).F(A)
 
 Where J = 2 * K

 externally provided memory `mem[]` should be allocated to hold `5*chk` reals
 FJN, Strasbourg, 22.02.2015 & Cambridge, 10.05.2019 -- 13.05.2019
 */

int Chain::reshape_calculate(const size_t ns, real cutcut,
                             real const* mag, real const* pri, real const* sec,
                             real* mem, size_t chk)
{
    real * sca = mem;
    real * val = mem+chk;
    real * dia = mem+chk*2;
    real * low = mem+chk*3;
    real * upe = mem+chk*4;
    
    real err0 = 0;
    /*
     Perform here the first iteration of Newton's method
     the formula is the same as below, with all `sca` equal to zero,
     and thus 'vec == dif'
     The system is symmetric, and we can use a faster factorization
     */
    for ( size_t i = 0; i < ns; ++i )
    {
        sca[i] = mag[i] - cutcut;
        err0 += abs_real(sca[i]);  // calculating the 1-norm
        dia[i] = mag[i] * 4;
        low[i] = pri[i] * (-2);  //accessing pri[ns-1], but not used
    }
    
    int info = 0;
    lapack::xpttrf(ns, dia, low, &info);
    if ( info ) {
        std::cerr << " reshape_local failed " << info << '\n';
        return 1;
    }
    lapack::xptts2(ns, 1, dia, low, sca, ns);
    
#if ( 0 )
    printf("\n ====err %20.16f", err0);
    printf("\n     sca "); VecPrint::print(std::cout, ns, sca, 3);
#endif
    size_t cnt = 0;
    while ( ++cnt < 8 )
    {
        assert_true( ns > 1 );
        // set the matrix elements and RHS of system,
        {
            real b = 1-2*sca[0], c = sca[1];
            //vec = b*dif[i] + c*dif[i+1];
            real D = b * mag[0] + c * pri[0];
            real U = b * pri[0] + c * mag[1];
            val[0] = b * D + c * U - cutcut;
            dia[0] = D * ( -2 );
            upe[0] = U;
        }
        real err = abs_real(val[0]);
        for( size_t i = 1; i+1 < ns; ++i )
        {
            real a = sca[i-1];
            real b = 1 - 2 * sca[i];
            real c = sca[i+1];
            /*
             vec = a*dif[i-1] + b*dif[i] + c*dif[i+1];
             val = vec.normSqr() - curSqr;
             low = derivate(val, sca[i-1])
             dia = derivate(val, sca[i]);
             upe = derivate(val, sca[i+1);
             */
            real L = a * mag[i-1] + b * pri[i-1] + c * sec[i-1];
            real D = a * pri[i-1] + b * mag[i  ] + c * pri[i  ];
            real U = a * sec[i-1] + b * pri[i  ] + c * mag[i+1];

            val[i] = ( a * L + b * D ) + ( c * U - cutcut );
            low[i] = L;
            dia[i] = D * ( -2 );
            upe[i] = U;
            err += abs_real(val[i]);
        }
        {
            real a = sca[ns-2];
            real b = 1 - 2 * sca[ns-1];
            //vec = a*dif[ns-2] + b*dif[ns-1];
            real L = a * mag[ns-2] + b * pri[ns-2];
            real D = a * pri[ns-2] + b * mag[ns-1];
            val[ns-1] = a * L + b * D - cutcut;
            low[ns-1] = L;
            dia[ns-1] = D * ( -2 );
        }
        err += abs_real(val[ns-1]);
#if ( 0 )
        printf("\n %3lu err %20.16f norm(val) %8.5e", cnt, err, blas::nrm2(ns, val));
        //printf("\n     val "); VecPrint::print(std::cout, ns, val, 6);
        //printf("\n     sca "); VecPrint::print(std::cout, ns, sca, 6);
#endif
        if ( err < REAL_EPSILON )
            return 0;
        if ( err > err0 )
            return 3;
        err0 = err;
#if ( 0 )
        printf("\n diff(L,U) = %f", blas::max_diff(ns-1, upe, low+1));
        printf("\n L"); VecPrint::print(std::cout, std::min(16UL, ns-1), low+1, 3);
        printf("\n D"); VecPrint::print(std::cout, std::min(16UL, ns  ), dia, 3);
        printf("\n U"); VecPrint::print(std::cout, std::min(16UL, ns-1), upe, 3);
#endif
#if ( 0 )
        real asy = 0, sup = 0;
        for ( size_t i = 0; i < ns-1; ++i )
        {
            sup = std::max(sup, abs_real(low[i+1]+upe[i]));
            asy += abs_real(low[i+1]-upe[i]);
        }
        printf("\n %3i diff(low-upe) %12.6f", cnt, 2*asy/sup);
#endif
        lapack::xgtsv(ns, 1, low+1, dia, upe, val, ns, &info);
        if ( info )
        {
            std::cerr << " LAPACK dgtsv failed " << info << '\n';
            return 2;
        }
        
        // update result following Newton's iteration
        for ( size_t u = 0; u < ns; ++u )
            sca[u] -= 0.5 * val[u];
        
        //printf("\n   ->sca "); VecPrint::print(std::cout, ns, sca, 3);
    }
    //printf("\n   >>err %20.16f", err);
    //printf("\n   >>sca "); VecPrint::print(std::cout, ns, sca, 3);
    return 4;
}


/// old version (2018)
int Chain::reshape_calculate_old(const size_t ns, real cutcut, const real* dif,
                                 real* mem, size_t chk)
{
    real * sca = mem;
    real * val = mem+chk;
    real * dia = mem+chk*2;
    real * low = mem+chk*3;
    real * upe = mem+chk*4;
    assert_true( ns > 1 );

    /*
     Perform here the first iteration of Newton's method
     the formula is the same as below, with all `sca` equal to zero,
     and thus 'vec == dif'
     The system is symmetric, and we can use a faster factorization
     */
    real err0 = 0;
    for ( size_t i = 0; i < ns; ++i )
    {
        Vector difA(dif+DIM*i), difB(dif+DIM*(i+1));
        real n = difA.normSqr();
        sca[i] = n - cutcut;
        err0 += abs_real(sca[i]);
        dia[i] = n * 4;
        low[i] = dot(difA, difB) * (-2);  //using undefined value
    }
    
    int info = 0;
    lapack::xpttrf(ns, dia, low, &info);
    if ( info ) {
        std::cerr << " reshape_local failed " << info << '\n';
        return 1;
    }
    lapack::xptts2(ns, 1, dia, low, sca, ns);

    //printf("\n ----err %20.16f", err0);
    //printf("\n     sca "); VecPrint::print(std::cout, ns, sca, 3);

    size_t cnt = 0;
    while ( ++cnt < 16 )
    {
        {
            Vector dif0(dif), dif1(dif+DIM);
            // set the matrix elements and RHS of system,
            Vector vec = (1-2*sca[0]) * dif0 + sca[1] * dif1;
            val[0] = vec.normSqr() - cutcut;
            dia[0] = dot(vec, dif0) * ( -2 );
            upe[0] = dot(vec, dif1);
        }
        real err = abs_real(val[0]);
        for ( size_t i = 1; i+1 < ns; ++i )
        {
            Vector difA(dif+DIM*(i-1)), difB(dif+DIM*i), difC(dif+DIM*(i+1));
            Vector vec = sca[i-1]*difA + (1-2*sca[i])*difB + sca[i+1]*difC;
            val[i] = vec.normSqr() - cutcut;
            low[i] = dot(vec, difA);
            dia[i] = dot(vec, difB) * ( -2 );
            upe[i] = dot(vec, difC);
            err += abs_real(val[i]);
        }
        {
            Vector difA(dif+DIM*(ns-2)), difB(dif+DIM*(ns-1));
            Vector vec = sca[ns-2]*difA + (1-2*sca[ns-1])*difB;
            val[ns-1] = vec.normSqr() - cutcut;
            low[ns-1] = dot(vec, difA);
            dia[ns-1] = dot(vec, difB) * ( -2 );
        }
        err += abs_real(val[ns-1]);
#if ( 0 )
        printf("\n %3i err %20.16f norm(val) %8.5f", cnt, err, blas::nrm2(ns, val));
        printf("\n     val "); VecPrint::print(std::cout, ns, val, 3);
        printf("\n     sca "); VecPrint::print(std::cout, ns, sca, 3);
#endif
        if ( err < 1e-10 )
            return 0;
        if ( err > err0 )
            return 3;
        err0 = err;
#if ( 0 )
        printf("\n diff(L,U) = %f", max_diff(ns-1, upe, low+1));
        printf("\n L"); VecPrint::print(std::cout, std::min(16UL, ns-1), low+1, 3);
        printf("\n D"); VecPrint::print(std::cout, std::min(16UL, ns  ), dia, 3);
        printf("\n U"); VecPrint::print(std::cout, std::min(16UL, ns-1), upe, 3);
#endif
#if ( 0 )
        real asy = 0, sup = 0;
        for ( size_t i = 0; i < ns-1; ++i )
        {
            sup = std::max(sup, abs_real(low[i+1]+upe[i]));
            asy += abs_real(low[i+1]-upe[i]);
        }
        printf("\n %3i diff(low-upe) %12.6f", cnt, 2*asy/sup);
#endif
        lapack::xgtsv(ns, 1, low+1, dia, upe, val, ns, &info);
        if ( info )
        {
            std::cerr << " lapack::dgtsv failed " << info << '\n';
            return 2;
        }

        // update result following Newton's iteration
        for ( size_t u = 0; u < ns; ++u )
            sca[u] -= 0.5 * val[u];

        //printf("\n   ->sca "); VecPrint::print(std::cout, ns, sca, 3);
    }
    //printf("\n   >>err %20.16f", err);
    //printf("\n   >>sca "); VecPrint::print(std::cout, ns, sca, 3);
    return 4;
}

/**
 old version 22.02.2015
 externally provided memory `mem[]` should be allocated to hold `5*chk` reals
 */
int Chain::reshape_calculate_alt(const size_t ns, real cutcut,
                                 const real* dif, real* mem, size_t chk)
{
    real * sca = mem;
    real * val = mem+chk;
    real * dia = mem+chk*2;
    real * low = mem+chk*3;
    real * upe = mem+chk*4;

    zero_real(ns, sca);

    real err0 = INFINITY;
    size_t cnt = 0;
    while ( ++cnt < 16 )
    {
        //printf("\n   %i sca  ", cnt); VecPrint::print(std::cout, ns, sca, 3);
        // calculate the matrix elements and RHS of system
        {
            Vector dif0(dif), dif1(dif+DIM);
            Vector vec = (1-2*sca[0]) * dif0 + sca[1]*dif1;
            val[0] = vec.normSqr() - cutcut;
            dia[0] = dot(vec, dif0) * ( -2 );
            upe[0] = dot(vec, dif1);
        }
        real err = abs_real(val[0]);
        for ( size_t i = 1; i+1 < ns; ++i )
        {
            Vector difA(dif+DIM*(i-1)), difB(dif+DIM*i), difC(dif+DIM*(i+1));
            Vector vec = sca[i-1]*difA + (1-2*sca[i])*difB + sca[i+1]*difC;
            val[i] = vec.normSqr() - cutcut;
            low[i] = dot(vec, difA);
            dia[i] = dot(vec, difB) * ( -2 );
            upe[i] = dot(vec, difC);
            err += abs_real(val[i]);
        }
        {
            Vector difA(dif+DIM*(ns-2)), difB(dif+DIM*(ns-1));
            Vector vec = sca[ns-2]*difA + (1-2*sca[ns-1])*difB;
            val[ns-1] = vec.normSqr() - cutcut;
            low[ns-1] = dot(vec, difA);
            dia[ns-1] = dot(vec, difB) * ( -2 );
        }
        err += abs_real(val[ns-1]);
#if ( 0 )
        printf("\n%4lu err %20.16f norm(val) %8.5f", cnt, err, blas::nrm2(ns, val));
        printf("\n     val "); VecPrint::print(std::cout, ns, val, 3);
        printf("\n     sca "); VecPrint::print(std::cout, ns, sca, 3);
#endif
        if ( err < 1e-10 )
            return 0;
        if ( err > err0 )
            return 3;
        err0 = err;
#if ( 0 )
        printf("\n diff(L,U) = %f", max_diff(ns-1, upe, low+1));
        printf("\n L"); VecPrint::print(std::cout, std::min(16UL, ns-1), low+1, 3);
        printf("\n D"); VecPrint::print(std::cout, std::min(16UL, ns  ), dia, 3);
        printf("\n U"); VecPrint::print(std::cout, std::min(16UL, ns-1), upe, 3);
#endif
        int info = 0;
        lapack::xgtsv(ns, 1, low+1, dia, upe, val, ns, &info);
        if ( info )
        {
            std::cerr << " lapack::dgtsv failed " << info << '\n';
            return 1;
        }
        
        // update result following Newton's iteration
        for ( size_t u = 0; u < ns; ++u )
            sca[u] -= 0.5 * val[u];
    }

#if ( 0 )
    printf("\n%2i err %e", cnt, err);
    printf("\n%2i sca  ", cnt); VecPrint::print(std::cout, ns, sca, 3);
    printf("\n");
#endif
    
    return 4;
}



/**
 Apply correction of magnitude 'sca' along the segment directions:
 
 except for the edges, this is:
 P[i] <- P[i] + sca[i] * ( P[i+1] - P[i] ) - sca[i-1] * ( P[i] - P[i-1] )
 
 src[] and dst[] should be arrays of size `DIM*(nbs+1)`
 sca[] should be an array of size `nbs`
 
 The code is similar to projectForcesD(), with an additional difference
 
 */
void Chain::reshape_apply(const size_t nbs, const real* src,
                          const real * sca, real* dst)
{
    assert_true( nbs > 1 );
    Vector B(src);
    Vector nxt(B);

    for ( size_t i = 0; i < nbs; ++i )
    {
        Vector C(src+DIM*(i+1));
        Vector add = sca[i] * ( C - B );
        (nxt+add).store(dst+DIM*i);
        nxt = C - add;
        B = C;
    }
    
    nxt.store(dst+DIM*nbs);
}


/**
 Apply correction (old version)
 */
void Chain::reshape_apply_alt(const size_t nbs, const real* dif, const real* src,
                              const real* sca, real* dst)
{
    assert_true( nbs > 1 );
    real nxt[DIM];
    
    for ( int d = 0; d < DIM; ++d )
        nxt[d] = src[d];
    
    for ( size_t i = 0; i < nbs; ++i )
    {
        for ( int d = 0; d < DIM; ++d )
        {
            real x = sca[i] * dif[DIM*i+d];
            dst[DIM*i+d] = nxt[d] + x;
            nxt[d] = src[DIM*(i+1)+d] - x;
        }
    }
    
    for ( int d = 0; d < DIM; ++d )
        dst[DIM*nbs+d] = nxt[d];
}


/*
 Try to re-establish the length of the segments, by moving points along
 the directions of the flanking segments
 
 here ns = nbSegments() and tmp_size >= ns
 `tmp` should be of size (5+3+3)*tmp_size
 */
int Chain::reshape_local(const size_t nbs, const real* src, real* dst,
                         real cut, real* tmp, size_t tmp_size)
{
    int res;
    assert_true( nbs > 1 );
    //auto rdt = __rdtsc();

    // alias some of the memory provided for the work:
    real * mag = tmp;
    real * pri = tmp + tmp_size;
    real * sec = tmp + tmp_size * 2;
    real * mem = tmp + tmp_size * 3;
    
    // calculate terms using 'dif[]':
    Vector A, B, C;
    A.load_diff(src);
    B.load_diff(src+DIM);
    mag[0] = A.normSqr();
    mag[1] = B.normSqr();
    pri[0] = dot(A, B);
    for ( size_t i = 2; i < nbs; ++i )
    {
        C.load_diff(src+DIM*i);
        mag[i] = C.normSqr();
        pri[i-1] = dot(B, C);
        sec[i-2] = dot(A, C);
        A = B;
        B = C;
    }
    // these terms should not be used:
    pri[nbs-1] = 0;
    sec[nbs-2] = 0;
    sec[nbs-1] = 0;
    /*
    printf("\n mag "); VecPrint::print(std::cout, nbs, mag, 3);
    printf("\n pri "); VecPrint::print(std::cout, nbs, pri, 3);
    printf("\n sec "); VecPrint::print(std::cout, nbs, sec, 3);
    */
    res = reshape_calculate(nbs, cut*cut, mag, pri, sec, mem, tmp_size);
    
#if ( 0 )
    real * dif = new_real(tmp_size*DIM);
    // calculate differences
    for ( size_t p = 0; p < DIM*nbs; ++p )
        dif[p] = src[p+DIM] - src[p];
    for ( size_t p = 0; p < DIM; ++p )
        dif[DIM*nbs+p] = 0;

    // checking against older code
    copy_real(nbs, mem, mag);
    res = reshape_calculate_alt(nbs, cut*cut, dif, mem, tmp_size);
    real errS = blas::max_diff(nbs, mem, mag);
    if ( abs_real(errS) > 1e-8 )
    {
        printf("\n Chain:err scalar %20.16f", errS);
        printf("\n mul "); VecPrint::print(std::cout, nbs, mag, 3);
        printf("\nLmul "); VecPrint::print(std::cout, nbs, mem, 3);
    }
#endif
    //printf("reshape %3i cycles: %16llu\n", ns, __rdtsc()-rdt);

    if ( res == 0 )
    {
        reshape_apply(nbs, src, mem, dst);
#if ( 0 )
        // checking against older code
        reshape_apply_alt(nbs, dif, src, mem, mag);
        //projectForcesD_(nbs, dif, src, mem, mag);  // similar calculation!
        size_t nbv = DIM*nbs+DIM;
        real errP = blas::max_diff(nbv, mag, dst);
        //printf("\n Chain:reshape errors %6lu %20.16f %20.16f", nbs, errS, errP);
        if ( abs_real(errP) > 1e-8 )
        {
            printf("\n Chain:reshape errors %6lu %20.16f %20.16f", nbs, errS, errP);
            printf("\n pts "); VecPrint::print(std::cout, nbv, dst, 2);
            printf("\nLpts "); VecPrint::print(std::cout, nbv, mag, 2);
        }
        free_real(dif);
#endif
    }
    return res;
}


/**
 The response of this method to a sudden perpendicular force is not ideal:
 For example, a force applied to the bottom of a vertical fibers leads
 to a 'L' configuration after one step of `solve()`.
 reshape() reduces the bottom leg of the 'L', by translating the entire vertical portion
 of the fiber, irrespective of the length of this section.
 */

#if ( 1 )   // 1 = optimized version of Chain::reshape_global()

/**
 Move the vertices relative to each other, such that when this is done,
 all segments have the same distance `fnCut` ( =segmentation() ).
 This is operation does not change the center of gravity of the fiber.

 
 NOTE: if two consecutive points overlap, there is no unique way to
 restore the constraints! We do nothing in that case, because most 
 likely, the Brownian motion will push the points appart soon.
 */

void Chain::reshape_global(const size_t ns, const real* src, real* dst, real cut)
{
    Vector inc(0,0,0), sum(0,0,0), seg;
    seg.load_diff(src);
    real   dis = seg.norm();
    
    // translation needed to restore first segment
    if ( dis > REAL_EPSILON )
        inc = ( cut/dis - 1.0 ) * seg;
    
    Vector(src).store(dst);

    for ( size_t i = 1; i < ns; ++i )
    {
        seg.load_diff(src+DIM*i);
        dis = seg.norm();
        
        //move the left point by off:
        (inc+Vector(src+DIM*i)).store(dst+DIM*i);
        //update the uniform motion of the points:
        sum += inc;
        
        //add to the translation needed to restore this segment
        if ( dis > REAL_EPSILON )
            inc += ( cut/dis - 1.0 ) * seg;
    }
    
    // move the last point by dp:
    (inc+Vector(src+DIM*ns)).store(dst+DIM*ns);
    
    // calculate uniform motion needed to conserve the center of gravity:
    sum = ( sum + inc ) / ( ns + 1 );
    
    // translate entire fiber uniformly:
    for ( size_t i = 0; i <= ns; ++i )
        sum.sub_to(dst+DIM*i);
}

#else

// ------------  old ( less optimal ) version:
/**
 Move the vertices relative to each other, such that when this is done,
 all segments have the same distance segmentation() = fnCut.
 This is operation does not change the center of gravity of the fiber.
 */

void Chain::reshape_global(const size_t ns, const real* src, real* dst, real cut)
{
    Vector off, sum(0,0,0);

    copy_real(DIM*(ns+1), src, dst);
    for ( size_t pp = 1; pp <= ns; ++pp )
    {
        off.load_diff(dst+DIM*(pp-1));
        real dis = off.norm();
        if ( dis > REAL_EPSILON )
        {
            off  *= ( cut/dis - 1.0 );
            for ( size_t qq = pp; qq <= ns; ++qq )
                off.add_to(dst+DIM*qq);
            sum += ( 1 + ns - pp ) * off;
        }
    }
    
    sum /= ns + 1;
    for ( size_t i = 0; i <= ns; ++i )
        sum.sub_to(dst+DIM*i);
}

#endif

/**
 Replace coordinates by the ones provided in `ptr`
 A reshape operation is done
 */
void Chain::getPoints(real const* ptr)
{
#if NEW_SKIP_PROJECTION
    if ( skipProjection )
        return copy_real(DIM*nPoints, ptr, pPos);
#endif
#if 0
    // use here static memory
    // Attention: this only works if using a single thread
    static size_t alc = 0;
    static real* mem = nullptr;
    
    if ( alc < allocated() )
    {
        alc = allocated();
        free_real(mem);
        mem = new_real(alc*8);
    }
#else
    // use here thread-local static memory
    thread_local static size_t alc = 0;
    
    auto delete_real = [](real * x)
    {
        //printf("> del real[%lu] %p %p\n", alc, x, pthread_self());
        free_real(x);
        alc = 0;
    };

    thread_local static std::unique_ptr<real, decltype(delete_real)> uptr(nullptr, delete_real);
    
    if ( alc < allocated() || uptr.get() == nullptr )
    {
        alc = allocated();
        free_real(uptr.release());
        uptr.reset(new_real(alc*8));
        //printf("> new real[%lu] %p %p\n", alc, uptr.get(), pthread_self());
    }
    real * mem = uptr.get();
#endif
    
#if ( DIM > 1 )
    if ( nPoints == 2 )
        reshape_two(ptr, pPos, fnCut);
    else if ( reshape_local(nbSegments(), ptr, pPos, fnCut, mem, allocated()) )
#endif
    {
        reshape_global(nbSegments(), ptr, pPos, fnCut);
        Cytosim::warn << "a crude method was used to reshape " << reference() << '\n';
    }
    //dump(std::cerr);
    //iDirValid = false;
}


/**
 Flip all the points, such that minus_end becomes plus_end and vice-versa.
 This does not affects Abscissa and the abscissa of center thus stays as it is.
*/
void Chain::flipChainPolarity()
{
    size_t ii = 0;
    size_t jj = lastPoint();
    
    while ( ii < jj )
    {
        Vector P(pPos+DIM*ii);
        Vector Q(pPos+DIM*jj);
        Q.store(pPos+DIM*ii);
        P.store(pPos+DIM*jj);
        ++ii;
        --jj;
    }
}


//========================================================================
//=====================GROWING/SHRINKING==================================
//========================================================================
#pragma mark -

/**
 The argument 'delta' can be positive or negative:
 - delta > 0 : elongation,
 - delta < 0 : shortening
 .
 
 Note 1: This works nicely if `delta` is small compared to segmentation().
 For large decrease in length, use cutM().
 
 Note 2: Unless the Chain is straight, the length of the segments after this
 will not exactly match `segmentation()`.
*/
void Chain::growM(const real delta)
{
    assert_true( length() + delta > REAL_EPSILON );
    real a = -delta / length();
    
    if ( delta > 0 )
    {
        size_t p = 0, n = nbSegments();
        Vector dp0 = diffPoints(0), dp1;
        movePoint(p, ( a * n ) * dp0);
        ++p;
        --n;
        
        if ( n > 0  &&  ( n & 1 ) )
        {
            dp1 = diffPoints(p);
            movePoint(p, ( a * n ) * dp0);
            dp0 = dp1;
            ++p;
            --n;
        }
        
        while ( n > 1 )
        {
            //assert_true( 0 == (p & 1) );
            dp1 = diffPoints(p);
            movePoint(p, ( a * n ) * dp0);
            ++p; --n;
            //assert_true( 1 == (p & 1) );
            dp0 = diffPoints(p);
            movePoint(p, ( a * n ) * dp1);
            ++p; --n;
        }
    }
    else if ( delta < 0 )
    {
        for ( size_t p = 0, n = nbSegments(); n > 0; ++p, --n )
            movePoint(p, ( a * n ) * diffPoints(p));
    }
    
    fnAbscissaM -= delta;
    setSegmentation(std::max(fnCut+delta/nbSegments(), REAL_EPSILON));
    postUpdate();
}

/**
 This extends the fiber by adding one segment at the MINUS_END.
 Thus `segmentation()` is not changed, and the existing points are not displaced.
 */
void Chain::addSegmentM()
{
    size_t pp = 1+nPoints;
    setNbPoints(pp);
    
    pp *= DIM;
    while ( --pp >= DIM )
        pPos[pp] = pPos[pp-DIM];
    
    for ( pp = 0; pp < DIM; ++pp )
        pPos[pp] += pPos[pp] - pPos[pp+2*DIM];
    
    fnAbscissaM -= fnCut;
    postUpdate();
}


/**
 The Fiber length is reduced by `delta` ( which must be >= 0 ).
 The portion of size `delta` near the MINUS_END is removed,
 the (fewer) vertices are recalculated.
 
 Note: after cutM(), the distance between the points is not exactly
 equal to segmentation(). This is true only if the fiber is straight.
 */
void Chain::cutM(const real delta)
{
    real len = length();
    assert_true( 0 <= delta );
    assert_true( delta < len );
    
    const size_t np = bestNumberOfPoints((len-delta)/fnSegmentation);
    const real cut = (len-delta) / (np-1);
    real* tmp = new_real(DIM*np);

    // calculate intermediate points:
    for ( size_t i=0; i+1 < np; ++i )
    {
        Vector w = interpolateM(delta+i*cut).pos();
        w.store(tmp+DIM*i);
    }

    // copy the position of plus-end:
    copy_real(DIM, pPos+DIM*lastPoint(), tmp+DIM*(np-1));
    
    setNbPoints(np);
    fnAbscissaM += delta;
    setSegmentation(cut);
    getPoints(tmp);
    free_real(tmp);
    postUpdate();
}


/**
 The argument 'delta' can be positive or negative:
 - delta > 0 : elongation,
 - delta < 0 : shortening
 .
 
 Note 1: This works nicely if `delta` is small compared to segmentation().
 For large decrease in length, use cutP().

 Note 2: Unless the Chain is straight, the length of the segments after this
 will not exactly match `segmentation()`.
 */
void Chain::growP(const real delta)
{
    assert_true( length() + delta > REAL_EPSILON );
    real a = delta / length();
    
    if ( delta > 0 )
    {
        size_t p = lastPoint();
        Vector dp0 = diffPoints(p-1), dp1;
        movePoint(p, ( a * p ) * dp0);
        --p;
        
        if ( p > 0  &&  ( p & 1 ) )
        {
            dp1 = diffPoints(p-1);
            movePoint(p, ( a * p ) * dp0);
            dp0 = dp1;
            --p;
        }
        
        while ( p > 1 )
        {
            //assert_true( 0 == (p & 1) );
            dp1 = diffPoints(p-1);
            movePoint(p, ( a * p ) * dp0);
            --p;
            //assert_true( 1 == (p & 1) );
            dp0 = diffPoints(p-1);
            movePoint(p, ( a * p ) * dp1);
            --p;
        }
    }
    else if ( delta < 0 )
    {
        for ( size_t p = lastPoint() ; p > 0 ; --p )
            movePoint(p, ( a * p ) * diffPoints(p-1));
    }
    
    setSegmentation(std::max(fnCut+delta/nbSegments(), REAL_EPSILON));
    fnAbscissaP += delta;
    postUpdate();
}


/**
 This extends the fiber by adding one segment at the PLUS_END.
 Thus `segmentation()` is not changed, and the existing points are not displaced.
 */
void Chain::addSegmentP()
{
    size_t pp = nPoints;
    setNbPoints(pp+1);
    
    real * psp = pPos + DIM * pp;
    for ( size_t dd = 0; dd < DIM; ++dd )
        psp[dd] = 2 * psp[dd-DIM] - psp[dd-2*DIM];
    
    fnAbscissaP += fnCut;
    postUpdate();
}


/**
 The Fiber length is reduced by `delta` ( which must be >= 0 ).
 The portion of size `delta` near the PLUS_END is removed,
 and the fewer vertices are recalculated.

 Note: after cutP(), the distance between the points is not exactly
 equal to segmentation(). This is true only if the fiber is straight.
*/
void Chain::cutP(const real delta)
{
    real len = length();
    assert_true( 0 <= delta );
    assert_true( delta < len );
    
    const size_t np = bestNumberOfPoints((len-delta)/fnSegmentation);
    const real cut = (len-delta) / (np-1);
    real* tmp = new_real(DIM*np);

    // copy minus end:
    copy_real(DIM, pPos, tmp);

    // calculate intermediate points:
    for ( size_t i = 1; i < np; ++i )
    {
        Vector w = interpolateM(i*cut).pos();
        w.store(tmp+DIM*i);
    }

    setNbPoints(np);
    setSegmentation(cut);
    fnAbscissaP -= delta;
    getPoints(tmp);
    free_real(tmp);
    postUpdate();
}

//------------------------------------------------------------------------------

void Chain::grow(FiberEnd end, const real delta)
{
    if ( end == PLUS_END )
        growP(delta);
    else if ( end == MINUS_END )
        growM(delta);
}


void Chain::adjustLength(real len, FiberEnd ref)
{
    assert_true( len > 0 );
    
    if ( ref == PLUS_END )
    {
        if ( len < length() )
            cutP(length()-len);
        else
            growP(len-length());
    }
    else if ( ref == MINUS_END )
    {
        if ( len < length() )
            cutM(length()-len);
        else
            growM(len-length());
    }
}


void Chain::truncateM(size_t p)
{
    Mecable::truncateM(p);
    fnAbscissaM = abscissaPoint(p);
    postUpdate();
}


void Chain::truncateP(size_t p)
{
    Mecable::truncateP(p);
    fnAbscissaP = abscissaPoint(p);
    postUpdate();
}


/**
 `fib` is added at the PLUS_END of `*this`
 
 The vertex are reinterpolated linearly, and the length of the
 segments will not fullfil exactly the constraints of segmentation.
 If this is a problem, Chain::reshape() should be called.
 
 `fib` should usually be destroyed afterward.
 */
void Chain::join(Chain const* fib)
{
    const real len1 = length();
    const real lenT = len1 + fib->length();
    const size_t ns = bestNumberOfPoints(lenT/fnSegmentation) - 1;
    const real cut = lenT / real(ns);
    
    real* tmp = new_real(DIM*(ns+1));

    // calculate new points into tmp[]:
    for ( size_t i = 1; i < ns; ++i )
    {
        Vector w;
        if ( i*cut < len1 )
            w = interpolateM(i*cut).pos();
        else
            w = fib->interpolateM(i*cut-len1).pos();
        
        w.store(tmp+DIM*i);
    }
    
    // copy position of PLUS_END:
    fib->posEndP().store(tmp+DIM*ns);

    setNbPoints(ns+1);
    setSegmentation(cut);
    fnAbscissaP = fnAbscissaM + cut * fnCut;
    getPoints(tmp);
    free_real(tmp);
    updateFiber();
}

//------------------------------------------------------------------------------
#pragma mark -


/**
 Returns the minimum and maximum distance between consecutive points
 */
void Chain::segmentationMinMax(real& mn, real& mx) const
{
    mn = diffPoints(0).normSqr();
    mx = mn;
    for ( size_t n = 1; n < lastPoint(); ++n )
    {
        real r = diffPoints(n).normSqr();
        mn = std::min(mn, r);
        mx = std::max(mx, r);
    }
    mn = sqrt(mn);
    mx = sqrt(mx);
}

/**
 Returns the average and variances of segment length
 */
void Chain::segmentationVariance(real& avg, real& var) const
{
    avg = 0;
    var = 0;
    size_t cnt = nbSegments();
    for ( size_t n = 0; n < cnt; ++n )
    {
        real r = diffPoints(n).norm();
        avg += r;
        var += r*r;
    }
    avg /= (real)cnt;
    var = var / (real)cnt - avg * avg;
}

/**
 Calculate the inverse of the radius of the circle containing the points A, B, C
 
     cos(angle) = scalar_product( AB, BC ) / ( |AB| * |BC| )
     sin(angle) = sqrt( 1 - cos(angle)^2 )
     2 * radius * sin(angle) = |AC|
     curvature = 2 * sin(angle) / |AC|
     curvature = 2 * sqrt( ( 1 - cos(angle)^2 ) / AC^2 )

 curvature is ZERO if A, B and C are aligned
 */
real curvature3(Vector const& A, Vector const& B, Vector const& C)
{
    Vector ab = B - A;
    Vector bc = C - B;
    real P = dot(ab, bc);
    real S = std::max(0.0, 1.0 - ( P * P ) / ( ab.normSqr() * bc.normSqr() ));
    real D = ( C - A ).normSqr();
    return 2.0 * sqrt( S / D );
}


real Chain::curvature(size_t p) const
{
    if ( p < 1 || lastPoint() <= p )
        return 0;
    return curvature3(posP(p-1), posP(p), posP(p+1));
}


/**
 The normalized bending energy is an integral over the curvilinear abscissa `s`:
 
     1/2 * sum( curvature(s)^2 ds )
 
 The curvature is calculated from the positions of the vertices:
 Given theta = angle between two consecutive segments,

     curvature = 1/R = 2 * sin(angle/2) / segmentation
 
 and since
 
     sin^2(angle/2) = ( 1 - cos(angle) ) / 2
 
 hence:
 
     1/2 * curvature^2 = ( 1 - cos(angle) ) / ( segmentation^2 )
 
 and finaly:
 
     1/2 * sum( curvature^2 * ds ) = sum( 1 - cos(angle) ) / segmentation
 
 */
real Chain::bendingEnergy0() const
{
    real e = 0;
    
    const size_t lsp = nPoints - 2;
    if ( lsp > 0 )
    {
        for ( size_t p = 0; p < lsp ; ++p )
        {
            Vector A = posP(p);
            Vector B = posP(p+1);
            Vector C = posP(p+2);
            e += dot(B - A, C - B);  // e += cos(angle) * segmentation^2
        }
        // e <- sum( 1 - cos(angle) )
        e = lsp - e / ( fnCut * fnCut );
        
        /*
         We correct the result, because we only considered (nPoints-2) junctions,
         and thus covered only a fraction of the total length of the filament
         */
        e *= ( lsp + 1 ) / ( fnCut * lsp );
    }
    
    return e;
}


real Chain::minCosinus() const
{
    real result;
    Vector dir1, dir2;
    
    size_t ps = nbSegments() % 2;
    if ( ps )
    {
        dir1   = diffPoints(0);
        result = fnCut * fnCut;
    }
    else
    {
        dir1   = diffPoints(1);
        result = dot(diffPoints(0), dir1);
        ps = 2;
    }
    
    for ( ; ps < nbSegments(); ps+=2 )
    {
        dir2 = diffPoints(ps);
        real s = dot(dir1, dir2);
        if ( s < result ) result = s;
        dir1 = diffPoints(ps+1);
        real t = dot(dir1, dir2);
        if ( t < result ) result = t;
    }
    
    return result / ( fnCut * fnCut );
}


/**
 Returns the minimum and maximum distance between consecutive points
 */
size_t Chain::nbKinks(real threshold) const
{
    size_t res = 0;
    threshold *= fnCut * fnCut;
    Vector d = diffPoints(0);
    
    for ( size_t n = 1; n < lastPoint(); ++n )
    {
        Vector r = diffPoints(n);
        if ( dot(d, r) < threshold )
            ++res;
        d = r;
    }
    return res;
}

/**
 This calculates the intersection between the support line of segment `s`,
 and the plane defined by <em> n.pos + a = 0 </em>

 @return scalar `x` specifying the intersection with the support line:
 - `x = 0` if intersection occurs at point 's'
 - `x in ]0, 1[` for intersections that are within the segment boundaries
 - `x = 1` if intersection occurs at point 's+1'
 - `x = INFINITY` if the segment is parallel to the plane
 .
 
 The abscissa of the intersection is `abscissaPoint(s+a)`.
 The position of the cut is `posPoint(s, a)`
 */

real Chain::planarIntersect(size_t s, Vector const& n, const real a) const
{
    assert_true( s < nbSegments() );
    
    real sca = dot(diffPoints(s), n);
    
    // if segment is parallel to plane, there is no intersection:
    if ( abs_real(sca) < REAL_EPSILON )
        return INFINITY;
    
    Vector pos = posP(s);
    
    if ( modulo )
        modulo->fold(pos);
    
    return - ( dot(pos, n) + a ) / sca;
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 Recalculate the vertices of the Chain for 'ns' segments.
 
 @todo 2d-order interpolation in Chain::resegment()
 
 Note: Unless the Chain is straight, the length of the segments after
 an interpolation will not exactly match `segmentation()`, and we therefore
 call reshape() here to correct for the problem.
 */
void Chain::resegment(size_t ns)
{
    //std::clog << reference() << " resegment " << ns << "\n";
    assert_true( ns > 0 );
    real cut = nbSegments() * segmentation() / ns;
    
    // calculate new intermediate points in tmp[]:
    Vector a = posP(0), b = posP(1);
    
    real h = 0;
    size_t p = 1;
    real* tmp = new_real(DIM*(ns+2));
    
    a.store(tmp);
    for ( size_t n = 1; n < ns; ++n )
    {
        h += cut;
        
        while ( h > fnCut )
        {
            h -= fnCut;
            a = b;
            ++p;
            assert_true(p<nbPoints());
            b.load(pPos+DIM*p);
        }
        
        Vector w = a + ( h / fnCut ) * ( b - a );
        w.store(tmp+DIM*n);
    }
    
    // copy coordinates of last point:
    a.load(pPos+DIM*lastPoint());
    a.store(tmp+DIM*ns);

    // resize filament:
    setNbPoints(ns+1);
    setSegmentation(cut);
    getPoints(tmp);
    free_real(tmp);
}


#if !CURVATURE_DEPENDENT_SEGMENTATION

/**
 A fiber is segmented as a function of its length.
 The number of segments `NS` is the one that minimizes the absolute value:

     | length / NS - segmentation |
 
 Where `segmentation` is the parameter. NS is such that:

     length / NS < 4/3 * segmentation
     length / NS > 2/3 * segmentation

 */
void Chain::adjustSegmentation()
{
    assert_true( fnSegmentation > REAL_EPSILON );
    
    size_t best = bestNumberOfPoints(nbSegments()*fnCut/fnSegmentation);
    
    if ( best != nPoints )
    {
        //std::clog << reference() << " resegment " << nPoints << " -> " << best << "\n";
#if ( 1 )
        resegment(best-1);
#else
        // copy current points in temporary array:
        real* tmp = new_real(DIM*allocated());
        copy_real(DIM*nPoints, pPos, tmp);
        // re-interpolate:
        setShape(tmp, nPoints, best);
        free_real(tmp);
#endif
    }
}


#else

static constexpr real RECUT_PRECISION = 0.05;

/**
 A fiber is segmented as a function of its curvature.
 Only one new segmentation is done every time step at maximum,
 to allow the solver to equilibrate the new vertices.
 */
void Chain::adjustSegmentation()
{
    LOG_ONCE("adjustSegmentation() is using the fiber curvature\n");
    
    const int upLimit = 8;
    size_t nbs = nbSegments();
    real len = length();
    
    assert_true( fnSegmentation > REAL_EPSILON );
    // use one segment for very short tubes
    if ( len <= fnSegmentation )
    {
        if ( nbs > 1 )
        {
            resegment(1);
            goto clear_counter;
        }
        return;
    }
    //use at least two segments if length exceeds 2*fnSegmentation
    if ( nbs < 2 )
    {
        if ( len > fnSegmentation )
        {
            resegment(2);
            goto clear_counter;
        }
        return;
    }
    //the segment-length should not exceed a upper limit
    if ( fnCut >= upLimit * fnSegmentation )
    {
        resegment(2*nbs);
        goto clear_counter;
    }
    
    // accumulate the error in variable autoCut, 
    // The error is in angle-square: 1-cos(angle) ~ (angle)^2
    autoCutCnt += 1.0;
    autoCutVal += minCosinus();
    
    // after accumulation of RECUT_PERIOD time steps, we check the result
    if ( ++autoCutIdx < RECUT_PERIOD )
        return;
    
    if ( autoCutCnt > 2 )
    {
        // calculate error accumulated over the last steps
        autoCutVal = ( 1.0 - autoCutVal/autoCutCnt ) * fnCut;
        
        if ( autoCutVal > RECUT_PRECISION )
        {
            //cut more finely if the error is large
            if ( fnCut > fnSegmentation )
                resegment(2*nbs);
        }
        else if ( autoCutVal < 0.1 * RECUT_PRECISION )
        {
            /*
             We want the future error with fewer points to be small.
             At constant curvature, the error scales like (length-of-rods)^2
             We can expect the error to raise by 4x, if we divide the number of segments by 2.
             For safety, we use 0.125 here, instead of 1/4, i.e. add another factor 2.
             */
            if ( nbs > 3 && fnCut < (upLimit/2)*fnSegmentation )
                resegment(nbs/2);
        }
    }
clear_counter:
    //reset the counter and error accumulator
    autoCutIdx = 0;
    autoCutCnt = 0;
    autoCutVal = 0;
}

#endif

//------------------------------------------------------------------------------
#pragma mark -

/**
 return the abscissa with respect to the ORIGIN.
 */
real Chain::abscissaEnd(const FiberEnd end) const
{
    switch( end )
    {
        case ORIGIN:    return 0;
        case PLUS_END:  return abscissaP();
        case MINUS_END: return abscissaM();
        case CENTER:    return abscissaC();
        default:        ABORT_NOW("invalid argument value"); return 0;
    }
}


/**
 returns the abscissa (from the ORIGIN) of a point that is specified
 by a distance from the given reference.
 */
real Chain::abscissaFrom(const real dis, const FiberEnd ref) const
{
    switch( ref )
    {
        case ORIGIN:     return dis;
        case PLUS_END:   return abscissaP() - dis;
        case MINUS_END:  return dis + abscissaM();
        case CENTER:     return dis + abscissaC();
        default:         ABORT_NOW("invalid argument value"); return 0;
    }
}

/**
 This uses values at [1], [2] and [3] of opt[key] to define an abscissa
 on the Fiber:
 
     attach = FIBER, ABSCISSA, REFERENCE, MODIFIER
 
 with
 
     ABSCISSA = REAL
     REFERENCE = { 'plus_end', 'minus_end', 'center' }  (default = 'origin')
     MODIFIER = { 'none', 'uniform', 'exponential' }  (default = 'none')
 
 All these parameters are optional.
 The abscissa is counted from the reference and towards the other end.
 The MODIFIER introduces a random component to the position.
 If ABSCISSA is not specified, this return a random abscissa.
 
 Example:
 
     new filament
     {
         attach1 = simplex, 0.0, minus_end
         attach2 = simplex, 0.0, plus_end
     }

*/
real Chain::someAbscissa(real dis, FiberEnd ref, int mod, real alpha) const
{
    const real len = length();
    real a = dis;
    
    switch ( mod )
    {
        case 0:
            break;
        case 1:  // random
            do {
                a = dis * RNG.preal();
            } while ( a > len );
            break;
        case 2:  // exponential
            do {
                a = dis * RNG.exponential();
            } while ( a > len );
            break;
        case 3:  // regular
            a *= alpha;
            break;
        case 7:
            return RNG.real_uniform(abscissaM(), abscissaP());
    }
    
    dis = abscissaFrom(a, ref);

    if ( !betweenMP(dis) )
        throw InvalidParameter("hand::abscissa is out of range");
        
    return dis;
}

/**
 The Fiber is partitionned by this function in three regions:
 - a MINUS_END part of length `lambda`
 - a PLUS_END part, also of length `lambda`
 - and a NO_END section in between
 .
 Note that a Fiber shorter than `2*lambda` does not have a central region,
 and is composed of PLUS_END and MINUS_END parts of equal size.
 */    
FiberEnd Chain::whichEndDomain(const real ab, const real lambda) const
{
    const real abs = ab - fnAbscissaM;
    const real len = length();
    
    if ( 2 * abs > len )
    {
        if ( abs >= len - lambda )
            return PLUS_END;
    }
    else
    {
        if ( abs <= lambda )
            return MINUS_END;
    }
    return NO_END;
}


//------------------------------------------------------------------------------
#pragma mark -


Mecapoint Chain::exactEnd(const FiberEnd end) const
{
    if ( end == MINUS_END )
        return Mecapoint(this, 0);
    else
    {
        assert_true( end == PLUS_END );
        return Mecapoint(this, lastPoint());
    }
}


Interpolation Chain::interpolateEnd(const FiberEnd end) const
{
    switch( end )
    {
        case MINUS_END:
            return Interpolation(this,0, 1, 0);
        case PLUS_END:
            return Interpolation(this, nPoints-2, nPoints-1, 1);
        case CENTER:
            return interpolateCenter();
        default:
            throw Exception("unexpected argument");
    }
}


Interpolation Chain::interpolateCenter() const
{
    size_t n = lastPoint() / 2;
    return Interpolation(this, n, n+1, 0.5 * ( lastPoint() - 2*n ));
}


/**
 return Interpolation corresponding to a distance `ab` from the MINUS_END
 The interpolation describes a position:
       X = P(r) * (1-a) + P(r+1) * a
 where
 - `r` is an integer: 0 <= r < lastPoint(),
 - `a` is a positive real coefficient: 0 <= a <= 1
 .
 When `ab` is above the PLUS_END, an interpolation of the last point is returned.
 
 */
Interpolation Chain::interpolateM(const real ab) const
{
    real a = std::max(ab, (real)0) / fnCut;
    //beyond the last point, we interpolate the PLUS_END
    size_t s = std::min((size_t)a, nPoints-2);
    return Interpolation(this, s, s+1, std::min(a-s, (real)1));
}


Interpolation Chain::interpolate(const real ab, const FiberEnd end) const
{
    switch( end )
    {
        case ORIGIN:
            return interpolate(ab);
            
        case MINUS_END:
            return interpolateM(ab);
            
        case CENTER:
            return interpolateM(ab + 0.5*length());
            
        case PLUS_END:  //this is counted from the plus towards the minus end
            return interpolateM(fnCut*nbSegments() - ab);
        
        default:
            ABORT_NOW("invalid argument value");
    }
    return interpolate(0);
}

//------------------------------------------------------------------------------
#pragma mark -

#if ( DIM > 1 )
Vector Chain::posM(const real ab) const
{
    // return MINUS_END
    if ( ab <= 0 )
        return posP(0);
    
    real a = ab / fnCut;
    size_t s = (size_t)a;
    
    // check if PLUS_END is reached:
    if ( s+1 < nPoints )
        return posPoint(s, a-s);
    else
        return posP(lastPoint());
}

Vector Chain::dirM(const real ab) const
{
    // at MINUS_END
    if ( ab <= 0 )
        return dirSegment(0);
    
    real a = ab / fnCut;
    size_t s = (size_t)a;
    
    // check if PLUS_END is reached
    if ( s+1 < nPoints )
        return dirSegment(s);
    else
        return dirSegment(lastSegment());
}
#endif


Vector Chain::posEnd(FiberEnd end) const
{
    if ( end == MINUS_END )
        return posEndM();
    else if ( end == PLUS_END )
        return posEndP();
    else
        return posM(abscissaFrom(0, end));
}


Vector Chain::dirEnd(const FiberEnd end) const
{
    if ( end == MINUS_END )
        return dirSegment(0);
    else if ( end == PLUS_END )
        return dirSegment(lastSegment());
    else
        return dirM(abscissaFrom(0, end));
}


/// force on the PLUS_END projected on the direction of elongation
real Chain::projectedForceEndM() const
{
    return -dot(netForce(0), dirSegment(0));
}

/// force on the PLUS_END projected on the direction of elongation
real Chain::projectedForceEndP() const
{
    size_t p = lastSegment();
    return dot(netForce(p+1), dirSegment(p));
}


/**
 The returned value is negative when the force antagonizes elongation,
 and this is true at both ends. 
 */
real Chain::projectedForceEnd(const FiberEnd end) const
{
    if ( end == PLUS_END )
        return projectedForceEndP();
    else
    {
        assert_true( end == MINUS_END );
        return projectedForceEndM();
    }
}


//------------------------------------------------------------------------------
#pragma mark -

int Chain::check(std::ostream& os, real len) const
{
    int res = 0;
    real mn, mx;
    segmentationMinMax(mn, mx);
    real dev = ( mx - mn ) / segmentation();
    real con = contourLength(pPos, nPoints);
    res = ( dev > 0.01 ) + ( abs_real( con - len ) > 0.1 );
    if ( res )
    {
        os << "chain " << std::setw(7) << reference() << '\n';
        os << "{\n";
        os << "    segmentation = " << segmentation() << '\n';
        os << "    segments_min = " << mn << '\n';
        os << "    segments_max = " << mx << '\n';
        os << "    length  = " << len << '\n';
        os << "    contour = " << con << '\n';
        os << "}" << std::endl;
    }
    return res;
}


/**
 Prints info on the length of Segments, which can be useful for debugging
 */
void Chain::dump(std::ostream& os) const
{
    real mn, mx;
    real len = length();
    segmentationMinMax(mn, mx);
    real con = contourLength(pPos, nPoints);

    os << "chain " << std::setw(7) << reference() << '\n';
    os << "{\n";
    os << "    segmentation = " << segmentation() << '\n';
    os << "    segment_min = " << mn << '\n';
    os << "    segment_max = " << mx << '\n';
    os << "    length  = " << len << '\n';
    os << "    contour = " << con << '\n';
    os << "}" << std::endl;
}


void Chain::write(Outputter& out) const
{
    assert_small( length1() - length() );
    out.writeUInt32(signature());
    out.writeFloat(length());
    out.writeFloat(fnSegmentation);
    out.writeFloat(fnAbscissaM);
    out.writeFloat(fnBirthTime);
    Mecable::write(out);
}


/**
 The fiber will be re-segmented if its current desired segmentation 
 does not match the one stored in the file.
 */
void Chain::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    //Cytosim::log << "  reading Chain at " << in.pos() << '\n';
    
    ObjectSignature s = in.readUInt32();
    if ( s ) signature(s);
    
    real len    = in.readFloat();
    real seg    = in.readFloat();
    fnAbscissaM = in.readFloat();
    
#ifdef BACKWARD_COMPATIBILITY
    if ( in.formatID() > 49 ) // 12.12.2018 moved birthTime
#endif
        fnBirthTime = in.readFloat();

    if ( len <= 0 )
        throw InvalidIO("invalid (negative) fiber length");

    if ( len > 1e6 )
        throw InvalidIO("excessive fiber length");
    
    if ( seg <= 1e-6 || seg > 1e6 )
        throw InvalidIO("invalid fiber segmentation");

    Mecable::read(in, sim, tag);
    
    if ( nPoints < 2 )
        throw InvalidIO("invalid fiber with 0 or 1 point");

#ifdef BACKWARD_COMPATIBILITY
    if ( in.formatID() <= 37 )
    {
        setSegmentation(len);
        len *= nbSegments();
    }
    else
#endif
        setSegmentation(len/nbSegments());

    fnAbscissaP = fnAbscissaM + len;

    // resegment if the parameter has changed:
    if ( fnSegmentation != seg )
        adjustSegmentation();
    
    //Mecable::write(std::cerr);
    
    // verify the length and segmentation:
    if ( in.vectorSize() == DIM )
        check(std::clog, len);
}

