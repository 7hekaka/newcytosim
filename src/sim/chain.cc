// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include <cmath>
#include <memory>
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


// defined in "mecafil_project.cc"
//extern void projectForcesD_(index_t, const real*, const real*, const real*, real*);


/**
 This returns N+1, where N is the integer that minimizes
     abs_real( length / N - segmentation ),
 */
index_t Chain::bestNumberOfPoints(const real ratio)
{
    real n = 1 + std::floor(ratio);
    return (index_t)( n + ( ratio*(n-0.5) >= n*(n-1) ));
}


real Chain::contourLength(const real pts[], size_t n_pts)
{
    real len = 0;
    Vector a(pts), b;
    for ( index_t n = 1; n < n_pts; ++n )
    {
        b.load(pts+DIM*n);
        len += (b-a).norm();
        a = b;
    }
    return len;
}


Chain::Chain()
{
    fnCut          = 0;
    fnSegmentation = 0;
    fnAbscissaM = 0;
    fnAbscissaP = 0;
#if FIBER_HAS_NORMAL
    fnNormal.set(0, 0, 1);
#endif
    cDeltaM = 0;
    cDeltaP = 0;
    needUpdate  = false;
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
    Vector delta = dir * ( fnCut / dir.norm() );
    //
    for ( index_t p = 0 ; p < nPoints; ++p )
        setPoint(p, pos + p * delta);
}


void Chain::setStraight(Vector const& pos, Vector const& dir, real len)
{
    assert_true( fnSegmentation > REAL_EPSILON );
    assert_true( len > 0 );
    index_t np = bestNumberOfPoints(len/fnSegmentation);
    setNbPoints(np);
    setSegmentation(len/(np-1));
    fnAbscissaP = fnAbscissaM + len;
    setStraight(pos, dir);
}



void Chain::setCurved(Vector dir, real rad, real len, real off)
{
    assert_true( fnSegmentation > REAL_EPSILON );
    assert_true( len > 0 );
    // ensure that `dir` is orthogonal to X:
    dir.XX = 0;
    dir.normalize();
    index_t np = bestNumberOfPoints(len/fnSegmentation);
    setNbPoints(np);
    real S = len / ( np - 1 );
    setSegmentation(S);
    fnAbscissaP = fnAbscissaM + len;
    real delta = 2 * std::asin(0.5*S/rad);
    real angle = -0.5 * delta * np;
    for ( index_t p = 0 ; p < np; ++p, angle += delta )
    {
        real S = rad * std::sin(angle);
        real C = rad * std::cos(angle) - off;
        setPoint(p, Vector(S, 0, 0) + dir*C);
    }
    updateFiber();
}


/**
 The filament is set as a random walk with given persistence length

 This return a filament with the center of gravity at zero
 and the average orientation aligned with (1, 0, 0)
 */
void Chain::setEquilibrated(real len, real persistence_length)
{
    index_t np = bestNumberOfPoints(len/fnSegmentation);
    assert_true( np > 1 );
    
    setNbPoints(np);
    setSegmentation(len/(np-1));
    fnAbscissaP = fnAbscissaM + len;
    
    real sigma = std::sqrt(2*fnCut/persistence_length);
    
    Vector pos(0,0,0);
    Vector dir(1,0,0);
    setPoint(0, pos);
    
    for ( index_t p = 1 ; p < np; ++p )
    {
        pos += fnCut * dir;
        setPoint(p, pos);
        //rotate dir in a random direction:
        real a = sigma * RNG.gauss();
        dir = std::cos(a) * dir + dir.randOrthoU(std::sin(a));
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
            translate(posMiddle()-posEndP());
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
 The MINUS_END and plus end will be set to the first and last points in `pts[]`,
 and intermediate points will be interpolated at regular intervals on `pts[]`.
 
 The length of the resulting fiber will match the sum of given segment lengths,
 and the segments will be approximately equal to each other.
 Thus reshape() is called eventually to equalize the segments.
 */
void Chain::setShape(const real pts[], size_t n_pts, index_t np)
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
    index_t p = 1;
    --np;
    
    for ( index_t n = 1; n < np; ++n )
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
    reshape();
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
    LOG_ONCE("WARNING: Cytosim was compiled without FIBER_HAS_NORMAL\n");
    return d.orthogonal();
#endif
}


//==============================================================================
#pragma mark - Reshape

/*
 This deals with Fiber having one segment only,
 for which the procedure is straightforward
 */
void Chain::reshape_two(const real* src, real* dst, real cut)
{
    real X = src[  DIM] - src[0];
#if ( DIM == 1 )
    real s = 0.5 - 0.5 * (cut/abs_real(X));
#elif ( DIM == 2 )
    real Y = src[1+DIM] - src[1];
    real n = std::sqrt( X * X + Y * Y );
    real s = 0.5 - 0.5 * (cut/n);
#else
    real Y = src[1+DIM] - src[1];
    real Z = src[2+DIM] - src[2];
    real n = std::sqrt( X * X + Y * Y + Z * Z );
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
 
 i.e. 'A[]' should fulfill a set of equations F[i] = 0, with:
 
     F[i] = ( Y[i+1] - Y[i] )^2 - cut^2
 
 Note that:
 
     Y[i+1] - Y[i] = A[i+1] * dif[i+1] + (1-2*A[i]) * dif[i] + A[i-1] * dif[i-1]
 
 Method: use all zeros as first guess for 'sca', and follow
 Newton's method to iteratively refine the multidimensional guess.
 
 In practice, we calculate `A_new` from `A` using the relationship:
 
     J(A) * ( A_new - A ) = -F(A)
 
 Where J is the Jacobian matrix: J[i,j] = dF[i] / dA[j]
 
 For this problem, J is square and tri-diagonal but not symmetric,
 and must be recalculated at each iteration. A factor 2 can be factorized:

     A_new = A - 1/2 inv(K).F(A)
 
 Where J = 2 * K

 externally provided memory `mem[]` should be allocated to hold `5*chk` reals
 FJN, Strasbourg, 22.02.2015 & Cambridge, 10.05.2019 -- 13.05.2019
 */

int Chain::reshape_calculate(const index_t ns, real target, index_t max_iter,
                             real const* mag, real const* pri, real const* sec,
                             real* mem, size_t chk)
{
    real * mul = mem;
    real * rhs = mem+chk;
    real * dia = mem+chk*2;
    real * low = mem+chk*3;
    real * upe = mem+chk*4;
    constexpr real two = -2;
    real err, err0 = 0;
    const real tol = DIM * ns * REAL_EPSILON;
    /*
     Perform here the first iteration of Newton's method:
     the formula is the same as below, with `mul == 0`, and thus 'vec == dif'
     The system is symmetric definite positive, and we can use a faster factorization
     The matrix is multiplied by -2, to yield a -0.5 in the update: mul[] -= 0.5 * rhs[]
     */
    for ( index_t i = 0; i < ns; ++i )
    {
        mul[i] = mag[i] - target;
        dia[i] = mag[i] * 4;
        low[i] = pri[i] * two;  //accessing pri[ns-1], but not used
        err0 += abs_real(mul[i]);  // calculating the 1-norm
    }
    
    int info = 0;
    lapack::xpttrf(ns, dia, low, &info);
    if ( info ) {
        fprintf(stderr, " reshape_local (dpttrf) failed %i\n", info);
        return 1;
    }
    lapack::xptts2(ns, 1, dia, low, mul, ns);
    if ( err0 < tol )
        return 0;
    index_t cnt = 0;
    //fprintf(stderr, "\n %3lu err %20.16f norm(rhs) %8.5e\n", cnt, err0, blas::nrm2(ns, mul));

    while ( ++cnt <= max_iter )
    {
        assert_true( ns > 1 );
        // set the matrix elements and RHS of system,
        {
            real b = mul[0] * 2 - 1, c = mul[1];
            //vec = c*dif[i+1] - b*dif[i];
            real D = c * pri[0] - b * mag[0];
            real U = c * mag[1] - b * pri[0];
            rhs[0] = c * U - b * D - target;
            dia[0] = D * two;
            upe[0] = U;
        }
        err = abs_real(rhs[0]);
        for( index_t i = 1; i+1 < ns; ++i )
        {
            real a = mul[i-1];
            real b = mul[i] * 2 - 1;
            real c = mul[i+1];
            /*
             vec = a*dif[i-1] - b*dif[i] + c*dif[i+1];
             rhs = vec.normSqr() - target;
             low = derivate(rhs, mul[i-1])
             dia = derivate(rhs, mul[i]);
             upe = derivate(rhs, mul[i+1);
             */
            real L = a * mag[i-1] - b * pri[i-1] + c * sec[i-1];
            real D = a * pri[i-1] - b * mag[i  ] + c * pri[i  ];
            real U = a * sec[i-1] - b * pri[i  ] + c * mag[i+1];

            rhs[i] = a * L - b * D + c * U - target;
            low[i] = L;
            dia[i] = D * two;
            upe[i] = U;
            err += abs_real(rhs[i]);
        }
        {
            real a = mul[ns-2];
            real b = mul[ns-1] * 2 - 1;
            //vec = a * dif[ns-2] - b * dif[ns-1];
            real L = a * mag[ns-2] - b * pri[ns-2];
            real D = a * pri[ns-2] - b * mag[ns-1];
            rhs[ns-1] = a * L - b * D - target;
            low[ns-1] = L;
            dia[ns-1] = D * two;
            err += abs_real(rhs[ns-1]);
        }
        //fprintf(stderr, " %3lu err %20.16f norm(rhs) %8.5e", cnt, err, blas::nrm2(ns, rhs));
        //VecPrint::head("> rhs", ns, rhs, 3);

        if ( err < tol )
            return 0; // all good!
        if ( err > err0 )
        {
            // the solution is worse than before...
            //VecPrint::head("> mul", ns, mul, 3);
            //fprintf(stderr, "\n         %20.16f", tol);
            if ( err > 128 * tol ) return 1;
        }
        err0 = err;
#if ( 0 )
        printf("\n diff(L,U) = %f", blas::difference(ns-1, upe, low+1));
        printf("\n L"); VecPrint::head(ns-1, low+1, 3);
        printf("\n D"); VecPrint::head(ns  , dia, 3);
        printf("\n_U"); VecPrint::head(ns-1, upe, 3);
        //printf("\n rhs "); VecPrint::head(ns, rhs, 6);
#endif
#if ( 0 )
        real asy = 0, sup = 0;
        for ( index_t i = 0; i < ns-1; ++i )
        {
            sup = std::max(sup, abs_real(low[i+1]+upe[i]));
            asy += abs_real(low[i+1]-upe[i]);
        }
        printf("\n %3i diff(low-upe) %12.6f", cnt, 2*asy/sup);
#endif
        lapack::xgtsv(ns, 1, low+1, dia, upe, rhs, ns, &info);
        if ( info )
        {
            fprintf(stderr, " reshape_local (dgtsv) failed %i\n", info);
            return 2;
        }
        
        // update result following Newton's iteration
        //blas::xaxpy(ns, -0.5, rhs, 1, mul, 1);
        for ( index_t u = 0; u < ns; ++u )
            mul[u] -= 0.5 * rhs[u];
    }
    //printf("\n   >> err %20.16f", err);
    //printf("\n   >>mul "); VecPrint::print(ns, mul, 3);
    return ( err > err0 );
}


/**
 Apply correction of magnitude 'mul' along the segment directions:
 
 except for the edges, this is:
 P[i] <- P[i] + mul[i] * ( P[i+1] - P[i] ) - mul[i-1] * ( P[i] - P[i-1] )
 
 Source `src[]` and destination `dst[]` should be arrays of size `DIM*(nbs+1)`
 Array `mul[]` should be of size `nbs` containing the scalar multipliers
  
 The code is similar to projectForcesD(), with an additional difference
 
 Note that if 'src == dst', the result are still correct
 */
void Chain::reshape_apply(const index_t nbs, const real* src,
                          const real * mul, real* dst)
{
    assert_true( nbs > 1 );
    Vector A(src);
    Vector B = A;

    for ( index_t i = 0; i < nbs; ++i )
    {
        Vector C(src+DIM*(i+1));
        Vector add = mul[i] * ( C - B );
        (A+add).store(dst+DIM*i);
        A = C - add;
        B = C;
    }
    
    A.store(dst+DIM*nbs);
}


/**
 Try to re-establish the length of the segments, by moving points along
 the directions of the flanking segments
 
 here ns = nbSegments() and mem_size >= ns
 `mem` should be allocated to hold 8 * mem_size
 */
int Chain::reshape_local(const index_t nbs, const real* src, real* dst,
                         real cut, index_t max_iter, real* mem, size_t mem_size)
{
    assert_true( nbs > 1 );

    // alias some of the memory provided for the work:
    real * mag = mem + mem_size * 5;
    real * pri = mem + mem_size * 6;
    real * sec = mem + mem_size * 7;
    
    // rescale factor to make matrix elements close to 1:
    const real alpha = 1.0 / ( cut * cut );
    // calculate terms using 'dif[]':
    Vector A, B, C;
    A.load_diff(src);
    B.load_diff(src+DIM);
    mag[0] = alpha * A.normSqr();
    mag[1] = alpha * B.normSqr();
    pri[0] = alpha * dot(A, B);
    for ( index_t i = 2; i < nbs; ++i )
    {
        C.load_diff(src+DIM*i);
        mag[i] = alpha * C.normSqr();
        pri[i-1] = alpha * dot(B, C);
        sec[i-2] = alpha * dot(A, C);
        A = B;
        B = C;
    }
    // these terms should not be used:
    pri[nbs-1] = 0;
    sec[nbs-2] = 0;
    sec[nbs-1] = 0;
#if ( 0 )
    VecPrint::print("mag", nbs, mag, 3);
    VecPrint::print("pri", nbs, pri, 3);
    VecPrint::print("sec", nbs, sec, 3);
#endif
    int res = reshape_calculate(nbs, 1.0, max_iter, mag, pri, sec, mem, mem_size);

    if ( res == 0 )
    {
        //VecPrint::edges(">>> mul", nbs, mem, 3);
        reshape_apply(nbs, src, mem, dst);
    }
    return res;
}


/**
 Move the vertices relative to each other, such that when this is done,
 all segments have the same distance `cut` (ie. parameter `segmentation`).
 This is operation does not change the center of gravity of the fiber.

 NOTE: if two consecutive points overlap, there is no unique way to
 restore the constraints! We do nothing in that case, because most
 likely, the Brownian motion will push the points apart soon.

 Applying this method after a sudden perpendicular force is not ideal:
 For example, a force applied to the bottom of a vertical fibers leads
 to a 'L' configuration after one step of `Meca::solve`.
 reshape_global() reduces the bottom leg of the 'L', by translating the entire
 vertical portion of the fiber, irrespective of the length of this section.
 */

void Chain::reshape_global(const index_t ns, const real* src, real* dst, const real cut)
{
    Vector inc(0,0,0), sum(0,0,0), seg;
    seg.load_diff(src);
    real dis = seg.norm();
    
    // translation needed to restore first segment
    if ( dis > REAL_EPSILON )
        inc = ( cut/dis - 1.0 ) * seg;
    
    Vector(src).store(dst);

    for ( index_t i = 1; i < ns; ++i )
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
    
    // move the last point by `inc`:
    (inc+Vector(src+DIM*ns)).store(dst+DIM*ns);
    
    // calculate uniform motion needed to conserve the center of gravity:
    sum = ( sum + inc ) / ( ns + 1 );
    
    // translate entire fiber uniformly:
    for ( index_t i = 0; i <= ns; ++i )
        sum.sub_to(dst+DIM*i);
}


/**
 Replace coordinates by the ones provided in `ptr`, and call 'reshape()':
 the points are displaced to bring adjacent points to a distance 'segmentation'
 This is done in a way that can be parallellize on multiple threads
 */
void Chain::getPoints(real const* ptr)
{
    constexpr index_t NVEC = 8U;
#if 0
    // using static memory
    // Attention: this only works if using a single thread
    static size_t alc = 0;
    static real* mem = nullptr;
    
    if ( alc < allocated() )
    {
        alc = allocated();
        free_real(mem);
        mem = new_real(alc*NVEC);
    }
#else
    // using per-thread static memory
    thread_local static size_t alc = 0;
    
    auto delete_real = [](real * x)
    {
        //printf("> del real[%lu] %p %p\n", alc, x, pthread_self());
        free_real(x);
        alc = 0;
    };

    thread_local static std::unique_ptr<real, decltype(delete_real)> uptr(nullptr, delete_real);
    
    if ( alc < allocated() || !uptr.get() )
    {
        alc = allocated();
        free_real(uptr.release());
        uptr.reset(new_real(alc*NVEC));
        //printf("> new real[%lu] %p %p\n", alc, uptr.get(), pthread_self());
    }
    real * mem = uptr.get();
#endif

#if ( DIM > 1 )
    if ( nPoints == 2 )
        reshape_two(ptr, pPos, fnCut);
    else if ( reshape_local(nbSegments(), ptr, pPos, fnCut, 16, mem, allocated()) )
#endif
    {
        //VecPrint::print("/ ", DIM*nPoints, pPos, 2);
        //VecPrint::print("L ", DIM*nPoints, ptr, 2);
        std::string doc = document(ptr);
        real mov = std::sqrt(sum_square(DIM*nPoints, pPos, ptr)) / nPoints;
        reshape_global(nbSegments(), ptr, pPos, fnCut);
#if ( DIM > 1 )
        Cytosim::log(" wild motion for ", doc, " displacement ", mov, '\n');
        //copy_real(DIM*nPoints, ptr, pPos);
#endif
    }
    //document(std::clog, ptr);
    //fprintf(stderr, "\n%3u >pos ", identity()); VecPrint::print(DIM*nPoints, pPos, 6, DIM);
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
    const index_t ns = nbSegments();
    
    if ( delta > 0 )
    {
        index_t p = 0, n = ns;
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
        for ( index_t p = 0, n = ns; n > 0; ++p, --n )
            movePoint(p, ( a * n ) * diffPoints(p));
    }
    
    cDeltaM = delta;
    fnAbscissaM -= delta;
    setSegmentation(length()/ns);
    postUpdate();
}

/**
 This extends the fiber by adding one segment at the minus end.
 Thus `segmentation()` is not changed, and the existing points are not displaced.
 */
void Chain::addSegmentM()
{
    index_t pp = 1+nPoints;
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
 A portion of size `delta` is removed at the minus end, and vertices are recalculated.
 The Fiber length is reduced by `delta` ( which must be >= 0 ).
 
 Note: after cutM(), the distance between adjacent vertices is not exactly
 equal to segmentation(). This is true only if the fiber is straight.
 */
void Chain::cutM(const real delta)
{
    real len = length();
    const index_t ns = bestNumberOfPoints((len-delta)/fnSegmentation) - 1;
    const real cut = (len-delta) / ns;
    assert_true( cut > 0 );
    real* tmp = new_real(DIM*(ns+1));

    // calculate intermediate points:
    for ( index_t i=0; i < ns; ++i )
    {
        Vector w = interpolateM(delta+i*cut).pos();
        w.store(tmp+DIM*i);
    }

    // copy the position of plus-end:
    copy_real(DIM, pPos+DIM*lastPoint(), tmp+DIM*ns);
    
    setNbPoints(ns+1);
    fnAbscissaM += delta;
    setSegmentation(cut);
    getPoints(tmp);
    free_real(tmp);
    postUpdate();
}


void Chain::truncateM(index_t p)
{
    Mecable::truncateM(p);
    fnAbscissaM = abscissaPoint(p);
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
    const index_t ns = nbSegments();
    
    if ( delta > 0 )
    {
        Vector dp0 = diffPoints(ns-1), dp1;
        movePoint(ns, ( a * ns ) * dp0);
        index_t p = ns-1;
        
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
        //printf("+ %12.6f\n", dirEndP().norm());
    }
    else if ( delta < 0 )
    {
        for ( index_t p = ns ; p > 0 ; --p )
            movePoint(p, ( a * p ) * diffPoints(p-1));
#if 0
        for ( index_t n = 0; n < ns; ++n )
            printf(" %8.5f", diffPoints(n).norm()/len);
        printf(" - \n");
#endif
    }
    cDeltaP = delta;
    fnAbscissaP += delta;
    setSegmentation(length()/ns);
    postUpdate();
}


/**
 This extends the fiber by adding one segment at the plus end.
 Thus `segmentation()` is not changed, and the existing points are not displaced.
 */
void Chain::addSegmentP()
{
    index_t pp = nPoints;
    setNbPoints(pp+1);
    
    real * psp = pPos + DIM * pp;
    for ( index_t dd = 0; dd < DIM; ++dd )
        psp[dd] = 2 * psp[dd-DIM] - psp[dd-2*DIM];
    
    fnAbscissaP += fnCut;
    postUpdate();
}


/**
 A portion of size `delta` is removed at the plus end, and vertices are recalculated.
 The Fiber length is reduced by `delta` ( which must be >= 0 ).

 Note: after cutP(), the distance between adjacent vertices is not exactly
 equal to segmentation(). This is true only if the fiber is straight.
*/
void Chain::cutP(const real delta)
{
    real len = length();
    const index_t np = bestNumberOfPoints((len-delta)/fnSegmentation);
    const real cut = (len-delta) / (np-1);
    assert_true( cut > 0 );
    real* tmp = new_real(DIM*np);

    // copy minus end:
    copy_real(DIM, pPos, tmp);

    // calculate intermediate points:
    for ( index_t i = 1; i < np; ++i )
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


void Chain::truncateP(index_t p)
{
    Mecable::truncateP(p);
    fnAbscissaP = abscissaPoint(p);
    postUpdate();
}

//------------------------------------------------------------------------------
#pragma mark -

real Chain::freshAssembly(const FiberEnd end) const
{
    if ( end == PLUS_END )
        return freshAssemblyP();
    if ( end == MINUS_END )
        return freshAssemblyM();
    ABORT_NOW("invalid argument value");
    return 0;
}


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


/**
 `fib` is added at the plus end of `*this`
 
 The vertex are reinterpolated linearly, and the length of the
 segments will not fullfil exactly the constraints of segmentation.
 If this is a problem, Chain::reshape() should be called.
 
 `fib` should usually be destroyed afterward.
 */
void Chain::join(Chain const* fib)
{
    const real len = length();
    const real lenT = len + fib->length();
    const index_t ns = bestNumberOfPoints(lenT/fnSegmentation) - 1;
    const real cut = lenT / real(ns);
    
    real* tmp = new_real(DIM*(ns+1));

    // calculate new points into tmp[]:
    for ( index_t i = 1; i < ns; ++i )
    {
        Vector w;
        if ( i*cut < len )
            w = interpolateM(i*cut).pos();
        else
            w = fib->interpolateM(i*cut-len).pos();
        w.store(tmp+DIM*i);
    }
    
    // copy position of plus end:
    fib->posEndP().store(tmp+DIM*ns);

    setNbPoints(ns+1);
    setSegmentation(cut);
    fnAbscissaP = fnAbscissaM + cut * ns;
    getPoints(tmp);
    free_real(tmp);
    updateFiber();
}


/**
 Flip all the points, such that minus_end becomes plus_end and vice-versa.
 This does not affects Abscissa and the abscissa of center thus stays as it is.
*/
void Chain::flipChainPolarity()
{
    index_t ii = 0;
    index_t jj = lastPoint();
    
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

//------------------------------------------------------------------------------
#pragma mark -

static real lengthSegmentSqr(const real a[])
{
    real x = square( a[DIM] - a[0] );
    for( int i = 1; i < DIM; ++i )
        x += square( a[DIM+i] - a[i] );
    return x;
}

/**
 Returns the minimum and maximum distance between consecutive points, for 'cnt' segments
 */
void Chain::computeMinMax(index_t cnt, real const* ptr, real& in, real& ax)
{
    real x = lengthSegmentSqr(ptr);
    in = x;
    ax = x;
    for ( index_t i = 1; i < cnt; ++i )
    {
        x = lengthSegmentSqr(ptr+DIM*i);
        in = std::min(in, x);
        ax = std::max(ax, x);
    }
    in = std::sqrt(in);
    ax = std::sqrt(ax);
}

/**
 Returns the average and variances of segment length
 */
void Chain::computeMeanVar(index_t cnt, real const* ptr, real off, real& mean, real& variance)
{
    Vector vec;
    double avg = 0, var = 0;
    for ( index_t n = 0; n < cnt; ++n )
    {
        vec.load_diff(ptr+DIM*n);
        real r = vec.norm() - off;
        avg += r;
        var += r*r;
    }
    avg /= cnt;
    variance = var - square(avg) * cnt;
    if ( cnt > 1 )
        variance /= real(cnt-1);
    mean = avg + off;
}


/**
 Calculate Menger curvature:
 the inverse of the radius of the circle that passes through A, B and C
 
 1/R = 4 * Area(triangle) / ( |AB|*|BC|*|AC| )
 
 curvature is ZERO if A, B and C are aligned
 */
[[ maybe_unused ]]
static real curvature3(Vector const& A, Vector const& B, Vector const& C)
{
    // cross-product = 2 * surface; hence S = 16 * surface^2
    real S = 4 * normSqr(cross( C - A, C - B ));
    return sqrt( S / ( normSqr(B-A) * normSqr(C-B) * normSqr(C-A) ));
}

/**
 Calculate Menger curvature:
 the inverse of the radius of the circle that passes through A, B and C
 
 1/R = 4 * Area(triangle) / ( |AB|*|BC|*|AC| )
 
 curvature is ZERO if A, B and C are aligned
 
 Thank you, Serge to point this out!
 */
[[ maybe_unused ]]
static real curvature3(Vector const& A, Vector const& B, Vector const& C, real seg)
{
    real H = norm( A + C - 2 * B );  // 2 * height_of_triangle
    return H / ( seg * seg );
}

/**
 This returns the curvature estimated at point `p`
 */
real Chain::curvature(index_t p) const
{
    if (( 0 < p ) & ( p < lastPoint() ))
        return curvature3(posP(p-1), posP(p), posP(p+1), fnCut);
    return 0;
}


/**
 The normalized bending energy is an integral over the curvilinear abscissa `s`:
 
     1/2 * sum( curvature(s)^2 ds )
 
 The curvature is calculated from the positions of the vertices:
 Given theta = angle between two consecutive segments,

     curvature = 1/R = 2 * std::sin(angle/2) / segmentation
 
 and since
 
     sin^2(angle/2) = ( 1 - std::cos(angle) ) / 2
 
 hence:
 
     1/2 * curvature^2 = ( 1 - std::cos(angle) ) / ( segmentation^2 )
 
 and finaly:
 
     1/2 * sum( curvature^2 * ds ) = sum( 1 - std::cos(angle) ) / segmentation
 
 */
real Chain::bendingEnergy0() const
{
    real e = 0;
    const index_t lsp = nPoints - 2;
    if ( lsp > 0 )
    {
        for ( index_t p = 0; p < lsp ; ++p )
        {
            Vector A = posP(p);
            Vector B = posP(p+1);
            Vector C = posP(p+2);
            e += dot(B - A, C - B);  // e += std::cos(angle) * segmentation^2
        }
        // e <- sum( 1 - std::cos(angle) )
        e = lsp - e / ( fnCut * fnCut );
        
        /*
         We correct the result, because we only considered (nPoints-2) junctions,
         and thus covered only a fraction of the total length of the filament
         */
        e *= ( lsp + 1 ) / ( fnCut * lsp );
    }
    return e;
}


real Chain::minCosine() const
{
    real res;
    Vector dir1, dir2;
    
    index_t ps = nbSegments() % 2;
    if ( ps )
    {
        dir1 = diffPoints(0);
        res = fnCut * fnCut;
    }
    else
    {
        dir1 = diffPoints(1);
        res = dot(diffPoints(0), dir1);
        ps = 2;
    }
    
    for ( ; ps < nbSegments(); ps += 2 )
    {
        dir2 = diffPoints(ps);
        res = std::min(res, dot(dir1, dir2));
        dir1 = diffPoints(ps+1);
        res = std::min(res, dot(dir1, dir2));
    }
    
    return res / ( fnCut * fnCut );
}


/**
 Returns the minimum and maximum distance between consecutive points
 */
index_t Chain::nbKinks(real threshold) const
{
    index_t res = 0;
    threshold *= fnCut * fnCut;
    Vector d = diffPoints(0);
    
    for ( index_t n = 1; n < lastPoint(); ++n )
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

real Chain::planarIntersect(index_t s, Vector const& n, const real a) const
{
    assert_true( s < nbSegments() );
    
    real sca = dot(diffPoints(s), n);
    
    // if segment is parallel to plane, there is no intersection:
    if ( abs_real(sca) < REAL_EPSILON )
        return INFINITY;
    
    Vector pos = posP(s);
    
    if ( modulo )
    {
        // calculate image that is closest to plane:
        Vector cen = n * ( -a / n.normSqr() );
        modulo->fold(pos, cen);
    }
    
    return - ( dot(pos, n) + a ) / sca;
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 Recalculate the vertices to get 'ns' segments.
 
 @todo 2d-order interpolation in Chain::resegment()
 
 Note: Unless the Chain is straight, the length of the segments after this
 new interpolation will not exactly match `segmentation()`, but this calls
 getPoints() here which calls reshape().
 */
void Chain::resegment(index_t ns)
{
    //std::clog << reference() << " resegment(" << ns << ")\n";
    assert_true( ns > 0 );
    // 'cut' is in unit of segments
    real cut = (real)nbSegments() / (real)ns;
    
    // calculate new intermediate points in tmp[]:
    Vector a = posP(0), b = posP(1);
    
    real h = 0;
    index_t p = 1;
    real* tmp = new_real(DIM*(ns+2));
    
    a.store(tmp);
    for ( index_t n = 1; n < ns; ++n )
    {
        h += cut;
        
        while ( h > 1.0 )
        {
            h -= 1.0;
            a = b;
            ++p;
            assert_true(p<nbPoints());
            b.load(pPos+DIM*p);
        }
        
        Vector w = a + h * ( b - a );
        w.store(tmp+DIM*n);
    }
    
    // copy coordinates of last point:
    a.load(pPos+DIM*lastPoint());
    a.store(tmp+DIM*ns);

    // resize filament:
    setNbPoints(ns+1);
    setSegmentation(length()/ns);
    getPoints(tmp);
    free_real(tmp);
}


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
    
    index_t best = bestNumberOfPoints(nbSegments()*fnCut/fnSegmentation);
    
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


void Chain::adjustSegmentation(real arg)
{
    if ( (float)fnSegmentation != (float)arg )
    {
        std::clog << reference() << " segmentation " << fnSegmentation << " -> " << arg << "\n";
        fnSegmentation = arg;
        adjustSegmentation();
        updateFiber();
    }
}


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
    {
        std::string str = "["+std::to_string(abscissaM())+" "+std::to_string(abscissaP())+"]";
        throw InvalidParameter("hand::abscissa is out of range "+str);
    }
    return dis;
}

/**
 The Fiber is partitionned by this function in three regions:
 - a MINUS_END part of length `lambda`
 - a PLUS_END part, also of length `lambda`
 - and a NO_END section in between
 .
 Note that a Fiber shorter than `2*lambda` does not have a central region,
 and is composed of plus end and minus end parts of equal size.
 */    
FiberEnd Chain::whichEndDomain(const real ab, const real lambda) const
{
    const real abM = ab - fnAbscissaM;
    const real abP = fnAbscissaP - ab;
    
    if ( abM > abP )
    {
        if ( abP <= lambda )
            return PLUS_END;
    }
    else
    {
        if ( abM <= lambda )
            return MINUS_END;
    }
    return NO_END;
}


//------------------------------------------------------------------------------
#pragma mark -

Interpolation Chain::interpolateEndM() const
{
    return Interpolation(this, real(0.0), 0);
}


Interpolation Chain::interpolateEndP() const
{
    return Interpolation(this, real(1.0), nPoints-2);
}


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
            return Interpolation(this, real(0.0), 0);
        case PLUS_END:
            return Interpolation(this, real(1.0), nPoints-2);
        case CENTER:
            return interpolateCenter();
        default:
            throw Exception("unexpected argument");
    }
}


Interpolation Chain::interpolateCenter() const
{
    index_t n = lastPoint() / 2;
    return Interpolation(this, real(0.5)*(lastPoint()-2*n), n);
}


/**
 return Interpolation corresponding to a distance `abs` from the minus end
 The interpolation describes a position:
       X = P(i) * (1-a) + P(i+1) * a
 where
 - `i` is an integer: 0 <= i < lastPoint(),
 - `a` is a positive real coefficient: 0 <= a <= 1
 .
 When `abs` is above the plus end, an interpolation of the last point is returned.
 
 */
Interpolation Chain::interpolateM(const real ab) const
{
    real a = max_real(ab*iCut, 0);
    index_t i(nbPoints()-2);
    i = std::min(i, (index_t)a);
    return Interpolation(this, std::min(a-i, real(1)), i);
}


Interpolation Chain::interpolateAbs(const real ab) const
{
    real a = max_real((ab-fnAbscissaM)*iCut, 0);
    unsigned i = std::min((unsigned)a, (unsigned)nPoints-2);
    return Interpolation(this, std::min(a-i, real(1)), i);
}


Interpolation Chain::interpolateFrom(const real ab, const FiberEnd end) const
{
    switch( end )
    {
        case ORIGIN:
            return interpolateM(ab-fnAbscissaM);
            
        case MINUS_END:
            return interpolateM(ab);
            
        case CENTER:
            return interpolateM(ab + 0.5*length());
            
        case PLUS_END:  //this is counted from the plus towards the minus end
            return interpolateM(fnCut*nbSegments() - ab);
        
        default:
            ABORT_NOW("invalid argument value");
    }
    return interpolateM(-fnAbscissaM);
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
    index_t s = static_cast<index_t>(a);
    
    // check if plus end is reached:
    if ( s+1 < nPoints )
        return midPoint(s, a-s);
    else
        return posP(lastPoint());
}

Vector Chain::dirM(const real ab) const
{
    // at MINUS_END
    if ( ab <= 0 )
        return dirSegment(0);
    
    index_t a = static_cast<index_t>( ab / fnCut );
    index_t t(nbPoints()-2);
    return dirSegment(std::min(t, a));
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


/// force on the minus end projected on the direction of elongation
real Chain::projectedForceEndM() const
{
    return -dot(netForce(0), dirSegment(0));
}

/// force on the plus end projected on the direction of elongation
real Chain::projectedForceEndP() const
{
    return dot(netForce(lastPoint()), dirSegment(lastSegment()));
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
#pragma mark - document

void Chain::briefdoc(std::ostream& os, real len, real con, real mn, real mx) const
{
    std::streamsize p = os.precision(3);
    os << reference();
    os << " seg " << segmentation() << ": " << mn << " +" << mx-mn;
    os << " len " << len << " " << std::showpos << con-len << std::noshowpos << " ";
    os.precision(p);
}


/**
 Prints info on the length of Segments, which can be useful for debugging
 */
void Chain::document(std::ostream& os, real len, real con, real mn, real mx) const
{
    os << "chain " << std::setw(7) << reference() << '\n';
    os << "{\n";
    os << "    segmentation = " << segmentation() << '\n';
    os << "    segment_min = " << mn << '\n';
    os << "    segment_max = " << mx << '\n';
    os << "    length  = " << len << '\n';
    os << "    contour = " << con << '\n';
    os << "}" << std::endl;
}


int Chain::checkLength(real const* ptr, std::ostream& os, real len) const
{
    real mn, mx;
    computeMinMax(nbSegments(), ptr, mn, mx);
    real dev = ( mx - mn ) / segmentation();
    real con = contourLength(ptr, nPoints);
    real err = abs_real( con - len ) / ( con + len );
    int res = ( dev > 0.02 ) + ( err > 0.02 );
    if ( res )
        document(os, len, con, mn, mx);
    return res;
}


/**
 Prints info on the length of Segments, which can be useful for debugging
 */
void Chain::document(std::ostream& os, real const* ptr) const
{
    real mn, mx;
    real len = length();
    computeMinMax(nbSegments(), ptr, mn, mx);
    real con = contourLength(ptr, nPoints);
    briefdoc(os, len, con, mn, mx);
}


std::string Chain::document(real const* ptr) const
{
    std::ostringstream ss;
    document(ss, ptr);
    return ss.str();
}


//------------------------------------------------------------------------------
#pragma mark - I/O

void Chain::write(Outputter& out) const
{
    assert_small( length1() - length() );
    out.writeUInt32(signature());
    out.writeFloat(length());
    out.writeFloat(fnSegmentation);
    out.writeFloat(fnAbscissaM);
    Mecable::write(out);
}


/**
 The fiber will be re-segmented if its current desired segmentation 
 does not match the one stored in the file.
 This saves the extension since last frame in cDeltaM, cDeltaP
 */
void Chain::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    //Cytosim::log(" Chain::read at ", in.pos(), '\n');
    ObjectSignature s = in.readUInt32();
    if ( s ) signature(s);
    
    float len = in.readFloat();
    float seg = in.readFloat(); // saved target value for segmentation
    float abs = in.readFloat();
    
#if BACKWARD_COMPATIBILITY < 60
    if ( in.formatID() > 49 && in.formatID() < 59 ) // birthTime [ 12.12.2018, 11.01.2023 ]
    {
        // we are discarding the birthTime information... for old formats only
        in.readFloat();
    }
#endif

    Mecable::read(in, sim, tag);
    
    if ( len <= 0 )
        throw InvalidIO("invalid (negative) fiber length");

    if ( len > 1e6 )
        throw InvalidIO("excessive fiber length");
    
    if ( seg <= 1e-9 || seg > 1e6 )
        throw InvalidIO("invalid fiber segmentation");

    if ( nPoints < 2 )
        throw InvalidIO("invalid fiber with 0 or 1 point");

    bool valid = ( fnAbscissaP > fnAbscissaM );
    if ( valid )
    {
        cDeltaM = +fnAbscissaM;
        cDeltaP = -fnAbscissaP;
    }
    
    fnAbscissaM = abs;
    fnAbscissaP = abs + len;
    
    if ( valid )
    {
        cDeltaM -= fnAbscissaM;
        cDeltaP += fnAbscissaP;
    }
    
    setSegmentation(len/nbSegments());  // set segments' length
    //Mecable::write(std::cerr);
    
    if ( seg != targetSegmentation() )
    {
        //std::clog << reference() << " segmentation " << seg << " -> " << targetSegmentation() << "\n";
        //targetSegmentation(seg); // keep value saved in file
        adjustSegmentation(); // resegment the fiber
    }
    
#ifndef NDEBUG
    // verify the length and segmentation:
    if ( in.vectorSize() == DIM )
        checkLength(pPos, std::clog, len);
#endif
}


/**
 Write angles between consecutive points instead of coordinates.
 This reduces size, without loosing information if the segments all have
 the same length, which should normally be the case.
 The angles can be stored on 2 bytes, as their range is known: [-PI, PI]
 There is a loss of precision, but overflow is not possible
 */
void Chain::writeAngles(Outputter& out) const
{
    //out.writeUInt32(signature());
    out.writeFloatBinary(length());
    //out.writeFloat(fnSegmentation);
    out.writeFloatBinary(fnAbscissaM);
    out.writeUInt16Binary(nPoints-1);
    // first point:
    out.writeFloatsBinary(pPos, DIM);
    // angles:
    for ( index_t i = DIM; i < DIM*nPoints; i += DIM )
    {
        // in 2D, there is only one angle:
        real x = pPos[i  ] - pPos[i-DIM  ];
#if ( DIM > 1 )
        real y = pPos[i+1] - pPos[i-DIM+1];
        real a = std::atan2(y, x);
#else
        real a = 0;
#endif
#if ( DIM >= 3 )
        real z = pPos[i+2] - pPos[i-DIM+2];
        // the second angle is always positive, in [0, PI]
        real b = std::atan2(std::sqrt(x*x+y*y), z);
        out.writeEulerAngles(a, b);
#else
        out.writeAngle(a);
#endif
    }
}


/**
 Read angles, and recalculate the coordinates of the points
 This might be slow due to trigonometry calls
 */
void Chain::readAngles(Inputter& in, Simul&, ObjectTag)
{
    //Cytosim::log(" Chain::readAngles at ", in.pos(), '\n');
    float len = in.readFloatBinary();

    if ( len <= 0 )
        throw InvalidIO("invalid (negative) fiber length");

    if ( len > 1e6 )
        throw InvalidIO("excessive fiber length");

    fnAbscissaM = in.readFloatBinary();
    fnAbscissaP = fnAbscissaM + len;
    index_t cnt = in.readUInt16();
    
    if ( cnt < 1 )
        throw InvalidIO("invalid fiber with no segment?");
    
    const float L = len / float(cnt);

    setNbPoints(cnt+1);
    // read first point:
    float pos[4] = { 0 };
    in.readFloats(pos, DIM);
    for ( int d = 0; d < DIM; ++d )
        pPos[d] = pos[d];
    if ( in.vectorSize() > 2 )
    {
        // 3D data with two angles per segment
        for ( index_t i = 1; i <= cnt; ++i )
        {
            float a, b;
            in.readEulerAngles(a, b);
            pos[0] += L * cosf(a) * sinf(b);
            pos[1] += L * sinf(a) * sinf(b);
            pos[2] += L * cosf(b);
            for ( index_t d = 0; d < DIM; ++d )
                pPos[DIM*i+d] = pos[d];
        }
    }
    else
    {
        // 2D case, one angle per segment
        for ( index_t i = 1; i <= cnt; ++i )
        {
            float a = in.readAngle();
            pos[0] += L * cosf(a);
            pos[1] += L * sinf(a);
            for ( index_t d = 0; d < DIM; ++d )
                pPos[DIM*i+d] = pos[d];
        }
    }
    setSegmentation(L);
    //checkLength(pPos, std::clog, len);
}


