// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "assert_macro.h"
#include "filament.h"
#include "iowrapper.h"
#include "messages.h"
#include "mecapoint.h"
#include "interpolation.h"
#include "fiber_site.h"
#include "exceptions.h"
#include "glossary.h"
#include "clapack.h"
#include "modulo.h"
#include "vector3.h"
#include "simul.h"

extern Modulo const* modulo;

/**
 This returns N+1, where N is the integer that minimizes
     fabs( length / N - segmentation ),
 */
unsigned Filament::bestNumberOfPoints(const real ratio)
{
    unsigned n = (int)ratio;
    
    if ( (2*n+1)*ratio >= 2*n*(n+1) )
        return n+2;
    
    return n+1;
}


real Filament::contourLength(const real* pts, unsigned n_pts)
{
    real len = 0;
    Vector a(pts), b;
    for ( unsigned n = 1; n < n_pts; ++n )
    {
        b.load(pts+DIM*n);
        len += (b-a).norm();
        a = b;
    }
    return len;
}


Filament::Filament()
{
    fnNormal.set(0, 0, 1);
    fnCut           = 0;
    fnSegmentation  = 0;
    fnAbscissaM     = 0;
    fnAbscissaP     = 0;
    fnBirthTime     = 0;
#if CURVATURE_DEPENDENT_SEGMENTATION
    fnCutError      = 0;
    // set different 'seeds' to desynchronize the adjustSegmentation()
    fnCutErrorIndex = RNG.pint(RECUT_PERIOD);
#endif
    needUpdate      = false;
}



//------------------------------------------------------------------------------
#pragma mark -

/**
This does not change length or segmentation
*/
void Filament::setStraight(Vector const& pos, Vector const& dir)
{
    assert_true( dir.norm() > 0.1 );
    // 'dir' is normalized for safety:
    Vector dpts = dir * ( fnCut / dir.norm() );
    //
    for ( unsigned p = 0 ; p < nPoints; ++p )
        setPoint( p, pos + p * dpts );
}


void Filament::setStraight(Vector const& pos, Vector const& dir, real len, const FiberEnd ref)
{
    assert_true( fnSegmentation > REAL_EPSILON );

    if ( len <= 0 )
        throw InvalidParameter("fiber:length must be > 0");

    unsigned np = bestNumberOfPoints(len/fnSegmentation);
    
    setNbPoints(np);
    setSegmentation(len/(np-1));

    switch( ref )
    {
        case MINUS_END:
            setStraight(pos, dir);
            break;
            
        case PLUS_END:
            setStraight(pos+dir*len, -dir);
            break;
            
        case CENTER:
            setStraight(pos-0.5*dir*len, dir);
            break;
            
        default:
            ABORT_NOW("invalid argument `ref`");
    }
    postUpdate();
}


/**
 This will set the Fiber with `np` points unless `np == 0`, in which case
 the number of points will be set automatically from fnSegmentation.
 pts[] should be of size DIM * n_pts and contain coordinates.

 The given set of points do not need to be equally distributed.
 The MINUS_END and PLUS_END will be set to the first and last points in `pts[]`,
 and intermediate points will be interpolated at regular intervals on `pts[]`.
 
 The length of the resulting fiber will be roughly equal to the sum of all segment lengths.
 However, the length of the segments will only be approximately equal to each other,
 and reshape() should be called to equalize them if necessary.
 */
void Filament::setShape(const real pts[], unsigned n_pts, unsigned np)
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
    
    a.load(pts);
    b.load(pts+DIM);
    setPoint(0, a);
    
    len = (b-a).norm();
    real h = 0;
    unsigned p = 1;
    --np;
    
    for ( unsigned n = 1; n < np; ++n )
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
    postUpdate();
}

/**
 The filament is set as a random walk with given persistence length

 This return a filament in a random direction, with the center of gravity at zero
 and the average orientation aligned with (1, 0, 0)
 */
void Filament::setEquilibrated(real len, real persistence_length)
{
    unsigned np = bestNumberOfPoints(len/fnSegmentation);
    assert_true( np > 1 );
    
    setNbPoints(np);
    setSegmentation(len/(np-1));
    
    real sigma = sqrt(2*fnCut/persistence_length);
    
    Vector dir(1,0,0);
    Vector vec(0,0,0);
    setPoint(0, vec);
    
    for ( unsigned p = 1 ; p < np; ++p )
    {
        vec += fnCut * dir;
        setPoint(p, vec);
        //rotate dir in a random direction:
        real a = sigma * RNG.gauss();
        real c = cos(a), s = sin(a);
        dir = c * dir + s * dir.randOrthoU(1);
    }
    
    // cancel out mean orientation and position:
    translate(-0.5*vec);
    if ( vec.normSqr() > 0.01 * fnCut )
    {
        Rotation rot = Rotation::rotationToVector(vec).transposed();
        rotate(rot);
    }
    postUpdate();
}


/**
 This adjusts the current `normal` or makes a new one if necessary
 (used for display)
 */
Vector3 Filament::adjustedNormal(Vector3 const& d) const
{
    if ( fnNormal.normSqr() < 0.8 || dot(fnNormal, d) > 0.5 )
        fnNormal = d.orthogonal(1.0);
    else
        fnNormal = d.orthogonal(fnNormal, 1.0);
    return fnNormal;
}


real Filament::age() const
{
    return simul().time() - fnBirthTime;
}


//===================================================================
#pragma mark -

/*
 This deals with Fiber having one segment only,
 for which the procedure is trivial
 */
void Filament::reshape_two(const real* src, real* dst, real cut)
{
#if ( DIM == 1 )
        Vector dif(src[1]-src[0]);
#elif ( DIM == 2 )
        Vector dif(src[2]-src[0], src[3]-src[1]);
#else
        Vector dif(src[3]-src[0], src[4]-src[1], src[5]-src[2]);
#endif
        dif *= 0.5 * ( 1 - cut/dif.norm() );
        
        dst[0    ] = src[0    ] + dif.XX;
        dst[  DIM] = src[  DIM] - dif.XX;
#if ( DIM > 1 )
        dst[1    ] = src[1    ] + dif.YY;
        dst[1+DIM] = src[1+DIM] - dif.YY;
#endif
#if ( DIM > 2 )
        dst[2    ] = src[2    ] + dif.ZZ;
        dst[2+DIM] = src[2+DIM] - dif.ZZ;
#endif
}


/**
 Apply correction
 */
void Filament::reshape_apply(const unsigned ns, const real* src, real* dst,
                             const real * sca, const Vector* dif)
{
    assert_true( ns > 1 );
    Vector d, e = sca[0] * dif[0];
    
    dst[0] = src[0] + e.XX;
#if ( DIM > 1 )
    dst[1] = src[1] + e.YY;
#endif
#if ( DIM > 2 )
    dst[2] = src[2] + e.ZZ;
#endif
    
    for ( unsigned pp = 1; pp < ns; ++pp )
    {
        d = sca[pp] * dif[pp];
        dst[DIM*pp  ] = src[DIM*pp  ] + d.XX - e.XX;
#if ( DIM > 1 )
        dst[DIM*pp+1] = src[DIM*pp+1] + d.YY - e.YY;
#endif
#if ( DIM > 2 )
        dst[DIM*pp+2] = src[DIM*pp+2] + d.ZZ - e.ZZ;
#endif
        e = d;
    }
    
    dst[DIM*ns  ] = src[DIM*ns  ] - d.XX;
#if ( DIM > 1 )
    dst[DIM*ns+1] = src[DIM*ns+1] - d.YY;
#endif
#if ( DIM > 2 )
    dst[DIM*ns+2] = src[DIM*ns+2] - d.ZZ;
#endif
}


/**
 Shorten segments to restore their length to 'cut'.
 We use a multidimensional Newton's method, to find iteratively the scalar
 coefficients that define the amount of displacement of each point.
 
 X[i] = vector of position
 We note 'dif' the differences between consecutive points:  dif[i] = X[i+1] - X[i]
 
 Given one scalar per segment: sca[i], the point is displaced as:
 Y[i] = X[i] + sca[i] * dif[i] - sca[i-1] * dif[i-1]
 except for the first and last points, for which there is only one term:
 Y[0] = X[0] + sca[  0] * dif[  0]
 Y[L] = X[L] - sca[L-1] * dif[L-1]
 
 We want 'sca' to restore the length of segments:
 ( Y[i+1] - Y[i] )^2 = cut^2
 
 i.e. 'sca' should fulfill a set of equalities F[i] = 0, with:
 F[i] = ( Y[i+1] - Y[i] )^2 - cut^2
 
 Method: use all zeros as first guess for 'sca', and apply a multidimensional
 Newton's method to iteratively refine the guess.
 
 In practice, we calculate `sca_next` from `sca` using the relationship:
 J(sca) * ( sca_next - sca ) = -F(sca)
 
 Where J is the Jacobian matrix: J[i,j] = dF[i] / dX[j]
 
 For this problem, J is square and tri-diagonal but not symmetric,
 and must be recalculated at each iteration.

 FJN, Strasbourg, 22 Feb 2015
 */
#if ( 1 )
int Filament::reshape_it(const unsigned ns, const real* src, real* dst, real cut)
{
    assert_true( ns > 1 );
    const real alphaSqr = cut * cut;
    unsigned cnt = 0;
    real err = 0;
    int info = 0;
    
    Vector * dif = new Vector[ns];
    
    // make a multiple of chunk to align memory:
    size_t chk = chunk_real(ns);

    //std::clog << "fiber::reshape_it allocates " << chk << std::endl;
    real * mem = new_real(chk*5);
    real * sca = mem;
    real * val = mem+chk;
    real * dia = mem+chk*2;
    real * low = mem+chk*3;
    real * upe = mem+chk*4;

    // calculate differences:
    if ( sizeof(Vector) == DIM * sizeof(real) )
    {
        real * dif_ = dif->data();
        const int end = DIM * ns;
        for ( int p = 0; p < end; ++p )
            dif_[p] = src[p+DIM] - src[p];
    }
    else
    {
        for ( unsigned p = 0; p < ns; ++p )
            dif[p] = diffPoints(src, p);
    }
    
    /*
     Perform here the first iteration of Newton's method
     the formula is the same as below, with all `sca` equal to zero,
     and thus 'vec == dif'
     The system is symmetric, and we can use a faster factorization
     */
    val[0] = dif[0].normSqr() - alphaSqr;
    dia[0] = 2 * dif[0].normSqr();
    for ( unsigned pp = 1; pp < ns; ++pp )
    {
        real n = dif[pp].normSqr();
        val[pp  ] = n - alphaSqr;
        dia[pp  ] = 2 * n;
        low[pp-1] = -dot(dif[pp], dif[pp-1]);
    }
    
    lapack::xpttrf(ns, dia, low, &info);
    if ( info ) {
        std::cerr << " reshape_it lapack::xpttrf failed " << info << std::endl;
        goto finish;
    }
    lapack::xptts2(ns, 1, dia, low, val, ns);
    
    err = 0;
    for ( unsigned pp = 0; pp < ns; ++pp )
    {
        sca[pp] = 0.5 * val[pp];
        err += fabs(val[pp]);
    }

#if ( 0 )
    printf("\n --- \n 0 sca ");
    for ( unsigned pp = 0; pp < ns; ++pp )
        printf(" %+6.4f", sca[pp]);
#endif

    while ( err > 1e-12 )
    {
        assert_true( ns > 1 );
        // set the matrix elements and RHS of system,
        // calculating 'vec' on the fly
        Vector vec0 = (1-2*sca[0])*dif[0] + sca[1]*dif[1];
        val[0] = vec0.normSqr() - alphaSqr;
        dia[0] = -2 * dot(vec0, dif[0]);
        upe[0] = dot(vec0, dif[1]);
        unsigned pp = 1;
        while ( pp+1 < ns )
        {
            vec0 = sca[pp-1]*dif[pp-1] + (1-2*sca[pp])*dif[pp] + sca[pp+1]*dif[pp+1];
            val[pp  ] = vec0.normSqr() - alphaSqr;
            dia[pp  ] = -2 * dot(vec0, dif[pp]);
            upe[pp  ] = dot(vec0, dif[pp+1]);
            low[pp-1] = dot(vec0, dif[pp-1]);
            ++pp;
        }
        assert_true( pp == ns-1 );
        vec0 = sca[pp-1]*dif[pp-1] + (1-2*sca[pp])*dif[pp];
        val[pp  ] = vec0.normSqr() - alphaSqr;
        dia[pp  ] = -2 * dot(vec0, dif[pp]);
        low[pp-1] = dot(vec0, dif[pp-1]);

#if ( 0 )
        real sum = 0;
        for ( unsigned pp = 0; pp < ns; ++pp )
            sum += val[pp];
        printf("\n   %i sum %8.5f", cnt, sum);
#endif
#if ( 0 )
        printf("\n   %i val  ", cnt);
        for ( unsigned pp = 0; pp < ns; ++pp )
            printf("%+6.4f ", val[pp]);
        printf("\n   %i upe  ", cnt);
        for ( unsigned pp = 0; pp+1 < ns; ++pp )
            printf("%+6.4f ", upe[pp]);
        printf("\n   %i dia  ", cnt);
        for ( unsigned pp = 0; pp < ns; ++pp )
            printf("%+6.4f ", dia[pp]);
        printf("\n   %i low  ", cnt);
        for ( unsigned pp = 0; pp < ns-1; ++pp )
            printf("%+6.4f ", low[pp]);
#endif
        
        lapack::xgtsv(ns, 1, low, dia, upe, val, ns, &info);
        if ( info )
        {
            std::cerr << " LAPACK dgtsv failed " << info << std::endl;
            goto finish;
        }

        // update `sca` and calculate residual error
        err = 0;
        for ( pp = 0; pp < ns; ++pp )
        {
            sca[pp] -= 0.5 * val[pp];
            err += fabs(val[pp]);
        }
        if ( ++cnt > 31 )
        {
            info = 1;
            goto finish;
        }
        
#if ( 0 )
        printf("\n %3i sca ", cnt);
        for ( pp = 0; pp < ns; ++pp )
            printf(" %+6.4f", sca[pp]);
#endif
    }
    
#if ( 0 )
    printf("\n%2i sca  ", cnt);
    for ( pp = 0; pp < ns; ++pp )
        printf("%+6.4f ", sca[pp]);
    printf("\n%2i err %e\n", cnt, err);
#endif
    
    //apply corrections:
    reshape_apply(ns, src, dst, sca, dif);
    
finish:
    free_real(mem);
    delete[] dif;
    
    return info;
}

#else

int Filament::reshape_it(const unsigned ns, const real* src, real* dst, real cut)
{
    assert_true( ns > 1 );
    const real alphaSqr = cut * cut;
    int info = 0;
    
    Vector * dif = new Vector[ns];
    Vector * vec = new Vector[ns];
    real * sca = new_real(ns);
    real * val = new_real(ns);
    real * dia = new_real(ns);
    real * low = new_real(ns);
    real * upe = new_real(ns);
    
    // calculate differences
    for ( unsigned pp = 0; pp < ns; ++pp )
    {
        dif[pp] = diffPoints(src, pp);
        sca[pp] = 0.0;
    }
    
    real err = 0;
    unsigned cnt = 0;
    do {
#if ( 0 )
        printf("\n   %i sca  ", cnt);
        for ( unsigned pp = 0; pp < ns; ++pp )
            printf("%+6.4f ", sca[pp]);
#endif

        // calculate all values of 'vec'
        vec[0] = (1-2*sca[0])*dif[0] + sca[1]*dif[1];
        for ( unsigned pp = 1; pp+1 < ns; ++pp )
            vec[pp] = sca[pp-1]*dif[pp-1] + (1-2*sca[pp])*dif[pp] + sca[pp+1]*dif[pp+1];
        vec[ns-1] = sca[ns-2]*dif[ns-2] + (1-2*sca[ns-1])*dif[ns-1];
        
        // calculate the matrix elements and RHS of system
        val[0] = vec[0].normSqr() - alphaSqr;
        dia[0] = -2 * ( vec[0] * dif[0] );
        for ( unsigned pp = 1; pp < ns; ++pp )
        {
            val[pp] = vec[pp].normSqr() - alphaSqr;
            dia[pp] = -2 * ( vec[pp] * dif[pp] );
            upe[pp-1] = vec[pp-1] * dif[pp];
            low[pp-1] = vec[pp] * dif[pp-1];
        }
        
#if ( 0 )
        printf("\n   %i val  ", cnt);
        for ( unsigned pp = 0; pp < ns; ++pp )
            printf("%+6.4f ", val[pp]);
        printf("\n   %i upe  ", cnt);
        for ( unsigned pp = 0; pp+1 < ns; ++pp )
            printf("%+6.4f ", upe[pp]);
        printf("\n   %i dia  ", cnt);
        for ( unsigned pp = 0; pp < ns; ++pp )
            printf("%+6.4f ", dia[pp]);
        printf("\n   %i low  ", cnt);
        for ( unsigned pp = 0; pp < ns-1; ++pp )
            printf("%+6.4f ", low[pp]);
#endif
        
        lapack::xgtsv(ns, 1, low, dia, upe, val, ns, &info);
        if ( info )
        {
            std::cerr << " LAPACK dgtsv failed " << info << std::endl;
            goto finish;
        }
        
        err = 0;
        for ( unsigned pp = 0; pp < ns; ++pp )
        {
            sca[pp] += -0.5 * val[pp];
            err += fabs(val[pp]);
        }
        if ( ++cnt > 32 )
        {
            info = 1;
            goto finish;
        }
    } while ( err > 0.0001 );

    
#if ( 0 )
    printf("\n%2i err %e", cnt, err);
    printf("\n%2i sca  ", cnt);
    for ( unsigned pp = 0; pp < ns; ++pp )
        printf("%+6.4f ", sca[pp]);
    printf("\n");
#endif
    
    //apply corrections:
    {
        Vector d, e = sca[0] * dif[0];
        dst[0] = src[0] + e.XX;
#if ( DIM > 1 )
        dst[1] = src[1] + e.YY;
#endif
#if ( DIM > 2 )
        dst[2] = src[2] + e.ZZ;
#endif
        for ( unsigned pp = 1; pp < ns; ++pp )
        {
            d = sca[pp] * dif[pp];
            dst[DIM*pp  ] = src[DIM*pp  ] + d.XX - e.XX;
#if ( DIM > 1 )
            dst[DIM*pp+1] = src[DIM*pp+1] + d.YY - e.YY;
#endif
#if ( DIM > 2 )
            dst[DIM*pp+2] = src[DIM*pp+2] + d.ZZ - e.ZZ;
#endif
            e = d;
        }
        dst[DIM*ns+0] = src[DIM*ns+0] - d.XX;
#if ( DIM > 1 )
        dst[DIM*ns+1] = src[DIM*ns+1] - d.YY;
#endif
#if ( DIM > 2 )
        dst[DIM*ns+2] = src[DIM*ns+2] - d.ZZ;
#endif
    }
    
finish:
    delete[] dif;
    delete[] vec;
    free_real(sca);
    free_real(val);
    free_real(dia);
    free_real(low);
    free_real(upe);
    return info;
}
#endif

/**
 The response of this method to a sudden perpendicular force is not ideal:
 For example, a force applied to the bottom of a vertical fibers leads
 to a 'L' configuration after one step of `solve()`.
 reshape() reduces the bottom leg of the 'L', by translating the entire vertical portion
 of the fiber, irrespective of the length of this section.
 */

#if ( 1 )   // 1 = optimized version of Filament::reshape_sure()

/**
 Move the vertices relative to each other, such that when this is done,
 all segments have the same distance `fnCut` ( =segmentation() ).
 This is operation does not change the center of gravity of the fiber.

 
 NOTE: if two consecutive points overlap, there is no unique way to
 restore the constraints! We do nothing in that case, because most 
 likely, the Brownian motion will push the points appart soon.
 */

void Filament::reshape_sure(const unsigned ns, real* vec, real cut)
{
    Vector dp(0,0,0), sum(0,0,0);
    Vector seg = diffPoints(vec, 0);
    real   dis = seg.norm();
    
    // translation needed to restore first segment
    if ( dis > REAL_EPSILON )
        dp = ( cut/dis - 1.0 ) * seg;
    
    for ( unsigned pp = 1; pp < ns; ++pp )
    {
        seg = diffPoints(vec, pp);
        dis = seg.norm();
        
        //move the left point by dp:
        dp.add_to(vec+DIM*pp);
        //update the uniform motion of the points:
        sum += dp;
        
        //add to the translation needed to restore this segment
        if ( dis > REAL_EPSILON )
            dp += ( cut/dis - 1.0 ) * seg;
    }
    
    //move the last point by dy[]:
    dp.add_to(vec+DIM*ns);
    
    // calculte a uniform motion to conserve the center of gravity:
    sum = ( sum + dp ) * ( -1.0 / ( ns + 1 ) );
    
    //translate the entire fiber uniformly:
    for ( unsigned pp = 0; pp <= ns; ++pp )
        sum.add_to(vec+DIM*pp);
}

#else

// ------------  old ( unoptimal ) version:
/**
 Move the vertices relative to each other, such that when this is done,
 all segments have the same distance segmentation() = fnCut.
 This is operation does not change the center of gravity of the fiber.
 */

void Filament::reshape_sure(const unsigned ns, real* vec, real cut)
{
    Vector dp, sum(0,0,0);
    
    for ( unsigned pp = 1; pp <= ns; ++pp )
    {
        dp       = diffPoints(vec, pp-1);
        real dis = dp.norm();
        if ( dis > REAL_EPSILON )
        {
            dp  *= ( cut/dis - 1.0 );
            for ( unsigned qq = pp; qq <= ns; ++qq )
                dp.add_to(vec+DIM*qq);
            sum += ( 1 + ns - pp ) * dp;
        }
    }
    
    sum *= ( -1.0 / (1+ns) );
    for ( unsigned pp = 0; pp <= ns; ++pp )
        sum.add_to(vec+DIM*pp);
}

#endif


void Filament::reshape()
{
    assert_true( nPoints > 1 );
#if ( DIM > 1 )
    if ( nPoints == 2 )
        reshape_two(pPos, pPos, fnCut);
    else if ( reshape_it(nbSegments(), pPos, pPos, fnCut) )
#endif
        reshape_sure(nbSegments(), pPos, fnCut);
}


void Filament::getPoints(const real * x)
{
#if ( DIM == 1 )
    Mecable::getPoints(x);
    reshape_sure(nbSegments(), pPos, fnCut);
#else
    if ( nPoints == 2 )
        reshape_two(x, pPos, fnCut);
    else if ( reshape_it(nbSegments(), x, pPos, fnCut) )
    {
        Mecable::getPoints(x);
        reshape_sure(nbSegments(), pPos, fnCut);
        //std::cerr << "A crude method was used to reshape " << reference() << '\n';
    }
#endif
    //dump(std::cerr);
}


/**
 Flip all the points. We do not change fnAscissa,
 and the abscissa of center thus stays as it is:
*/
void Filament::flipPolarity()
{
    unsigned ii = 0;
    unsigned jj = lastPoint();
    
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
 
 Note 2: Unless the Filament is straight, the length of the segments after this
 will not exactly match `segmentation()`.
*/
void Filament::growM(const real delta)
{
    assert_true( length() + delta > 0 );
    real a = -delta / length();
    
    if ( delta > 0 )
    {
        unsigned p = 0, n = nbSegments();
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
        for ( unsigned p = 0, n = nbSegments(); n > 0; ++p, --n )
            movePoint(p, ( a * n ) * diffPoints(p));
    }
    
    fnAbscissaM -= delta;
    setSegmentation(fnCut-a*fnCut);
    postUpdate();
}

/**
 This extends the fiber by adding one segment at the MINUS_END.
 Thus `segmentation()` is not changed, and the existing points are not displaced.
 */
void Filament::addSegmentM()
{
    unsigned pp = 1+nPoints;
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
void Filament::cutM(const real delta)
{
    real len = length();
    assert_true( 0 <= delta );
    assert_true( delta < len );
    
    const unsigned np = bestNumberOfPoints((len-delta)/fnSegmentation);
    const real cut = (len-delta) / (np-1);
    real* tmp = new_real(DIM*np);

    // calculate intermediate points into tmp[]:
    for ( unsigned pp=0; pp+1 < np; ++pp )
    {
        Vector w = interpolateM(delta+pp*cut).pos();
        w.store(tmp+DIM*pp);
    }
    
    // copy the position of plus-end into tmp[]:
    const unsigned lp = lastPoint();
    for ( unsigned d = 0 ; d < DIM; ++d )
        tmp[DIM*(np-1)+d] = pPos[DIM*lp+d];
    
    setNbPoints(np);
    
    // copy calculated points to pPos[]
    for ( unsigned pp = 0; pp < DIM*np; ++pp )
        pPos[pp] = tmp[pp];
    
    free_real(tmp);
    fnAbscissaM += delta;
    setSegmentation(cut);
    postUpdate();
}


/**
 The argument 'delta' can be positive or negative:
 - delta > 0 : elongation,
 - delta < 0 : shortening
 .
 
 Note 1: This works nicely if `delta` is small compared to segmentation().
 For large decrease in length, use cutP().

 Note 2: Unless the Filament is straight, the length of the segments after this
 will not exactly match `segmentation()`.
 */
void Filament::growP(const real delta)
{
    assert_true( length() + delta > 0 );
    real a = delta / length();
    
    if ( delta > 0 )
    {
        unsigned p = lastPoint();
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
        for ( unsigned p = lastPoint() ; p > 0 ; --p )
            movePoint(p, ( a * p ) * diffPoints(p-1));
    }
    
    setSegmentation(fnCut+a*fnCut);
    postUpdate();
}


/**
 This extends the fiber by adding one segment at the PLUS_END.
 Thus `segmentation()` is not changed, and the existing points are not displaced.
 */
void Filament::addSegmentP()
{
    unsigned pp = nPoints;
    setNbPoints(pp+1);
    
    real * psp = pPos + pp * DIM;
    for ( unsigned int dd = 0; dd < DIM; ++dd )
        psp[dd] = 2 * psp[dd-DIM] - psp[dd-2*DIM];
    
    postUpdate();
}


/**
 The Fiber length is reduced by `delta` ( which must be >= 0 ).
 The portion of size `delta` near the PLUS_END is removed,
 and the fewer vertices are recalculated.

 Note: after cutP(), the distance between the points is not exactly
 equal to segmentation(). This is true only if the fiber is straight.
*/
void Filament::cutP(const real delta)
{
    real len = length();
    assert_true( 0 <= delta );
    assert_true( delta < len );
    
    const unsigned np = bestNumberOfPoints((len-delta)/fnSegmentation);
    const real cut = (len-delta) / (np-1);
    real* tmp = new_real(DIM*np);
    
    // calculate intermediate points into tmp[]:
    for ( unsigned pp = 1; pp < np; ++pp )
    {
        Vector w = interpolateM(pp*cut).pos();
        w.store(tmp+DIM*pp);
    }
    
    setNbPoints(np);
    
    // copy calculated points to pPos[]
    // point at minus-end has not changed
    for ( unsigned pp = DIM; pp < DIM*np; ++pp )
        pPos[pp] = tmp[pp];
    
    free_real(tmp);
    setSegmentation(cut);
    postUpdate();
}

//------------------------------------------------------------------------------

void Filament::grow(FiberEnd end, const real delta)
{
    if ( end == PLUS_END )
        growP(delta);
    else if ( end == MINUS_END )
        growM(delta);
}


void Filament::adjustLength(real len, FiberEnd ref)
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


void Filament::truncateM(const unsigned int p)
{
    Mecable::truncateM(p);
    fnAbscissaM = abscissaPoint(p);
    postUpdate();
}


void Filament::truncateP(const unsigned int p)
{
    Mecable::truncateP(p);
    postUpdate();
}


/**
 `fib` is attached at the PLUS_END of `*this`
 
 The vertex are reinterpolated linearly, and the length of the
 segments will not fullfil the constraints of segmentation.
 If this is a problem, Filament::reshape() should be called.
 
 `fib` should usually be destroyed afterward.
 */
void Filament::join(Filament const* fib)
{
    const real len1 = length();
    const real lenT = len1 + fib->length();
    const unsigned np = bestNumberOfPoints(lenT/fnSegmentation) - 1;
    const real cut = lenT / real(np);
    
    // save position of PLUS_END:
    Vector ppe = fib->posEndP();
    
    real* tmp = new_real(DIM*np);
    
    // calculate new points into tmp[]:
    for ( unsigned pp = 1; pp < np; ++pp )
    {
        Vector w;
        if ( pp*cut < len1 )
            w = interpolateM(pp*cut).pos();
        else
            w = fib->interpolateM(pp*cut-len1).pos();
        
        w.store(tmp+DIM*pp);
    }
    
    setNbPoints(np+1);
    
    // copy point back in place:
    for ( unsigned int pp = DIM; pp < DIM*np; ++pp )
        pPos[pp] = tmp[pp];
    
    ppe.store(pPos+DIM*np);
    
    free_real(tmp);
    setSegmentation(cut);
    postUpdate();
}

//------------------------------------------------------------------------------
#pragma mark -


/**
 Returns the minimum and maximum distance between consecutive points
 */
void Filament::segmentationMinMax(real& mn, real& mx) const
{
    mn = diffPoints(0).norm();
    mx = mn;
    for ( unsigned n = 1; n < lastPoint(); ++n )
    {
        real r = diffPoints(n).norm();
        if ( r > mx )
            mx = r;
        if ( r < mn )
            mn = r;
    }
}

/**
 Returns the average and variances of segment length
 */
void Filament::segmentationVariance(real& avg, real& var) const
{
    avg = 0;
    var = 0;
    unsigned cnt = nbSegments();
    for ( unsigned n = 0; n < cnt; ++n )
    {
        real r = diffPoints(n).norm();
        avg += r;
        var += r*r;
    }
    avg /= cnt;
    var = var / cnt - avg * avg;
}

/**
 Calculate the inverse of the radius of the circle containing the points A, B, C
 
     cos(angle) = scalar_product( AB, BC ) / ( |AB| * |BC| )
     sin(angle) = sqrt( 1 - cos(angle)^2 )
     2 * radius * sin(angle) = |AC|
     curvature = 2 * sin(angle) / |AC|
     curvature = 2 * sqrt( ( 1 - cos(angle)^2 ) / AC^2 )

 */
real curvature3(Vector const& A, Vector const& B, Vector const& C)
{
    Vector ab = B - A;
    Vector bc = C - B;
    real P = dot(ab, bc);
    real S = 1.0 - ( P * P ) / ( ab.normSqr() * bc.normSqr() );
    real D = ( C - A ).normSqr();
    return 2.0 * sqrt( S / D );
}


real Filament::curvature(unsigned p) const
{
    assert_true( 0 < p && p < lastPoint() );
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
real Filament::bendingEnergy0() const
{
    real e = 0;
    
    const unsigned lsp = nPoints - 2;
    if ( lsp > 0 )
    {
        for ( unsigned p = 0; p < lsp ; ++p )
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


real Filament::minCosinus() const
{
    real result;
    Vector dir1, dir2;
    
    unsigned ps = nbSegments() % 2;
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
unsigned Filament::nbKinks(real threshold) const
{
    threshold *= fnCut * fnCut;
    unsigned res = 0;
    Vector d = diffPoints(0);
    
    for ( unsigned n = 1; n < lastPoint(); ++n )
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
 The position of the cut is `interpolatePoints(s, s+1, a)`
 */

real Filament::planarIntersect(unsigned s, Vector const& n, const real a) const
{
    assert_true( s < nbSegments() );
    
    real sca = dot(diffPoints(s), n);
    
    // if segment is parallel to plane, there is no intersection:
    if ( -REAL_EPSILON < sca  &&  sca < REAL_EPSILON )
        return INFINITY;
    
    Vector pos = posP(s);
    
    if ( modulo )
        modulo->fold(pos);
    
    return - ( dot(pos, n) + a ) / sca;
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 Recalculate the vertices of the Filament for 'ns' segments.
 
 @todo 2d-order interpolation in Filament::resegment()
 */
void Filament::resegment(unsigned ns)
{
    assert_true( ns > 0 );
    
    real cut = nbSegments() * fnCut / ns;
    
    // calculate new intermediate points in tmp[]:
    real* tmp = new_real(DIM*ns);
    Vector a = posP(0), b = posP(1);
    
    real h = 0;
    unsigned p = 1;
    
    for ( unsigned n = 1; n < ns; ++n )
    {
        h += cut;
        
        while ( h > fnCut )
        {
            h -= fnCut;
            a = b;
            ++p;
            assert_true(p<nPoints);
            b.load(pPos+DIM*p);
        }
        
        Vector w = a + ( h / fnCut ) * ( b - a );
        w.store(tmp+DIM*n);
    }
    
    // save index of PLUS_END
    p = DIM*lastPoint();
    
    // resize array:
    setNbPoints(ns+1);

    // move coordinates of last point
    for ( unsigned d = 0; d < DIM; ++d )
        pPos[DIM*ns+d] = pPos[p+d];

    // copy calculated coordinates back into pPos
    for ( unsigned d = DIM; d < DIM*ns; ++d )
        pPos[d] = tmp[d];

    free_real(tmp);
    setSegmentation(cut);
    /*
     Note: Unless the Filament is straight, the length of the segments after the
     interpolation will not exactly match `segmentation()`, and we therefore
     call reshape() to correct for the problem.
     */
    reshape();
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
void Filament::adjustSegmentation()
{
    assert_true( fnSegmentation > REAL_EPSILON );
    
    unsigned best = bestNumberOfPoints(nbSegments()*fnCut/fnSegmentation);
    
    if ( best != nPoints )
    {
        //std::clog << reference() << " resegment " << nPoints << " -> " << best << "\n";
#if ( 1 )
        resegment(best-1);
#else
        unsigned np = nPoints;
        // copy current points in temporary array:
        real* tmp = new_real(DIM*np);
        for ( int n = 0; n < DIM*np; ++n )
            tmp[n] = pPos[n];
        // re-interpolate:
        setShape(tmp, np, best);
        free_real(tmp);
#endif
    }
}


#else

static const real RECUT_PRECISION = 0.05;

/**
 A fiber is segmented as a function of its curvature.
 Only one new segmentation is done every time step at maximum,
 to allow the solver to equilibrate the new vertices.
 */
void Filament::adjustSegmentation()
{
    PRINT_ONCE("adjustSegmentation() is using the fiber curvature\n");
    
    const int upLimit = 8;
    unsgiend nbs = nbSegments();
    real len = nbs * fnCut;
    
    assert_true( fnSegmentation > REAL_EPSILON );
    //one segment for very short tubes
    if ( len <= fnSegmentation )
    {
        if ( nbs > 1 )
        {
            resegment(1);
            fnCutErrorIndex = 0;
        }
        return;
    }
    //at least two segments if length exceeds 2*fnSegmentation
    if ( nbs < 2 )
    {
        if ( len > fnSegmentation )
        {
            resegment(2);
            fnCutErrorIndex = 0;
        }
        return;
    }
    //the segment-length should not exceed a upper limit
    if ( fnCut >= upLimit * fnSegmentation )
    {
        resegment(2*nbs);
        fnCutErrorIndex = 0;
        return;
    }
    
    // accumulate the error in variable fnCutError, 
    // The error is in angle-square: 1-cos(angle) ~ (angle)^2
    if ( fnCutErrorIndex == 0 )
        fnCutError  =  1.0 - minCosinus();
    else
        fnCutError +=  1.0 - minCosinus();
    
    // after accumulation of RECUT_PERIOD time steps, we check the result
    if ( ++fnCutErrorIndex >= RECUT_PERIOD )
    {
        //we scale the error accumulated over the last steps
        fnCutError *= fnCut / ( RECUT_PRECISION * RECUT_PERIOD );
        
        if ( fnCutError > 1.0 )
        {
            //cut more finely if the error is large
            if ( fnCut > fnSegmentation )
                resegment(2*nbs);
        }
        else if ( fnCutError < 0.1 )
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
        
        //reset the counter and error accumulator
        fnCutErrorIndex = 0;
    }
}

#endif

//------------------------------------------------------------------------------
#pragma mark -

/**
 return the abscissa with respect to the ORIGIN.
 */
real Filament::abscissaEnd(const FiberEnd end) const
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
real Filament::abscissaFrom(const real dis, const FiberEnd ref) const
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
         attach1 = protein, 0.0, minus_end
         attach2 = protein, 0.0, plus_end
     }

*/
real Filament::someAbscissa(std::string const& key, Glossary& opt, real alpha) const
{
    const real len = length();
    real abs = len;

    if ( opt.set(abs, key, 1) )
    {
        FiberEnd ref = ORIGIN;
        opt.set(ref, key, 2, {{"plus_end", PLUS_END}, {"minus_end", MINUS_END}, {"center", CENTER}});
        
        int mod = 0;
        if ( opt.set(mod, key, 3, {{"off", 0}, {"uniform", 1}, {"exponential", 2}, {"regular", 3}}) )
        {
            if ( mod == 1 )
            {
                real a;
                do {
                    a = abs * RNG.preal();
                } while ( a > len );
                abs = a;
            }
            else if ( mod == 2 )
            {
                real a;
                do {
                    a = abs * RNG.exponential();
                } while ( a > len );
                abs = a;
            }
            else if ( mod == 3 )
                abs *= alpha;
        }
        
        abs = abscissaFrom(abs, ref);

        if ( !betweenMP(abs) )
            throw InvalidParameter("hand::abscissa is out of range");
        
        return abs;
    }

    // abscissa is set randomly:
    return RNG.real_uniform(abscissaM(), abscissaP());
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
FiberEnd Filament::whichEndDomain(const real ab, const real lambda) const
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


Mecapoint Filament::exactEnd(const FiberEnd end) const
{
    if ( end == MINUS_END )
        return Mecapoint(this, 0);
    else
    {
        assert_true( end == PLUS_END );
        return Mecapoint(this, lastPoint());
    }
}


Interpolation Filament::interpolateEnd(const FiberEnd end) const
{
    if ( end == MINUS_END )
        return interpolateEndM();
    else
    {
        assert_true( end == PLUS_END );
        return interpolateEndP();
    }
}


Interpolation Filament::interpolateCenter() const
{
    unsigned int n = lastPoint() / 2;
    if ( 2*n == lastPoint() )
        return Interpolation(this, n, n+1, 0);
    else
        return Interpolation(this, n, n+1, 0.5);
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
Interpolation Filament::interpolateM(const real ab) const
{
    if ( ab <= 0 )
        return Interpolation(this, 0, 1, 0.0);
    
    real a = ab / fnCut;
    unsigned s = (unsigned)a;
    
    //beyond the last point, we interpolate the PLUS_END
    if ( s+1 < nPoints )
        return Interpolation(this, s, s+1, a-s);
    else
        return Interpolation(this, nPoints-2, nPoints-1, 1.0);
}


Interpolation Filament::interpolate(const real ab, const FiberEnd end) const
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
Vector Filament::posM(const real ab) const
{
    // return MINUS_END
    if ( ab <= 0 )
        return posP(0);
    
    real a = ab / fnCut;
    unsigned s = (unsigned)a;
    
    // check if PLUS_END is reached:
    if ( s+1 < nPoints )
        return interpolatePoints(s, s+1, a-s);
    else
        return posP(lastPoint());
}

Vector Filament::dirM(const real ab) const
{
    // at MINUS_END
    if ( ab <= 0 )
        return dirSegment(0);
    
    real a = ab / fnCut;
    unsigned s = (unsigned)a;
    
    // check if PLUS_END is reached
    if ( s+1 < nPoints )
        return dirSegment(s);
    else
        return dirSegment(lastSegment());
}
#endif


Vector Filament::posEnd(FiberEnd end) const
{
    if ( end == MINUS_END )
        return posEndM();
    else if ( end == PLUS_END )
        return posEndP();
    else
        return posM(abscissaFrom(0, end));
}


Vector Filament::dirEnd(const FiberEnd end) const
{
    if ( end == MINUS_END )
        return dirSegment(0);
    else if ( end == PLUS_END )
        return dirSegment(lastSegment());
    else
        return dirM(abscissaFrom(0, end));
}



/// force on the PLUS_END projected on the direction of elongation
real Filament::projectedForceEndM() const
{
    return -dot(netForce(0), dirSegment(0));
}

/// force on the PLUS_END projected on the direction of elongation
real Filament::projectedForceEndP() const
{
    unsigned p = lastSegment();
    return dot(netForce(p+1), dirSegment(p));
}


/**
 The returned value is negative when the force antagonizes elongation,
 and this is true at both ends. 
 */
real Filament::projectedForceEnd(const FiberEnd end) const
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

int Filament::checkLength(real len) const
{
    assert_small( length() - len );
    real con = contourLength(pPos, nPoints);
    if ( fabs( con - len ) > 0.1 )
    {
        Cytosim::log << "Warning: length of " << reference() << " is " << con << " but " << len << " was expected" << std::endl;
        return 1;
    }
    return 0;
}


int Filament::checkSegments() const
{
    real mn, mx;
    segmentationMinMax(mn, mx);
    real d = mx - mn;
    if ( d > 1e-3 )
    {
         std::clog << "Fiber " << reference() << " segments [ " << std::fixed << mx << " " << std::fixed << mx;
         std::clog << " ] deviate from segmentation " << segmentation() << " by " << d << std::endl;
         return 1;
    }
    return 0;
}


/**
 Prints info on the length of Segments, which can be useful for debugging
 */
void Filament::dump(std::ostream& os) const
{
    os << "Fiber " << std::setw(7) << reference();
    os << "  " << std::left << std::setw(6) << fnCut << " {";
    
#if ( 1 )
    real mn, mx;
    segmentationMinMax(mn, mx);
    real p = 100 * ( mx - mn ) / mn;
    os.precision(4);
    os << " " << mn << " + " << std::setw(7) << std::fixed << p << " %";
#else
    for ( size_t pp = 0; pp < lastPoint(); ++pp )
    {
        real p = 100 * ( diffPoints(pp).norm() / fnCut - 1.0 );
        os << " " << std::left << std::setw(5) << p;
    }
#endif
    os << " }" << std::endl;
}


void Filament::write(Outputter& out) const
{
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
void Filament::read(Inputter & in, Simul& sim, ObjectTag tag)
{
    //Cytosim::log << "  reading Filament at " << in.pos() << '\n';
    
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

    // resegment if the sementation parameter has changed:
    if ( fnSegmentation != seg )
        adjustSegmentation();
    
    //Mecable::write(std::cerr);
    
    // verify the length and segmentation:
    if ( in.vectorSize() == DIM )
    {
        checkLength(len);
        checkSegments();
    }
}

