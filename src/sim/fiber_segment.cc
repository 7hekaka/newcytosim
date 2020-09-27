// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "fiber_segment.h"
#include "space.h"
#include "fiber.h"
#include "simul.h"
#include "modulo.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------
//---------------- DISTANCE FROM A POINT TO A SECTION OF FIBER -----------------
//------------------------------------------------------------------------------

/**
 W is projected on the line supporting this FiberSegment
 The function calculates:
 - abs <- the signed distance from pos1() to the projection of W
 - dis <- the distance between W and its projection
 .
 
 This uses FiberSegment::lenInv() that should return 1.0 / segment_length
 */
real FiberSegment::projectPoint0(Vector W, real& dis) const
{
    Vector A = pos1();
    W -= A;
    
#if GRID_HAS_PERIODIC
    if ( modulo )
        modulo->fold(W);
#endif
    
    // project with the scalar product:
    real abs = dot(W, pos2()-A) * lenInv();
    
    // calculate distance to projection
#if ( DIM == 1 )
    dis = 0;
#else
    dis = W.normSqr() - abs * abs;
#endif

    return abs;
}


/*
 Assuming that the information in iDir[] is up-to-date
 */
real FiberSegment::projectPoint1(Vector W, real& dis) const
{
    //assert_true(fib_->iDirValid);
    W -= pos1();
    
#if GRID_HAS_PERIODIC
    if ( modulo )
        modulo->fold(W);
#endif
    
    // project with the scalar product:
    real abs = dot(W, dirS());
    
    // calculate distance to projection
#if ( DIM == 1 )
    dis = 0;
#else
    dis = W.normSqr() - abs * abs;
#endif

    return abs;
}


/**
 W is projected on the line that supports this FiberSegment
 The function calculates:
 - abs <- the signed distance from pos1() to the projection of W
 - dis <- the distance between W and its projection
 .
 
 It is assumed here that len() returns the distance between the two points of the FiberSegment
 Attention: `dis` may NOT be set if ( abs < 0 ) or ( abs > len() )
 */
real FiberSegment::projectPoint(Vector const& W, real& dis) const
{
    Vector A = pos1();
    Vector AW = W - A;
    
#if GRID_HAS_PERIODIC
    if ( modulo )
        modulo->fold(AW);
#endif
    
    // project with the scalar product:
    real abs = dot(AW, pos2()-A) * lenInv();
    
    // test boundaries of filament:
    if ( abs < 0 )
    {
        if ( isFirst() )
            dis = distanceSqr(W, A);
    }
    else if ( abs > len() )
    {
        if ( isLast() )
            dis = distanceSqr(W, pos2());
    }
    else
    {
#if ( DIM == 1 )
        dis = 0;
#else
        dis = AW.normSqr() - abs * abs;
#endif
    }
    return abs;
}


/**
 This may be faster than projectPoint(), but it does not work with periodic boundaries
 */
real FiberSegment::projectPointF(const real w[], real& dis) const
{
    assert_true( !modulo );
    
    const real * p = fib_->addrPoint(pti_);
    
    real dX = p[DIM  ] - p[0];
    real aX = w[0]     - p[0];
#if ( DIM > 1 )
    real dY = p[DIM+1] - p[1];
    real aY = w[1]     - p[1];
#endif
#if ( DIM > 2 )
    real dZ = p[DIM+2] - p[2];
    real aZ = w[2]     - p[2];
#endif
    
    // project with the scalar product:
#if ( DIM == 1 )
    real abs = ( dX * aX ) * lenInv();
#elif ( DIM == 2 )
    real abs = ( dX * aX + dY * aY ) * lenInv();
#elif ( DIM == 3 )
    real abs = ( dX * aX + dY * aY + dZ * aZ ) * lenInv();
#endif
    
    // test boundaries of segment:
    if ( abs < 0 )
    {
        if ( isFirst() )
            dis = distanceSqr(w, pos1());
    }
    else if ( abs > len() )
    {
        if ( isLast() )
            dis = distanceSqr(w, pos2());
    }
    else
    {
#if   ( DIM == 1 )
        dis = 0;
#elif ( DIM == 2 )
        dis = aX * aX + aY * aY - abs * abs;
#elif ( DIM == 3 )
        dis = aX * aX + aY * aY + aZ * aZ - abs * abs;
#endif
        
#if ( 0 )
        // verify that the results are identical to projectPoint()
        real d = dis;
        real a = projectPoint(Vector(w), d);
        assert_small(a-abs);
        assert_small(d-dis);
#endif
    }
    return abs;
}


//------------------------------------------------------------------------------
//---------------- DISTANCE TO ANOTHER SECTION OF A FIBER ----------------------
//------------------------------------------------------------------------------

/**
 Evaluate the minimal distance between two segments.
 
 This finds the positions P1, P2 on the two supporting lines that are closest
 to each other. P1 belongs to `*this`, and P2 to `seg`. These points are defined
 by their abscissa (abs1, abs2) for which the values [0, segment_length] match
 the edges of the segments.
 
 @return `distance(P1, P2)^2', the square of the distance between the lines.
 @sets arguments `abs1` and `abs2` to the corresponding points
 
 The calling function must check if 'abs1' and 'abs2' are within [0, segment_length]
 to do any meaningful work.

 If the segments are parallel, P1 and P2 are set to the middle of the overlapping
 section between the two segments.
 */

real FiberSegment::shortestDistance(FiberSegment const& seg, real& abs1, real& abs2) const
{
    const real len1 = len();
    const real len2 = seg.len();

    Vector off = seg.pos1() - pos1();
#if 0
    assert_true(fib_->iDirValid);
    Vector d11 = dirS();
    Vector d22 = seg.dirS();
#else
    Vector d11 = dir();
    Vector d22 = seg.dir();
#endif
    
#if GRID_HAS_PERIODIC
    if ( modulo )
        modulo->fold(off);
#endif
    
    real C = dot(d11, d22);  // cosinus of angle
    real S = 1.0 - C * C;    // sinus squared

    if ( S > REAL_EPSILON )
    {
        real iS = 1.0 / S;
        // This deals with the general case of non-parallel lines
        real d1off = dot(d11, off) * iS;
        real d2off = dot(d22, off) * iS;
        
        //abs1 = dot(d11-beta*d22, off) / scal;
        abs1 = d1off - C * d2off;
        
        //abs2 = dot(beta*d11-d22, off) / scal;
        abs2 = C * d1off - d2off;
#if 0
        // check that identified line path is orthogonal to both segments:
        Vector p1 = pos1() + abs1 * d11;
        Vector p2 = seg.pos1() + abs2 * d22;
        real n1 = dot(d11, p2-p1);
        real n2 = dot(d22, p2-p1);
        real res1 = std::sqrt(( off - d11 * abs1 ).normSqr() - abs2 * abs2);
        real res2 = std::sqrt(( off + d22 * abs2 ).normSqr() - abs1 * abs1);
        printf("shortestDistance %+9.6f %+9.6f   %6.3f   %6.3f ", n1, n2, res1, res2);
#endif
#if ( DIM > 2 )
        //real res = ( off + abs2 * d22 - abs1 * d11 ).normSqr();
        real res = ( off - d11 * abs1 ).normSqr() - abs2 * abs2;
        //real res = ( off + d22 * abs2 ).normSqr() - abs1 * abs1;
        //printf("dis %03u:%02lu", fiber()->identity(), point());
        //printf(" %03u:%02lu   %6.3f", seg.fiber()->identity(), seg.point(), std::sqrt(res));
        return res;
#else
        return 0;
#endif
    }

    /*
     This deals with the case where the two segments are almost parallel:
     S ~= 0; C ~= +/- 1
     m1 = projection of seg.pos1() on this segment
     p1 = projection of seg.pos2()
     */
    const real d1off = dot(d11, off);
    
    real m1 = d1off;
    real p1 = d1off + C * len2;

    // clamp inside segment and use mid-point
    abs1 = 0.5 * ( std::min(len1, std::max(m1, p1)) + std::max((real)0, std::min(m1, p1)));
    
    real m2 = -dot(d22, off);
    real p2 = m2 + C * len1;

    // clamp inside segment and use mid-point
    abs2 = 0.5 * ( std::min(len2, std::max(m2, p2)) + std::max((real)0, std::min(m2, p2)));

    // return distance between lines
    return off.normSqr() - d1off * d1off;
}


void FiberSegment::print(std::ostream& os) const
{
    if ( fiber() )
        os << "(" << fiber()->reference() << " " << std::setw(3) << point() << ":)";
    else
        os << "(null)";
}


std::ostream& operator << (std::ostream& os, FiberSegment const& arg)
{
    arg.print(os);
    return os;
}
