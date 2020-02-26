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
 W is projected on the line that supports this FiberSegment
 The function calculates:
 - abs <- the signed distance from pos1() to the projection of W
 - dis <- the distance between W and its projection
 .
 
 It is assumed here that len() returns the distance between the two points of the FiberSegment
 Attention: `dis` is not set if ( abs < 0 ) or ( abs > len() )
 */
real FiberSegment::projectPoint0(Vector aw, real& dis) const
{
    assert_true( fib_ );
    
    Vector A = pos1();
    aw -= A;
    
    if ( modulo )
        modulo->fold(aw);
    
    // project with the scalar product:
    real abs = dot(aw, pos2()-A) * lenInv();
    
    // calculate distance to projection
#if ( DIM == 1 )
    dis = 0;
#else
    dis = aw.normSqr() - abs * abs;
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
real FiberSegment::projectPoint(Vector const& w, real& dis) const
{
    assert_true( fib_ );
    
    Vector A = pos1();
    Vector aw = w - A;
    
    if ( modulo )
        modulo->fold(aw);
    
    // project with the scalar product:
    real abs = dot(aw, pos2()-A) * lenInv();
    
    // test boundaries of filament:
    if ( abs < 0 )
    {
        if ( isFirst() )
            dis = distanceSqr(w, A);
    }
    else if ( abs > len() )
    {
        if ( isLast() )
            dis = distanceSqr(w, pos2());
    }
    else
    {
#if ( DIM == 1 )
        dis = 0;
#else
        dis = aw.normSqr() - abs * abs;
#endif
    }
    return abs;
}


/**
 This may be faster than projectPoint(), but it does not work with periodic boundaries
 */
real FiberSegment::projectPointF(const real w[], real& dis) const
{
    assert_true( fib_ );
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
 
 This finds the positions P1, P2 connecting the two supporting lines in the shortest
 way possible. P1 belongs to *this, and P2 to `seg`. These points are defined by
 their abscissa (abs1, abs2) for which the values {0, segment_length} match the
 edges of the segments.
 
 The function then evaluates if the points P1, P2 are inside or outside their 
 respective segments.
 
 @return `INFINITY` if `P1` or `P2` are outside their respective segment.
 
 @return `distance squared` if P1 and P2 are both inside their segments.
 In this case, `dis` is set to be the square of the distance (P1, P2)
 Hence `dis`, `abs1` and `abs2` are only set valid if the return value is `1`.

 If the segments are parallel, the mid-point of the overlapping section is returned.
 */

real FiberSegment::shortestDistance(FiberSegment const& seg, real& abs1, real& abs2) const
{
    const real len1 = len();
    const real len2 = seg.len();

    Vector pos = pos1();
    Vector off = seg.pos1();
    Vector d11 = ( pos2() - pos ) * lenInv();
    Vector d22 = ( seg.pos2() - off ) * seg.lenInv();
    off -= pos; // off = seg.pos1() - pos1()
    
    if ( modulo )
        modulo->fold(off);
    
    real beta = dot(d11, d22);
    real scal = 1.0 - beta * beta;
    
    if ( scal > REAL_EPSILON )
    {
        // This is the general case of non-parallel lines:
        real d1off = dot(d11, off) / scal;
        real d2off = dot(d22, off) / scal;
        
        //abs1 = dot(d11-beta*d22, off) / scal;
        abs1 = d1off - beta * d2off;
        
        //abs2 = dot(beta*d11-d22, off) / scal;
        abs2 = beta * d1off - d2off;
        
#if 0
        // check that identified line path is orthogonal to both segments:
        Vector p1 = pos1() + abs1 * dir();
        Vector p2 = seg.pos1() + abs2 * seg.dir();
        real n1 = dot(d11, p2-p1);
        real n2 = dot(d22, p2-p1);
        printf("shortestDistance %+9.6f %+9.6f\n", n1, n2);
#endif

        if ( abs1 < 0 | len1 <= abs1 | abs2 < 0 | len2 <= abs2 )
            return INFINITY;
        
#if ( DIM > 2 )
        //dis = ( off + abs2 * d22 - abs1 * d11 ).normSqr();
        //dis = ( off - d11 * abs1 ).normSqr() - abs2 * abs2;
        return ( off + d22 * abs2 ).normSqr() - abs1 * abs1;
#else
        return 0;
#endif
    }
    else
    {
        /*
         This deals with the case where the two segments are almost parallel:
         beta ~ +/- 1
         m1 = projection of seg.pos1() on this segment
         p1 = projection of seg.pos2()
         */
        const real d1off = dot(d11, off);
        
        real m1 = d1off;
        real p1 = d1off + beta * len2;
        
        if (( m1 < 0 & p1 < 0 ) | ( m1 > len1 & p1 > len1 ))
            return INFINITY;

        // clamp inside segment and take mid-point:
        abs1 = 0.5 * ( clamp(m1, 0, len1) + clamp(p1, 0, len1) );

        real m2 = -dot(d22, off);
        real p2 = m2 + beta * len1;
        
        if (( m2 < 0 & p2 < 0 ) | ( m2 > len2 & p2 > len2 ))
            return INFINITY;

        // clamp inside segment and take mid-point:
        abs2 = 0.5 * ( clamp(m2, 0, len2) + clamp(p2, 0, len2) );
        
        // return distance between
        return off.normSqr() - d1off * d1off;
    }
    
    return 1;
}


void FiberSegment::print(std::ostream& os) const
{
    if ( fiber() )
        //os << "(" << fiber()->reference() << " seg " << point() << ":" << point()+1 << ")";
        os << "(f" << fiber()->identity() << " " << point() << ":" << point()+1 << ")";
    else
        os << "(null)";
}


std::ostream& operator << (std::ostream& os, FiberSegment const& obj)
{
    obj.print(os);
    return os;
}
