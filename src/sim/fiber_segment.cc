// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "fiber_segment.h"
#include "map.h"  // for GRID_HAS_PERIODIC
#include "space.h"
#include "fiber.h"
#include "modulo.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------
//---------------- DISTANCE FROM A POINT TO A SECTION OF FIBER -----------------
//------------------------------------------------------------------------------

/**
 W is projected on the line supporting this FiberSegment
 The function calculates:
 - abs <- the signed distance from pos1() to the projection of W
 - dis <- the SQUARE of the distance between W and its projection
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
#elif ( DIM == 2 )
    dis = W.normSqr() - abs * abs;
#else
    dis = ( W.XX * W.XX + W.YY * W.YY ) + ( W.ZZ * W.ZZ - abs * abs );
#endif

    return abs;
}


/**
 W is projected on the line that supports this FiberSegment
 The function calculates:
 - abs <- the signed distance from pos1() to the projection of W
 - dis <- the SQUARE of the distance between W and its projection
 .
 
 It is assumed here that len() returns the distance between the two points of the FiberSegment
 Attention: `dis` may NOT be set if ( abs < 0 ) or ( abs > len() )
 */
real FiberSegment::projectPoint(Vector W, real& dis) const
{
    Vector A = pos1();
    Vector D = diff();
    W -= A;
    
#if GRID_HAS_PERIODIC
    if ( modulo )
        modulo->fold(W);
#endif
    
    // project with the scalar product:
    real abs = dot(W, D) * lenInv();
    
    // test boundaries of filament:
    if ( abs < 0 )
    {
        if ( isFirst() )
            dis = W.normSqr();
    }
    else if ( abs > len() )
    {
        if ( isLast() )
            dis = (W-D).normSqr();
    }
    else
    {
#if ( DIM == 1 )
        dis = 0;
#elif ( DIM == 2 )
        dis = W.normSqr() - abs * abs;
#else
        // must optimize explicitely, since the compiler cannot reorder additions
        dis = ( W.XX * W.XX + W.YY * W.YY ) + ( W.ZZ * W.ZZ - abs * abs );
#endif
    }
    return abs;
}


/**
 This may be faster than projectPoint(), but it does not work with periodic boundaries
 The function calculates:
 - abs <- the signed distance from pos1() to the projection of W
 - dis <- the SQUARE of the distance between W and its projection
 .
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
#else
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
        dis = aX * aX + ( aY * aY - abs * abs );
#else
        dis = ( aX * aX + aY * aY ) + ( aZ * aZ - abs * abs );
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

 @sets arguments `abs1` and `abs2` to be the abscissa of the corresponding points
 @returns `distance(P1, P2)^2`, the SQUARE of the distance BETWEEN THE LINES, and
 not the segments. In 2D, `dis2 = 0` unless the lines are parallel.

 If the segments are parallel, P1 and P2 are set to the middle of the overlapping
 section between the two segments.

 The calling function must check if 'abs1' and 'abs2' are within [0, segment_length]
 to do anything meaningful.
 */

real FiberSegment::shortestDistanceSqr(FiberSegment const& seg, real& abs1, real& abs2) const
{
    Vector off = seg.pos1() - pos1();
    Vector d11 = dir();
    Vector d22 = seg.dir();
    
#if GRID_HAS_PERIODIC
    if ( modulo )
        modulo->fold(off);
#endif
    
#if ( DIM > 2 )
    Torque N = cross(d11, d22);
    real S = N.normSqr();  // sine^2
#else
    real S = square(cross(d11, d22));  // sine^2
#endif
    
    if ( abs_real(S) > 128 * REAL_EPSILON )
    {
        // This deals with the general case of non-parallel lines
        real iS = 1 / S;    // 1.0 / sine^2
#if ( DIM > 2 )
        real DD = square( dot(off, N) ) * iS;
#endif
        /*
         We do not necessarily need to calculate the projections `abs1 , abs2`
         if the distance is greater than the threshold above which nothing
         is do be done with these results...
         */
        real C = dot(d11, d22);  // cosine of angle
        real d1off = dot(d11, off);
        real d2off = dot(d22, off);

        abs1 = ( d1off - C * d2off ) * iS;
        abs2 = ( C * d1off - d2off ) * iS;
#if 0
        // check that identified line path is orthogonal to both segments:
        Vector p1 = pos1() + abs1 * d11;
        Vector p2 = seg.pos1() + abs2 * d22;
        real n1 = dot(d11, p2-p1);
        real n2 = dot(d22, p2-p1);
        printf("belowDistance %+9.6f %+9.6f :", n1, n2);
        // check different formula for distance betwen lines:
        real res0 = ( off + abs2 * d22 - abs1 * d11 ).norm();
        real res1 = sqrt( D * iS );  // 1.0 / N.normSqr() == iS
        real res2 = ( p2 - p1 ).norm();
        printf("%6.4f  %6.4f  %6.4f\n", res0, res1, res2);
#endif
#if ( DIM > 2 )
        // direction N of the shortest path is orthogonal to both lines, thus
        // distance between lines = dot(off, N) / N.norm()
        return DD;
#endif
        return 0;
    }

    /*
     This deals with the case where the two segments are almost parallel:
     S ~= 0; C ~= +/- 1
     m1 = projection of seg.pos1() on this segment
     p1 = projection of seg.pos2()
     */
    const real d1off = dot(d11, off);
    const real len1 = len(), len2 = seg.len();
    
    // the distance between lines
    real DD = off.normSqr() - d1off * d1off;
    //real E = ( off + abs2 * d22 - abs1 * d11 ).normSqr();

    real C = dot(d11, d22);  // cosine of angle
    real m1 = d1off;
    real p1 = d1off + C * len2;
    real m2 = -dot(d22, off);
    real p2 = m2 + C * len1;

    // clamp inside segment and use mid-point
    abs1 = 0.5 * ( min_real(len1, max_real(m1, p1)) + max_real(0, min_real(m1, p1)));
    
    // clamp inside segment and use mid-point
    abs2 = 0.5 * ( min_real(len2, max_real(m2, p2)) + max_real(0, min_real(m2, p2)));
    
    // return distance between lines
    return DD;
}


void FiberSegment::print(std::ostream& os) const
{
    if ( fiber() )
        os << "(" << fiber()->reference() << " " << point() << ")";
    else
        os << "(null)";
}

std::string FiberSegment::to_string() const
{
    std::ostringstream oss;
    print(oss);
    return oss.str();
}

std::ostream& operator << (std::ostream& os, FiberSegment const& arg)
{
    arg.print(os);
    return os;
}
