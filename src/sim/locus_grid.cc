// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "assert_macro.h"
#include "locus_grid.h"
#include "mecapoint.h"
#include "simd_float.h"
#include "fiber_segment.h"
#include "exceptions.h"
#include "messages.h"
#include "modulo.h"
#include "space.h"
#include "solid.h"
#include "meca.h"
#include "simd.h"

//------------------------------------------------------------------------------

index_t LocusGrid::setGrid(Vector inf, Vector sup, real min_width)
{
    assert_true( min_width > 0 );
    
    index_t cnt[3] = { 1, 1, 1 };
    for ( int d = 0; d < DIM; ++d )
    {
        // minimum number of cells in dimension 'd'
        int n = std::floor(( sup[d] - inf[d] ) / min_width);
        
        if ( n < 0 )
            throw InvalidParameter("invalid space:boundaries");
        
        if ( modulo && modulo->isPeriodic(d) )
        {
            assert_small( modulo->period_[d] - sup[d] + inf[d] );
            // adjust the grid to match the edges
            cnt[d] = std::max(1, n);
            pGrid.setPeriodic(d, true);
        }
        else
        {
            //extend in any dimension that is not periodic, adjusting cell size to min_width
            cnt[d] = n + 1;
            real w = 0.5 * cnt[d] * min_width;
            real m = 0.5 * ( inf[d] + sup[d] );
            inf[d] = m - w;
            sup[d] = m + w;
        }
    }
    
    pGrid.setDimensions(inf, sup, cnt);
    return pGrid.nbCells();
}


void LocusGrid::createCells()
{
    if ( pGrid.createCells() )
        pGrid.printSummary(Cytosim::log, "   LocusGrid");

    // create side regions suitable for pairwise interactions:
    if ( !pGrid.hasRegions() )
        pGrid.createSideRegions(1);
}


size_t LocusGrid::capacity() const
{
    size_t res = 0;
    for ( index_t i = 0; i < pGrid.nbCells(); ++i )
        res += pGrid[i].capacity();
    return res;
}


//------------------------------------------------------------------------------
#pragma mark - Check two Objects: P = Point; L = Line segment

/// used to check distance of two particles (dd = square of distance) against threshold L
static inline bool below(const real dd, const real L) { return ( dd < L*L ) && ( dd > REAL_EPSILON ); }


/**
 This is used to check two spherical objects:
 Solid/Bead/Sphere or the terminal vertex (the tips) of a Fiber
 
 The force is applied if the objects are closer to the maximum
 of their specified range + radius.
 */
void LocusGrid::checkPP(BigPoint const& aa, BigPoint const& bb) const
{
    //std::clog << "   PP- " << bb.mec_ << " " << aa.mec_ << '\n';
    assert_true( aa.obj_->tag() != Fiber::TAG );
    assert_true( bb.obj_->tag() != Fiber::TAG );

    Vector vab = bb.cen() - aa.cen();
    const float ran = aa.rad_ + bb.rad_;

    if ( modulo )
        modulo->fold(vab);
    real ab2 = vab.normSqr();
    
    if ( below(ab2, ran) )
        meca_.addLongLink1(aa.vertex1(), bb.vertex1(), vab, ab2, ran, push_);
}


/**
 This is used to check a segment of a fiber against a spherical object:
 Solid/Bead/Sphere or Fiber-tip.
 
 The force is applied if the objects are closer than the maximum of the two range + radius,
 and if the center of the sphere projects inside the segment.
 */
void LocusGrid::checkPL(BigPoint const& aa, BigLocus const& bb) const
{
    //std::clog << "   PL- " << aa.vertex1() << " " << bb.segment() << '\n';
    assert_true( aa.obj_->tag() != Fiber::TAG );
    assert_true( bb.obj_->tag() == Fiber::TAG );

    const float ran = aa.rad_ + bb.rad_;
    
    // determine projection of `aa` on segment `bb`:
    real ab2 = INFINITY;
    FiberSegment seg = bb.segment();
    real abs = seg.projectPoint0(aa.cen(), ab2);
    
    if ( 0 <= abs )
    {
        if ( abs <= bb.seg() )
        {
            // the point projects inside the segment
            if ( below(ab2, ran) )
                meca_.addSideSlidingLink(seg, abs, aa.vertex1(), ran, push_);
        }
        else
        {
            // the point projects right past the segment, and we check the fiber tip
            if ( bb.isLast() )
            {
                Vector vab = bb.pos2() - aa.cen();
                if ( modulo )
                    modulo->fold(vab);
                ab2 = vab.normSqr();
                if ( below(ab2, ran) )
                    meca_.addLongLink1(aa.vertex1(), bb.vertex2(), vab, ab2, ran, push_);
            }
        }
    }
    else
    {
        // here abs < 0, and thus `bb` projects left past the segment
        Vector vab = bb.pos1() - aa.cen();
        if ( modulo )
            modulo->fold(vab);
        ab2 = vab.normSqr();
        /*
         This code handles the interactions will all the joints of a fiber:
         interact with the node only if this projects on the previous segment
         or if this is the terminal point of a fiber.
         */
        if ( below(ab2, ran) && ( bb.isFirst() || dot(vab, bb.prevDiff()) <= 0 ))
            meca_.addLongLink1(aa.vertex1(), bb.vertex1(), vab, ab2, ran, push_);
    }
}


/**
 This is used to check vertex 1 of segment 'aa' against vertex 1 of segment 'bb'
 
 The interaction is applied only if the vertex projects 'inside' the segment.
 */
void LocusGrid::checkLLP1(BigLocus const& aa, BigLocus const& bb, float ran, real abs, Vector const& vab) const
{
    //std::clog << "   LLP1 " << aa.vertex1() << " " << bb.vertex1() << '\n';

    //const float ran = aa.rad_ + bb.rad_;
    if ( aa.isFirst() )
    {
        /*
         Check the projection of aa.vertex1(),
         on the segment represented by 'bb'
         */
        if ( &bb < &aa  &&  bb.isFirst() )
        {
            real ab2 = vab.normSqr();
            if ( below(ab2, ran)  &&  dot(vab, bb.diff()) >= 0 )
                meca_.addLongLink1(aa.vertex1(), bb.vertex1(), vab, ab2, ran, push_);
        }
    }
    else
    {
        /*
         Check the projection to the segment located before 'aa',
         and interact if 'bb.vertex1()' falls on the right side of it
         */
        if ( dot(vab, aa.prevDiff()) >= 0 )
        {
            real ab2 = vab.normSqr();
            if ( below(ab2, ran) )
                meca_.addLongLink1(aa.vertex1(), bb.vertex1(), vab, ab2, ran, push_);
        }
    }
}

/**
 This is used to check vertex 1 of segment 'aa' against vertex 2 of segment 'bb'
 and vertex 2 of segment 'aa' against vertex 2 of segment 'bb'

 The interaction is applied only if the vertex projects 'inside' the segment.
 */
void LocusGrid::checkLLP2(BigLocus const& aa, BigLocus const& bb, float ran, real abs, Vector const& vab) const
{
    //std::clog << "   LLP2 " << aa.vertex1() << " " << bb.vertex2() << '\n';
    
    //const real ran = aa.rad_ + bb.rad_;
    if ( abs < 0 )
    {
        if ( aa.isFirst() )
        {
            assert_true(bb.isLast());
            real ab2 = vab.normSqr();

            if ( below(ab2, ran)  && dot(vab, bb.diff()) <= 0 )
                meca_.addLongLink1(aa.vertex1(), bb.vertex2(), vab, ab2, ran, push_);
        }
        else
        {
            if ( dot(vab, aa.prevDiff()) >= 0 )
            {
                real ab2 = vab.normSqr();
                if ( below(ab2, ran) )
                    meca_.addLongLink1(aa.vertex1(), bb.vertex2(), vab, ab2, ran, push_);
            }
        }
    }
    else if (( &bb < &aa ) & aa.isLast() )
    {
        /*
         Check the projection of aa.vertex2(),
         on the segment represented by 'bb'
         */
        assert_true(abs > aa.seg());
        assert_true(bb.isLast());
        
        Vector dab = bb.pos2() - aa.pos2();
        if ( modulo )
            modulo->fold(dab);
        real ab2 = dab.normSqr();
        
        if ( below(ab2, ran)  &&  dot(dab, bb.diff()) <= 0 )
            meca_.addLongLink1(aa.vertex2(), bb.vertex2(), dab, ab2, ran, push_);
    }
}


/** Older version */
/**
 This is used to check a segment of a fiber against the vertex of a fiber

 The interaction is applied only if the vertex projects 'inside' the segment.
 */
void LocusGrid::checkLL1(BigLocus const& aa, BigLocus const& bb) const
{
    //std::clog << "   LL1 " << aa.segment() << " " << bb.vertex1() << '\n';

    const float ran = aa.rad_ + bb.rad_;

    // get position of bb.vertex1() with respect to segment 'aa'
    real dis2 = INFINITY;
    FiberSegment seg = aa.segment();
    real abs = seg.projectPoint0(bb.pos1(), dis2);
    
    if ( below(dis2, ran) &&  aa.within(abs) )
    {
        /*
         bb.vertex1() projects inside segment 'aa'
         */
        meca_.addSideSlidingLink(seg, abs, bb.vertex1(), ran, push_);
    }
    else if ( abs < 0 )
    {
        if ( aa.isFirst() )
        {
            /*
             Check the projection of aa.vertex1(),
             on the segment represented by 'bb'
             */
            if ( &bb < &aa  &&  bb.isFirst() )
            {
                Vector vab = bb.pos1() - aa.pos1();
                if ( modulo )
                    modulo->fold(vab);
                real ab2 = vab.normSqr();
                if ( below(ab2, ran)  &&  dot(vab, bb.diff()) >= 0 )
                    meca_.addLongLink1(aa.vertex1(), bb.vertex1(), vab, ab2, ran, push_);
            }
        }
        else
        {
            Vector vab = bb.pos1() - aa.pos1();
            if ( modulo )
                modulo->fold(vab);
            /*
             Check the projection to the segment located before 'aa',
             and interact if 'bb.vertex1()' falls on the right side of it
             */
            if ( dot(vab, aa.prevDiff()) >= 0 )
            {
                real ab2 = vab.normSqr();
                if ( below(ab2, ran) )
                    meca_.addLongLink1(aa.vertex1(), bb.vertex1(), vab, ab2, ran, push_);
            }
        }
    }
}


/**
 This is used to check a segment of a fiber against the terminal vertex of a fiber

 The interaction is applied only if the vertex projects 'inside' the segment.
 */
void LocusGrid::checkLL2(BigLocus const& aa, BigLocus const& bb) const
{
    //std::clog << "   LL2 " << aa.segment() << " " << bb.vertex2() << '\n';

    const real ran = aa.rad_ + bb.rad_;
    
    // get position of bb.vertex2() with respect to segment 'aa'
    real dis2 = INFINITY;
    FiberSegment seg = aa.segment();
    real abs = seg.projectPoint0(bb.pos2(), dis2);
    
    if ( aa.within(abs) )
    {
        /*
         bb.vertex2() projects inside segment 'aa'
         */
        if ( below(dis2, ran) )
            meca_.addSideSlidingLink(seg, abs, bb.vertex2(), ran, push_);
    }
    else if ( abs < 0 )
    {
        /*
         Check the projection to the segment located before 'aa',
         and interact if 'bb.vertex1()' falls on the right side of it
         */
        Vector vab = bb.pos2() - aa.pos1();
        
        if ( modulo )
            modulo->fold(vab);
        if ( aa.isFirst() )
        {
            assert_true(bb.isLast());
            real ab2 = vab.normSqr();

            if ( below(ab2, ran)  && dot(vab, bb.diff()) <= 0 )
                meca_.addLongLink1(aa.vertex1(), bb.vertex2(), vab, ab2, ran, push_);
        }
        else
        {
            if ( dot(vab, aa.prevDiff()) >= 0 )
            {
                real ab2 = vab.normSqr();
                if ( below(ab2, ran) )
                    meca_.addLongLink1(aa.vertex1(), bb.vertex2(), vab, ab2, ran, push_);
            }
        }
    }
    else if (( &bb < &aa ) & aa.isLast() )
    {
        /*
         Check the projection of aa.vertex2(),
         on the segment represented by 'bb'
         */
        assert_true(abs > aa.seg());
        assert_true(bb.isLast());
        
        Vector vab = bb.pos2() - aa.pos2();
        if ( modulo )
            modulo->fold(vab);
        real ab2 = vab.normSqr();
        
        if ( below(ab2, ran)  &&  dot(vab, bb.diff()) <= 0 )
            meca_.addLongLink1(aa.vertex2(), bb.vertex2(), vab, ab2, ran, push_);
    }
}


#if ( DIM >= 3 )
/**
 This is used to check two FiberSegment, each representing a segment of a Fiber.
 The segments are tested for intersection in 3D.
 */
void LocusGrid::checkLL(BigLocus const& aa, BigLocus const& bb) const
{
    //std::clog << aa.segment().to_string() << " " << bb.segment().to_string() << "\n";
    assert_true( aa.obj_->tag() == Fiber::TAG );
    assert_true( bb.obj_->tag() == Fiber::TAG );

    const float ran = aa.rad_ + bb.rad_;

    /*
     in 3D, check the shortest distance between two segments, and if close
     enough, use the result to build an interaction.
     The code below is from FiberSegment::shortestDistanceSqr(),
     but having it here inline allows to use the results also for checkLL1()
     The operations could be done in single precision, or 2x2 using SIMD
     */
    real a, b;

    //std::clog << to_string() << " " << seg.to_string() << "\n";
    Vector a1 = aa.pos1();
    Vector b1 = bb.pos1();

    Vector daa = ( aa.pos2() - a1 ) * aa.segInv();
    Vector dbb = ( b1 - bb.pos2() ) * bb.segInv(); // sign inverted on purpose
    Vector off = b1 - a1;
    if ( modulo )
        modulo->fold(off);
    // direction axis of the shortest path is orthogonal to both lines
    Vector axis = cross(daa, dbb);
    real C = dot(daa, dbb);  // cosine of angle between lines

    real iS = dot(axis, axis);  // sine^2
    real D = dot(off, axis);
    real m1 = dot(off, daa);  // = abscissa on segment A of projection of b1
    real m2 = dot(off, dbb);  // = abscissa on segment B of projection of a1
    real dis2;
    
    if ( abs_real(iS) > 128 * REAL_EPSILON )
    {
        // This deals with the general case of non-parallel lines
        iS = 1 / iS;    // 1.0 / sine^2
        a = ( m1 - C * m2 ) * iS;
        b = ( m2 - C * m1 ) * iS;
        dis2 = ( D * D ) * iS;
    }
    else
    {
        // nearly parallel lines, which is unlikely in practice
        const real len1 = aa.seg();
        const real len2 = bb.seg();
        real p1 = m1 - C * len2;
        real p2 = m2 - C * len1;
        // clamp inside segment and use mid-point
        a = 0.5 * ( min_real(len1, max_real(m1, p1)) + max_real(0, min_real(m1, p1)));
        // clamp inside segment and use mid-point
        b = 0.5 * ( min_real(len2, max_real(m2, p2)) + max_real(0, min_real(m2, p2)));
        dis2 = off.normSqrSubtracted(m1);
    }
    
    if ( dis2 < ran*ran )
    {
        dis2 = off.normSqr();
        
        if ( aa.within(a) && bb.within(b) )
        {
            // Since axis is orthogonal to daa, we know the norm of the cross-product:
            Vector leg = cross(daa, axis) * ( copysign(ran, D) * aa.segInv() * sqrt(iS) );
            meca_.addSideSlidingLink3D(aa.interpolation(a), leg, bb.interpolation(b), daa, push_);
        }
        
        /* If the shortest distance between the lines is greater than 'ran', then
          the vertices associated with this segment will also be too far to interact */

        //std::clog << "   LL " << aa.segment() << " " << bb.segment() << '\n';
        // the distance between bb.pos1() and segment `aa` is dis2-m1*m1:
        if ( below(dis2-m1*m1, ran) && aa.within(m1) )
        {
            Vector leg = cross(daa, off).normalized( ran * aa.segInv() );
            meca_.addSideSlidingLink3D(aa.interpolation(m1), leg, bb.vertex1(), daa, push_);
        }
        else if ( m1 < 0 )
            checkLLP1(aa, bb, ran, m1, off);

        if ( bb.isLast() )
        {
            // consider bb.vertex2() against segment 'aa'
            Vector a1b2 = bb.pos2() - a1;
            real mm = dot(daa, a1b2); // projection of b2 on segment 'aa'
            if ( aa.within(mm) )
            {
                // the distance between b2 and segment `aa` is a1b2.normSqr()-mm*mm:
                if ( below(a1b2.normSqr()-mm*mm, ran) )
                {
                    Vector leg = cross(daa, a1b2).normalized( ran * aa.segInv() );
                    meca_.addSideSlidingLink3D(aa.interpolation(mm), leg, bb.vertex2(), daa, push_);
                }
            }
            else
                checkLLP2(aa, bb, ran, mm, a1b2);
        }

        // the distance between aa.pos1() and segment `bb` is dis2-m2*m2:
        if ( below(dis2-m2*m2, ran) && bb.within(m2) )
        {
            Vector leg = cross(dbb, off).normalized( ran * bb.segInv() ); // two minus sign cancel out
            meca_.addSideSlidingLink3D(bb.interpolation(m2), leg, aa.vertex1(), -dbb, push_);
        }
        else if ( m2 < 0 )
            checkLLP1(bb, aa, ran, m2, -off);

        if ( aa.isLast() )
        {
            // consider aa.vertex2() against segment 'bb'
            Vector a2b1 = b1 - aa.pos2(); // sign inverted on purpose
            real mm = dot(dbb, a2b1); // projection of a2 on segment 'bb'
            if ( bb.within(mm) )
            {
                if ( below(a2b1.normSqr()-mm*mm, ran) )
                {
                    Vector leg = cross(dbb, a2b1).normalized( ran * bb.segInv() ); // two minus sign cancel out
                    meca_.addSideSlidingLink3D(bb.interpolation(mm), leg, aa.vertex2(), -dbb, push_);
                }
            }
            else
                checkLLP2(bb, aa, ran, mm, -a2b1);
        }
    }
}

#elif ( DIM == 2 )

/**
 This is used to check two FiberSegment, each representing a segment of a Fiber.
 This is the 2D implementation.
 */
void LocusGrid::checkLL(BigLocus const& aa, BigLocus const& bb) const
{
    assert_true( aa.obj_->tag() == Fiber::TAG );
    assert_true( bb.obj_->tag() == Fiber::TAG );

    const float ran = aa.rad_ + bb.rad_;
    // get position of bb.vertex1() with respect to segment 'aa'
    Vector a1 = aa.pos1();
    Vector b1 = bb.pos1();

    Vector daa = ( aa.pos2() - a1 ) * aa.segInv();
    Vector dbb = ( b1 - bb.pos2() ) * bb.segInv(); // sign inverted on purpose
    Vector off = b1 - a1;
    if ( modulo )
        modulo->fold(off);
    real m1 = dot(off, daa);
    real m2 = dot(off, dbb);
    real dis2 = off.normSqr();

    //std::clog << "   LL " << aa.segment() << " " << bb.segment() << '\n';
    // the distance between bb.pos1() and segment `aa` is dis2-m1*m1:
    if ( below(dis2-m1*m1, ran) && aa.within(m1) )
    {
        real leg = std::copysign(ran*aa.segInv(), cross(daa, off));
        meca_.addSideSlidingLink2D(aa.interpolation(m1), leg, bb.vertex1(), daa, push_);
    }
    else if ( m1 < 0 )
        checkLLP1(aa, bb, ran, m1, off);
    
    if ( bb.isLast() )
    {
        // consider bb.vertex2() against segment 'aa'
        Vector a1b2 = bb.pos2() - a1;
        real mm = dot(daa, a1b2); // projection of b2 on segment 'aa'
        if ( aa.within(mm) )
        {
            // the distance between b2 and segment `aa` is a1b2.normSqr()-mm*mm:
            if ( below(a1b2.normSqr()-mm*mm, ran) )
            {
                real leg = std::copysign(ran*aa.segInv(), cross(daa, a1b2));
                meca_.addSideSlidingLink2D(aa.interpolation(mm), leg, bb.vertex2(), daa, push_);
            }
        }
        else
            checkLLP2(aa, bb, ran, mm, a1b2);
    }
    
    // the distance between aa.pos1() and segment `bb` is dis2-m2*m2:
    if ( below(dis2-m2*m2, ran) && bb.within(m2) )
    {
        real leg = std::copysign(ran*bb.segInv(), cross(dbb, off)); // two minus sign cancel out
        meca_.addSideSlidingLink2D(bb.interpolation(m2), leg, aa.vertex1(), -dbb, push_);
    }
    else if ( m2 < 0 )
        checkLLP1(bb, aa, ran, m2, -off);

    if ( aa.isLast() )
    {
        // consider aa.vertex2() against segment 'bb'
        Vector a2b1 = b1 - aa.pos2(); // sign inverted on purpose
        real mm = dot(dbb, a2b1); // projection of a2 on segment 'bb'
        if ( bb.within(mm) )
        {
            if ( below(a2b1.normSqr()-mm*mm, ran) )
            {
                real leg = std::copysign(ran*bb.segInv(), cross(dbb, a2b1)); // two minus sign cancel out
                meca_.addSideSlidingLink2D(bb.interpolation(mm), leg, aa.vertex2(), -dbb, push_);
            }
        }
        else
            checkLLP2(bb, aa, ran, mm, -a2b1);
    }
}

#elif ( DIM == 1 )

/**
 This is used to check two FiberSegment, each representing a segment of a Fiber.
 This is the 2D implementation.
 */
void LocusGrid::checkLL(BigLocus const&, BigLocus const&) const
{
    throw Exception("Steric is meaningless in 1D");
}

#else

/**
 This is an older and simpler version, but some calculations are duplicated in 3D:
 checkLL1() will project point1 of segment 1 on segment 2,
 but this calculation is already done in shortestDistanceSqr()
 */
void LocusGrid::checkLL(BigLocus const& aa, BigLocus const& bb) const
{
    assert_true( aa.obj_->tag() == Fiber::TAG );
    assert_true( bb.obj_->tag() == Fiber::TAG );

#if ( DIM >= 3 )
    
    const real ran = aa.rad_ + bb.rad_;
    
    FiberSegment as = aa.segment();
    FiberSegment bs = bb.segment();

    /* in 3D, check the shortest distance between two segments, and if close
     enough, use the result to build an interaction */
    real a, b;
    real dis2 = as.shortestDistanceSqr(bs, a, b);
    
    if ( dis2 < ran*ran )
    {
        if ( as.within(a) && bs.within(b) )
            meca_.addSideSlidingLink(as, a, Interpolation(bs, b), ran, push_);
    }
    else
    {
        /* If the shortest distance between the lines is greater than 'ran', then
          the vertices associated with this segment will also be too far to interact */
        return;
    }
#endif

    //std::clog << "   LL " << aa.segment() << " " << bb.segment() << '\n';
    checkLL1(aa, bb);
    
    if ( aa.isLast() )
        checkLL2(bb, aa);
    
    checkLL1(bb, aa);
    
    if ( bb.isLast() )
        checkLL2(aa, bb);
}

#endif

//------------------------------------------------------------------------------
#pragma mark - Selections of pairs excluded from Sterics

/*
 In general, these test will only exclude relatively rare pairs from interacting,
 and thus are less stringent than BigVector::near(): they should be tested after.
 */

/// excluding two spheres when they are from the same Solid
static inline bool not_adjacentPP(BigPoint const& a, BigPoint const& b)
{
    return a.obj_ != b.obj_;
}


/// excluding Fiber and Solid from the same Aster
static inline bool not_adjacentPL(BigPoint const& a, BigLocus const& b)
{
    //a->mec_->Buddy::print(std::clog);
    //b->obj_->Buddy::print(std::clog);
    return ! b.obj_->isBuddy(a.obj_);
}


/// excluding segments that are adjacent on the same fiber, or protofilaments from Tubule
static inline bool not_adjacentLL(BigLocus const& a, BigLocus const& b)
{
#if FIBER_HAS_FAMILY
    Fiber const* fibA = static_cast<Fiber const*>(a.obj_);
    Fiber const* fibB = static_cast<Fiber const*>(b.obj_);
    return (( fibA->family_ != fibB->family_ )
            || (( a.vix_ > 1 + b.vix_ ) | ( b.vix_ > 1 + a.vix_ )));
#else
    return (( a.obj_ != b.obj_ )
            || (( a.vix_ > 1 + b.vix_ ) | ( b.vix_ > 1 + a.vix_ )));
#endif
    // we cannot use abs() above because `vix_` is unsigned
}

//------------------------------------------------------------------------------
#pragma mark - Check all possible object pairs from two Cells

/**
 This will consider once all pairs of objects from the given lists
 */
void LocusGrid::setSterics0(BigLocusList const& list) const
{
    BigLocus const* mid = list.middle();
    //printf("LG:setSterics0  %lu:%lu\n", mid-list.begin(), list.end()-mid);
    
    for ( BigLocus const* ii = list.begin(); ii < mid; ++ii )
    {
        for ( BigLocus const* jj = ii+1; jj < mid; ++jj )
            if ( not_adjacentLL(*ii, *jj) )
                checkLL(*ii, *jj);
        
        for ( BigPoint const* kk = mid; kk < list.end(); ++kk )
            if ( not_adjacentPL(*kk, *ii) )
                checkPL(*kk, *ii);
    }

    for ( BigPoint const* ii = mid; ii < list.end(); ++ii )
    {
        for ( BigPoint const* jj = ii+1; jj < list.end(); ++jj )
            if ( not_adjacentPP(*ii, *jj) )
                checkPP(*ii, *jj);
    }
}


/**
 This will consider once all pairs of objects from the given lists,
 assuming that the list are different and no object is repeated
 */
void LocusGrid::setSterics0(BigLocusList const& list1,
                            BigLocusList const& list2) const
{
    assert_true( &list1 != &list2 );
    BigLocus const* mid1 = list1.middle();
    BigLocus const* mid2 = list2.middle();

    for ( BigLocus const* ii = list1.begin(); ii < mid1; ++ii )
    {
        for ( BigLocus const* jj = list2.begin(); jj < mid2; ++jj )
            if ( not_adjacentLL(*ii, *jj)  )
                checkLL(*ii, *jj);

        for ( BigPoint const* kk = mid2; kk < list2.end(); ++kk )
            if ( not_adjacentPL(*kk, *ii) )
                checkPL(*kk, *ii);
    }

    for ( BigPoint const* ii = mid1; ii < list1.end(); ++ii )
    {
        for ( BigLocus const* jj = list2.begin(); jj < mid2; ++jj )
            if ( not_adjacentPL(*ii, *jj) )
                checkPL(*ii, *jj);

        for ( BigPoint const* jj = mid2; jj < list2.end(); ++jj )
            if ( not_adjacentPP(*ii, *jj) )
                checkPP(*ii, *jj);
    }
}


/**
 This will consider once all pairs of objects from the given lists.
 Compared to `setSterics0()`, this performs additional tests to exclude
 objects that are too far apart to interact, based on BigVector::near()
 */
void LocusGrid::setStericsT(BigLocusList const& list) const
{
    BigLocus const* mid = list.middle();

    for ( BigLocus const* ii = list.begin(); ii < mid; ++ii )
    {
        const BigVector pos = ii->pos_;
        
        for ( BigLocus const* jj = ii+1; jj < mid; ++jj )
            if ( pos.near(jj->pos_) && not_adjacentLL(*ii, *jj) )
                checkLL(*ii, *jj);
        
        for ( BigPoint const* jj = mid; jj < list.end(); ++jj )
            if ( pos.near(jj->pos_) && not_adjacentPL(*jj, *ii) )
                checkPL(*jj, *ii);
    }

    for ( BigPoint const* ii = mid; ii < list.end(); ++ii )
    {
        const BigVector pos = ii->pos_;

        for ( BigPoint const* jj = ii+1; jj < list.end(); ++jj )
            if ( pos.near(jj->pos_) && not_adjacentPP(*ii, *jj) )
                checkPP(*ii, *jj);
    }
}


/**
 This will consider once all pairs of objects from the given lists,
 assuming that the list are different and no object is repeated.

 Compared to `setSterics0()`, this performs additional tests to exclude
 objects that are too far apart to interact, based on BigVector::near()
*/
void LocusGrid::setStericsT(BigLocusList const& list1,
                            BigLocusList const& list2) const
{
    assert_true( &list1 != &list2 );
    BigLocus const* mid1 = list1.middle();
    BigLocus const* mid2 = list2.middle();
    
    //std::clog << std::setw(4) << list1.size() << " vs " << std::setw(4) << list2.size() << "\n";
    /*
     The tests pos.near(jj->pos_) can be calculated using SIMD instructions
     */
    for ( BigLocus const* ii = list1.begin(); ii < mid1; ++ii )
    {
        BigVector pos = ii->pos_;
        if ( modulo )
            modulo->fold_float(pos, list2[0].pos_);
        for ( BigLocus const* jj = list2.begin(); jj < mid2; ++jj )
            if ( pos.near(jj->pos_) && not_adjacentLL(*ii, *jj)  )
                checkLL(*ii, *jj);

        for ( BigPoint const* jj = mid2; jj < list2.end(); ++jj )
            if ( pos.near(jj->pos_) && not_adjacentPL(*jj, *ii) )
                checkPL(*jj, *ii);
    }

    for ( BigPoint const* ii = mid1; ii < list1.end(); ++ii )
    {
        BigVector pos = ii->pos_;
        if ( modulo )
            modulo->fold_float(pos, list2[0].pos_);
        for ( BigLocus const* jj = list2.begin(); jj < mid2; ++jj )
            if ( pos.near(jj->pos_) && not_adjacentPL(*ii, *jj) )
                checkPL(*ii, *jj);

        for ( BigPoint const* jj = mid2; jj < list2.end(); ++jj )
            if ( pos.near(jj->pos_) && not_adjacentPP(*ii, *jj) )
                checkPP(*ii, *jj);
    }
}

#if ( DIM == 3 ) && USE_SIMD

/// BitField must be a 64 bit integer type
typedef uint64_t BitField;

/**
 Evaluate 4 pos.near(jj->pos_) using SIMD instructions
 @return a 4-bit integer where each bit represents the result of one test
 */
inline BitField four_near_bits(vec4f const& xyzr, BigLocus const* src)
{
    vec4f tt = sub4f(xyzr, loadu4f(src[0].pos_));
    vec4f yy = sub4f(xyzr, loadu4f(src[1].pos_));
    vec4f uu = sub4f(xyzr, loadu4f(src[2].pos_));
    vec4f rr = sub4f(xyzr, loadu4f(src[3].pos_));
    // transpose 4x4 data matrix:
    vec4f xx = unpacklo4f(tt, yy);
    tt = unpackhi4f(tt, yy);
    vec4f zz = unpacklo4f(uu, rr);
    rr = unpackhi4f(uu, rr);
    yy = movelh4f(xx, zz);
    xx = movehl4f(zz, xx);
    uu = fmadd4f(yy, yy, mul4f(xx, xx)); // x*x + y*y
    zz = movelh4f(tt, rr);  // zz, rr
    rr = movehl4f(rr, tt);
    // calculate test:
    tt = fnmadd4f(zz, zz, mul4f(rr, rr)); // r*r - z*z
    return lower_mask4f(uu, tt);  // x*x + y*y < r*r - z*z
}

inline BitField two_near_bits(vec4f const& xyzr, BigLocus const* src)
{
    vec4f tt = sub4f(xyzr, loadu4f(src[0].pos_));
    vec4f yy = sub4f(xyzr, loadu4f(src[1].pos_));
    // transpose 4x4 data matrix:
    vec4f xx = unpacklo4f(tt, yy);
    tt = unpackhi4f(tt, yy);
    yy = movelh4f(xx, setzero4f());
    xx = movehl4f(setzero4f(), xx);
    xx = fmadd4f(yy, yy, mul4f(xx, xx)); // x*x + y*y
    yy = movelh4f(tt, setzero4f());  // zz, rr
    vec4f rr = movehl4f(setzero4f(), tt);
    // calculate test:
    tt = fnmadd4f(yy, yy, mul4f(rr, rr)); // r*r - z*z
    return lower_mask4f(xx, tt);  // x*x + y*y < r*r - z*z
}

inline BitField one_near_bit(vec4f const& xyzr, BigLocus const* src)
{
    vec4f xx = sub4f(xyzr, loadu4f(src[0].pos_));
    xx = mul4f(xx, xx);
    return xx[0] + xx[1] < xx[3] - xx[2];  // x*x + y*y < r*r - z*z
}

/**
Evaluate `cnt` pos.near(jj->pos_) using SIMD instructions
@return a bitfield representing the result of all tests
*/
BitField compute_near_bits(vec4f const& xyzr, BigLocus const* start, size_t cnt)
{
    BitField res = 0;
    unsigned shift = 0;
    BigLocus const* ptr = start;
    BigLocus const* end;
#if 0
    /*
     Unrolling can help since all four_near_bits() are independent,
     and the calculations can be executed out-of-order efficiently,
     but this is effective only for list size > 32
     */
    end = start + ( cnt & ~15UL );
    while ( ptr < end )
    {
        BitField t = four_near_bits(xyzr, ptr) << shift;
        BitField u = four_near_bits(xyzr, ptr+4) << (shift+4);
        BitField v = four_near_bits(xyzr, ptr+8) << (shift+8);
        BitField w = four_near_bits(xyzr, ptr+12) << (shift+12);
        res |= (( t | u ) | ( v | w ));
        shift += 16;
        ptr += 16;
    }
#endif
    end = start + ( cnt & ~7UL );
    while ( ptr < end )
    {
        BitField t = four_near_bits(xyzr, ptr) << shift;
        BitField u = four_near_bits(xyzr, ptr+4) << (shift+4);
        res |= ( t | u );
        shift += 8;
        ptr += 8;
    }
    
    end = start + ( cnt & ~3UL );
    if ( ptr < end )
    {
        BitField t = four_near_bits(xyzr, ptr);
#if 0
        BitField a = two_near_bits(xyzr, ptr) + (two_near_bits(xyzr, ptr+2)<<2);
        if ( t != a ) printf("[%lX %lX]", t, a);
#endif
        res |= t << shift;
        shift += 4;
        ptr += 4;
    }
    
    end = start + ( cnt & ~1UL );
    if ( ptr < end )
    {
        BitField t = two_near_bits(xyzr, ptr);
        res |= t << shift;
        shift += 2;
        ptr += 2;
    }

    end = start + cnt;
    while ( ptr < end )
    {
        //printf("at %lu remaining %lu\n", ptr-start, end-ptr);
#if 0
        BitField t = four_near_bits(xyzr, ptr); // we may load invalid data
        unsigned i = end - ptr;
        assert_true( i < 8 );
        // a mask to clear the bits past the end:
        BitField k = ~( ~0LU << i );
        res |= ( t & k ) << shift;
        shift += 4;
        ptr += 4;
#else
        BitField t = one_near_bit(xyzr, ptr); // we may load invalid data
        res |= t << shift;
        ++shift;
        ++ptr;
#endif
    }

    //printf(" near_bits %2lu : %lu\n", cnt, res);
    return res;
}


/**
 Set bitL corresponding to BigLocus in first part of list,
 and bitP corresponding to BigPoints in second part of list
*/
void compute_near_bits(BitField& bitL, BitField& bitP, vec4f const& xyzr, BigLocusList const& list, size_t start)
{
    size_t cnt = std::min(list.size()-start, 64UL);
    BitField bits = compute_near_bits(xyzr, list.begin()+start, cnt);
    size_t nloc = list.num_locus();
    size_t shift = nloc - std::min(start, nloc);
    if ( shift < 64 )
    {
        BitField mask = ~0UL << shift;
        bitL = bits & ~mask;
        bitP = bits & mask;
    } else {
        bitL = bits;
        bitP = 0;
    }
}


void LocusGrid::setStericsX(BigLocusList const& list) const
{
    //printf(" stericsX: %2lu+%2lu (%lu)\n", list.num_locus(), list.num_points(), list.capacity());
    BigLocus const* mid = list.middle();
    BitField bitP, bitL;
    
    for ( BigLocus const* ii = list.begin(); ii < mid; ++ii )
    {
        BigVector pos = ii->pos_;
        /* The radius is negated, such that it gets added by compute_near_bits() */
        vec4f xyzr { pos.xx, pos.yy, pos.zz, -pos.rr };
        /* In most situations, the list size would be < 64 and one round would
         be sufficient, but in all generality we must handle larger list size */
        for ( size_t offset = 0; offset < list.size(); offset += 64 )
        {
            size_t start = offset + ( ii - list.begin() );
            compute_near_bits(bitL, bitP, xyzr, list, start);
            //printf(" LL %i : ", __builtin_popcount(bitL));
            BigLocus const* blp = list.begin() + start;
            while ( bitL )
            {
                //printf(" %llX", bitL);
                int b = __builtin_ctzl(bitL);
                if ( not_adjacentLL(*ii, blp[b]) )
                    checkLL(*ii, blp[b]);
                // flip bit in position 'b':
                bitL ^= 1UL << b;
            }
            //printf(" LP ");
            while ( bitP )
            {
                //printf(" %llX", bitP);
                int b = __builtin_ctzl(bitP);
                if ( not_adjacentPL(blp[b], *ii) )
                    checkPL(blp[b], *ii);
                // flip bit in position 'b':
                bitP ^= 1UL << b;
            }
            //printf("\n");
        }
    }
    
    for ( BigPoint const* ii = mid; ii < list.end(); ++ii )
    {
        BigVector pos = ii->pos_;
        /* The radius is negated, such that it gets added by compute_near_bits() */
        vec4f xyzr { pos.xx, pos.yy, pos.zz, -pos.rr };
        /* In most situations, the list size would be < 64 and one round would
         be sufficient, but in all generality we must handle larger list size */
        for ( size_t offset = 0; offset < list.size(); offset += 64 )
        {
            size_t start = offset + ( ii - list.begin() );
            size_t cnt = std::min(list.size()-start, 64UL);
            bitP = compute_near_bits(xyzr, list.begin()+start, cnt);
            //printf(" PP ");
            BigLocus const* blp = list.begin() + start;
            while ( bitP )
            {
                //printf(" %llX", bitP);
                int b = __builtin_ctzl(bitP);
                if ( not_adjacentPP(*ii, blp[b]) )
                    checkPP(*ii, blp[b]);
                bitP ^= 1UL << b;
            }
            //printf("\n");
        }
    }
}

/**
 This will consider once all pairs of objects from the given lists,
 assuming that the list are different and no object is repeated.

 Compared to `setSterics0()`, this performs additional tests to exclude
 objects that are too far apart to interact, based on BigVector::near()
 
 The same approach can be used for periodic boundary conditions, if:
 - BigVector should be folded to their cannonical representation
 - distance should be calculated adding an offset, for cells that
   are accross a periodic boundary. Note that this offset is defined per
   cell pairs, and not per object pair: just need to update `xyzr` below.
 .
 
 This code relies on '__builtin_ctzl(x)' which gives the index of the first non-zero bit:
 Returns the number of trailing 0-bits in x, starting at the least significant bit position.
 If x is 0, the result is undefined.
*/
void LocusGrid::setStericsX(BigLocusList const& list1,
                            BigLocusList const& list2) const
{
    assert_true( &list1 != &list2 );
#if 0
    {
        index_t l1 = list1.num_locus(), p1 = list1.num_points();
        index_t l2 = list2.num_locus(), p2 = list2.num_points();
        printf(" stericsX: %2lu+%2lu  :  %2lu+%2lu\n", l1, p1, l2, p2);
    }
#endif
    BigLocus const* mid1 = list1.middle();
    BitField bitP, bitL;
    
    for ( BigLocus const* ii = list1.begin(); ii < mid1; ++ii )
    {
        BigVector pos = ii->pos_;
        if ( modulo )
            modulo->fold_float(pos, list2[0].pos_);
        /* The radius is negated, such that it gets added by compute_near_bits() */
        vec4f xyzr { pos.xx, pos.yy, pos.zz, -pos.rr };
        /* In most situations, the list size would be < 64 and one round would
         be sufficient, but in all generality we must handle larger list size */
        for ( size_t offset = 0; offset < list2.size(); offset += 64 )
        {
            compute_near_bits(bitL, bitP, xyzr, list2, offset);
            //printf(" LL");
            BigLocus const* blp = list2.begin() + offset;
#if 0
            // verify all tests:
            for ( size_t b = 0; b < std::min(64UL, list2.size()-offset); ++b )
            {
                BigVector vec = (ii+b+offset)->pos_;
                if ( modulo ) modulo->fold_float(vec, pos);
                bool n = pos.near(vec);
                bool p = (bitL+bitP) & ( 1UL << b );
                float d = square(pos.XX-vec.XX) + square(pos.YY-vec.YY) + square(pos.ZZ-vec.ZZ);
                // SIMD and scalar results can differ, near the edges (x+y ~ 0) :
                if ( n != p )
                {
                    printf("%8.2f |  %8.2f %8.2f | ", list2[0].pos_.XX, pos.XX, (ii+offset+b)->pos_.XX);
                    printf("!near %2lu+%2lu: %u%u (%f)\n", offset, b, n, p, d);
                }
            }
#endif
            while ( bitL )
            {
                //printf(" %lX", bitL);
                int b = __builtin_ctzl(bitL);
                if ( not_adjacentLL(*ii, blp[b]) )
                    checkLL(*ii, blp[b]);
                // flip bit in position 'b':
                bitL ^= 1UL << b;
            }
            //printf(" LP ");
            while ( bitP )
            {
                //printf(" %lX", bitP);
                int b = __builtin_ctzl(bitP);
                if ( not_adjacentPL(blp[b], *ii) )
                    checkPL(blp[b], *ii);
                // flip bit in position 'b':
                bitP ^= 1UL << b;
            }
            //printf("\n");
        }
    }
    
    for ( BigPoint const* ii = mid1; ii < list1.end(); ++ii )
    {
        BigVector pos = ii->pos_;
        if ( modulo )
            modulo->fold_float(pos, list2[0].pos_);
        /* The radius is negated, such that it gets added by compute_near_bits() */
        vec4f xyzr { pos.xx, pos.yy, pos.zz, -pos.rr };
        /* In most situations, the list size would be < 64 and one round would
         be sufficient, but in all generality we must handle larger list size */
        for ( size_t offset = 0; offset < list2.size(); offset += 64 )
        {
            compute_near_bits(bitL, bitP, xyzr, list2, offset);
            //printf(" PL ");
            BigLocus const* blp = list2.begin() + offset;
            while ( bitL )
            {
                //printf(" %lX", bitL);
                int b = __builtin_ctzl(bitL);
                if ( not_adjacentPL(*ii, blp[b]) )
                    checkPL(*ii, blp[b]);
                bitL ^= 1UL << b;
            }
            //printf(" PP ");
            while ( bitP )
            {
                //printf(" %lX", bitP);
                int b = __builtin_ctzl(bitP);
                if ( not_adjacentPP(*ii, blp[b]) )
                    checkPP(*ii, blp[b]);
                bitP ^= 1UL << b;
            }
            //printf("\n");
        }
    }
}
#endif

//------------------------------------------------------------------------------
#pragma mark - Check all pairs of Cells

#if ( MAX_STERIC_PANES == 1 )

/**
 Check interactions between objects contained in the grid:
 Scan all cells to examine all object pairs (ii, jj) only once.
 This version can handle periodic boundary conditions
 */
void LocusGrid::setSterics0() const
{
    for ( index_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int const* region;
        int nR = pGrid.getRegion(region, inx);
        
        BigLocusList& base = cell_list(inx);
        setSterics0(base);
        
        for ( int reg = 0; reg < nR; ++reg )
            setSterics0(base, cell_list(inx+region[reg]));
    }
}


/** This calls setStericsT() */
void LocusGrid::setSterics() const
{
    for ( index_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int const* region;
        int nR = pGrid.getRegion(region, inx);
        
        BigLocusList& base = cell_list(inx);
        if ( base.size() > 0 )
        {
            setSterics0(base);
#if ( DIM == 3 ) && USE_SIMD
            for ( int reg = 0; reg < nR; ++reg )
            {
                BigLocusList& side = cell_list(inx+region[reg]);
                if ( base.size() < side.size() )
                    setStericsX(base, side);
                else
                {
                    if ( side.size() > 0 )
                        setStericsX(side, base);
                }
            }
#else
            for ( int reg = 0; reg < nR; ++reg )
            {
                BigLocusList& side = cell_list(inx+region[reg]);
                if ( side.size() > 0 )
                    setSterics0(base, side);
            }
#endif
        }
    }
}


void LocusGrid::setStericsT() const
{
    //std::clog << "----" << '\n';
    if ( pGrid.isPeriodic() )
        setSterics0();
    else
        setStericsT();
}

#else

/**
 Check interactions between objects contained in the pane `pan`:
 Scan all cells to examine all object pairs (ii, jj) only once.
 This version can handle periodic boundary conditions
 */
void LocusGrid::setSterics0(index_t pan) const
{
    for ( index_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int const* region;
        int nR = pGrid.getRegion(region, inx);
       
        BigLocusList& base = cell_list(inx, pan);
        setSterics0(base);
        
        for ( int reg = 0; reg < nR; ++reg )
            setSterics0(base, cell_list(inx+region[reg], pan));
    }
}


void LocusGrid::setStericsT(index_t pan) const
{
    assert_false( pGrid.isPeriodic() );
    for ( index_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
         int const* region;
         int nR = pGrid.getRegion(region, inx);
        
         BigLocusList& base = cell_list(inx, pan);
         setStericsT(base);
         
         for ( int reg = 0; reg < nR; ++reg )
             setStericsT(base, cell_list(inx+region[reg], pan));
    }
}


/**
 Check interactions between the FatPoints contained in Panes `pan` and `bim`,
 where ( pan1 != pan2 )
 */
void LocusGrid::setSterics0(index_t pan, index_t bim) const
{
    assert_true(pan != bim);
    for ( index_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int const* region;
        int nR = pGrid.getRegion(region, inx);
        
        BigLocusList& base1 = cell_list(inx, pan);
        BigLocusList& base2 = cell_list(inx, bim);

        setSterics0(base1, base2);

        for ( int reg = 0; reg < nR; ++reg )
        {
            setSterics0(base1, cell_list(inx+region[reg], bim));
            setSterics0(base2, cell_list(inx+region[reg], pan));
        }
    }
}


void LocusGrid::setStericsT(index_t pan, index_t bim) const
{
    assert_false( pGrid.isPeriodic() );
    assert_true(pan != bim);
    for ( index_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int const* region;
        int nR = pGrid.getRegion(region, inx);
        
        BigLocusList& base1 = cell_list(inx, pan);
        BigLocusList& base2 = cell_list(inx, bim);

        setStericsT(base1, base2);

        for ( int reg = 1; reg < nR; ++reg )
        {
            setStericsT(base1, cell_list(inx+region[reg], bim));
            setStericsT(base2, cell_list(inx+region[reg], pan));
        }
    }
}


void LocusGrid::setSterics(index_t pan) const
{
    if ( pGrid.isPeriodic() )
        setSterics0(pan);
    else
        setStericsT(pan);
}


void LocusGrid::setSterics(index_t pan1, index_t pan2) const
{
    if ( pGrid.isPeriodic() )
        setSterics0(pan1, pan2);
    else
        setStericsT(pan1, pan2);
}

#endif


//------------------------------------------------------------------------------
#pragma mark - Display

#ifdef DISPLAY

#include "gym_view.h"
#include "gym_draw.h"
#include "gym_cap.h"

void drawBoundaries(Map<DIM> const&, float);

void LocusGrid::drawGrid() const
{
#if ( DIM <= 3 )
    gym::ref_view();
    gym::disableLighting();
    gym::color(1,0,0);
    drawBoundaries(pGrid, 0.5f);
#endif
}
#endif

