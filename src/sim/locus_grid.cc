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
#include "meca.h"
#include "simd.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------

size_t LocusGrid::setGrid(Space const* spc, real min_width)
{
    assert_true( min_width > 0 );
    Vector inf, sup;
    spc->boundaries(inf, sup);
    
    size_t cnt[3] = { 1, 1, 1 };
    for ( size_t d = 0; d < DIM; ++d )
    {
        // minimum number of cells in dimension 'd'
        int n = std::floor(( sup[d] - inf[d] ) / min_width);
        
        if ( n < 0 )
            throw InvalidParameter("invalid space:boundaries");
        
#if GRID_HAS_PERIODIC
        assert_true( modulo == spc->getModulo() );
        if ( modulo  &&  modulo->isPeriodic(d) )
        {
            assert_small( modulo->period_[d] - sup[d] + inf[d] );
            // adjust the grid to match the edges
            cnt[d] = std::max(1, n);
            pGrid.setPeriodic(d, true);
        }
        else
#endif
        {
            //extend in any dimension that is not periodic, adjusting cell size to min_width
            cnt[d] = n + 1;
            real w = cnt[d] * 0.5 * min_width;
            real m = inf[d] + sup[d];
            inf[d] = m - w;
            sup[d] = m + w;
        }
    }
    
    pGrid.setDimensions(inf, sup, cnt);
    return pGrid.nbCells();
}


void LocusGrid::createCells()
{
    pGrid.createCells();
    
    //Create side regions suitable for pairwise interactions:
    pGrid.createSideRegions(1);
    
    //report the grid size used
    pGrid.printSummary(Cytosim::log, "   LocusGrid");
}


void LocusGrid::delimit() const
{
    for ( size_t i = 0; i < pGrid.nbCells(); ++i )
        pGrid[i].delimit();
}


size_t LocusGrid::capacity() const
{
    size_t res = 0;
    for ( size_t i = 0; i < pGrid.nbCells(); ++i )
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
    const real ran = aa.rad_ + bb.rad_;

#if GRID_HAS_PERIODIC
    if ( modulo )
        modulo->fold(vab);
#endif
    real ab2 = vab.normSqr();
    
    if ( below(ab2, ran) )
        meca.addLongLink1(aa.vertex1(), bb.vertex1(), vab, ab2, ran, push);
}


/**
 This is used to check a segment of a fiber against a spherical object:
 Solid/Bead/Sphere or Fiber-tip.
 
 The force is applied if the objects are closer than the maximum of the two range + radius,
 and if the center of the sphere projects inside the segment.
 */
void LocusGrid::checkPL(BigPoint const& aa, BigLocus const& bb) const
{
    //std::clog << "   PL- " << aa.mec_ << " " << bb.obj_ << '\n';
    assert_true( aa.obj_->tag() != Fiber::TAG );
    assert_true( bb.obj_->tag() == Fiber::TAG );

    const real ran = aa.rad_ + bb.rad_;
    
    // determine projection of `aa` on segment `bb`:
    real ab2 = INFINITY;
    real abs = bb.segment().projectPoint0(aa.cen(), ab2);
    
    if ( 0 <= abs )
    {
        if ( abs <= bb.len() )
        {
            // the point projects inside the segment
            if ( below(ab2, ran) )
                meca.addSideSlidingLink(bb.segment(), abs, aa.vertex1(), ran, push);
        }
        else
        {
            // the point projects right past the segment, and we check the fiber tip
            if ( bb.isLast() )
            {
                Vector vab = bb.pos2() - aa.cen();
#if GRID_HAS_PERIODIC
                if ( modulo )
                    modulo->fold(vab);
#endif
                ab2 = vab.normSqr();
                if ( below(ab2, ran) )
                    meca.addLongLink1(aa.vertex1(), bb.vertex2(), vab, ab2, ran, push);
            }
        }
    }
    else
    {
        // here abs < 0, and thus `bb` projects left past the segment
        Vector vab = bb.pos1() - aa.cen();
#if GRID_HAS_PERIODIC
        if ( modulo )
            modulo->fold(vab);
#endif
        ab2 = vab.normSqr();
        /*
         This code handles the interactions will all the joints of a fiber:
         interact with the node only if this projects on the previous segment
         or if this is the terminal point of a fiber.
         */
        if ( below(ab2, ran) && ( bb.isFirst() || dot(vab, bb.prevDiff()) <= 0 ))
            meca.addLongLink1(aa.vertex1(), bb.vertex1(), vab, ab2, ran, push);
    }
}


/**
 This is used to check a segment of a fiber against another segment of fiber,
 not including the terminal vertex of fibers.
 
 The interaction is applied only if the vertex projects 'inside' the segment.
 */
void LocusGrid::checkLL1(BigLocus const& aa, BigLocus const& bb) const
{
    //std::clog << "   LL1 " << aa.obj_ << " " << bb.vertex1() << '\n';
    assert_true( aa.obj_->tag() == Fiber::TAG );
    assert_true( bb.obj_->tag() == Fiber::TAG );

    const real ran = aa.rad_ + bb.rad_;
    
    // get position of bb.vertex1() with respect to segment 'aa'
    real dis2 = INFINITY;
    FiberSegment seg = aa.segment();
    real abs = seg.projectPoint0(bb.pos1(), dis2);
    
    if ( below(dis2, ran) &&  ((0 <= abs) & (abs <= aa.len())) )
    {
        /*
         bb.vertex1() projects inside segment 'aa'
         */
        meca.addSideSlidingLink(seg, abs, bb.vertex1(), ran, push);
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
                
#if GRID_HAS_PERIODIC
                if ( modulo )
                    modulo->fold(vab);
#endif
                real ab2 = vab.normSqr();
                if ( below(ab2, ran)  &&  dot(vab, bb.diff()) >= 0 )
                    meca.addLongLink1(aa.vertex1(), bb.vertex1(), vab, ab2, ran, push);
            }
        }
        else
        {
            /*
             Check the projection to the segment located before 'aa',
             and interact if 'bb.vertex1()' falls on the right side of it
             */
            Vector vab = bb.pos1() - aa.pos1();
            
#if GRID_HAS_PERIODIC
            if ( modulo )
                modulo->fold(vab);
#endif
            if ( dot(vab, aa.prevDiff()) >= 0 )
            {
                real ab2 = vab.normSqr();
                if ( below(ab2, ran) )
                    meca.addLongLink1(aa.vertex1(), bb.vertex1(), vab, ab2, ran, push);
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
    //std::clog << "   LL2 " << aa.obj_ << " " << bb.vertex2() << '\n';
    assert_true( aa.obj_->tag() == Fiber::TAG );
    assert_true( bb.obj_->tag() == Fiber::TAG );

    const real ran = aa.rad_ + bb.rad_;
    
    // get position of bb.vertex2() with respect to segment 'aa'
    real dis2 = INFINITY;
    FiberSegment seg = aa.segment();
    real abs = seg.projectPoint0(bb.pos2(), dis2);
    
    if ((0 <= abs) & (abs <= aa.len()))
    {
        /*
         bb.vertex2() projects inside segment 'aa'
         */
        if ( below(dis2, ran) )
            meca.addSideSlidingLink(seg, abs, bb.vertex2(), ran, push);
    }
    else if ( abs < 0 )
    {
        /*
         Check the projection to the segment located before 'aa',
         and interact if 'bb.vertex1()' falls on the right side of it
         */
        Vector vab = bb.pos2() - aa.pos1();
        
#if GRID_HAS_PERIODIC
        if ( modulo )
            modulo->fold(vab);
#endif
        if ( aa.isFirst() )
        {
            assert_true(bb.isLast());
            real ab2 = vab.normSqr();

            if ( below(ab2, ran)  && dot(vab, bb.diff()) <= 0 )
                meca.addLongLink1(aa.vertex1(), bb.vertex2(), vab, ab2, ran, push);
        }
        else
        {
            if ( dot(vab, aa.prevDiff()) >= 0 )
            {
                real ab2 = vab.normSqr();
                if ( below(ab2, ran) )
                    meca.addLongLink1(aa.vertex1(), bb.vertex2(), vab, ab2, ran, push);
            }
        }
    }
    else if (( &bb < &aa ) & aa.isLast() )
    {
        /*
         Check the projection of aa.vertex2(),
         on the segment represented by 'bb'
         */
        assert_true(abs > aa.len());
        assert_true(bb.isLast());
        
        Vector vab = bb.pos2() - aa.pos2();
#if GRID_HAS_PERIODIC
        if ( modulo )
            modulo->fold(vab);
#endif
        real ab2 = vab.normSqr();
        
        if ( below(ab2, ran)  &&  dot(vab, bb.diff()) <= 0 )
            meca.addLongLink1(aa.vertex2(), bb.vertex2(), vab, ab2, ran, push);
    }
}


/**
 This is used to check two FiberSegment, that each represent a segment of a Fiber.
 The segments are tested for intersection in 3D.
 */
void LocusGrid::checkLL(BigLocus const& aa, BigLocus const& bb) const
{
    assert_true( aa.obj_->tag() == Fiber::TAG );
    assert_true( bb.obj_->tag() == Fiber::TAG );

#if ( DIM >= 3 )
    
    const real ran = aa.rad_ + bb.rad_;
    
    /* in 3D, check the shortest distance between two segments, and if close
     enough, use the result to build an interaction */
    FiberSegment as = aa.segment();
    FiberSegment bs = bb.segment();
    real a, b;
    /* We do not need to calculate `a` and `b` if the distance 'dis'
     is greater than 'ran' since nothing will be done in that case... */
    real dis2 = as.shortestDistanceSqr(bs, a, b);
    
    if ( below(dis2, ran) & as.within(a) & bs.within(b) )
        meca.addSideSlidingLink(as, a, Interpolation(bs, b), ran, push);
    
#endif

    //std::clog << "   LL " << aa.obj_ << " " << bb.obj_ << '\n';
    checkLL1(aa, bb);
    
    if ( aa.isLast() )
        checkLL2(bb, aa);
    
    checkLL1(bb, aa);
    
    if ( bb.isLast() )
        checkLL2(aa, bb);
}

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
    return b.obj_->buddy() != a.obj_->buddy();
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
 objects that are too far appart to interact, based on BigVector::near()
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
 objects that are too far appart to interact, based on BigVector::near()
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
        const BigVector pos = ii->pos_;
        
        for ( BigLocus const* jj = list2.begin(); jj < mid2; ++jj )
            if ( pos.near(jj->pos_) && not_adjacentLL(*ii, *jj)  )
                checkLL(*ii, *jj);

        for ( BigPoint const* jj = mid2; jj < list2.end(); ++jj )
            if ( pos.near(jj->pos_) && not_adjacentPL(*jj, *ii) )
                checkPL(*jj, *ii);
    }

    for ( BigPoint const* ii = mid1; ii < list1.end(); ++ii )
    {
        const BigVector pos = ii->pos_;

        for ( BigLocus const* jj = list2.begin(); jj < mid2; ++jj )
            if ( pos.near(jj->pos_) && not_adjacentPL(*ii, *jj) )
                checkPL(*ii, *jj);

        for ( BigPoint const* jj = mid2; jj < list2.end(); ++jj )
            if ( pos.near(jj->pos_) && not_adjacentPP(*ii, *jj) )
                checkPP(*ii, *jj);
    }
}

#if ( DIM == 3 ) && USE_SIMD

/**
 Evaluate 4 pos.near(jj->pos_) using SIMD instructions
 @return a 4-bit integer where each bit represents the result of one test
 */
inline int four_near(vec4f const& xyzr, BigLocus const* src)
{
    vec4f tt = sub4f(xyzr, loadu4f(src[0].pos_.data()));
    vec4f yy = sub4f(xyzr, loadu4f(src[1].pos_.data()));
    vec4f uu = sub4f(xyzr, loadu4f(src[2].pos_.data()));
    vec4f rr = sub4f(xyzr, loadu4f(src[3].pos_.data()));
    // transpose 4x4 data matrix:
    vec4f xx = unpacklo4f(tt, yy);
    tt = unpackhi4f(tt, yy);
    vec4f zz = unpacklo4f(uu, rr);
    rr = unpackhi4f(uu, rr);
    yy = movelh4f(xx, zz);
    xx = movehl4f(zz, xx);
    uu = mul4f(yy, yy);
    zz = movelh4f(tt, rr);  // yy, rr
    rr = movehl4f(rr, tt);
    tt = mul4f(zz, zz);
    // calculate test:
    uu = add4f(mul4f(xx, xx), uu); // x*x + y*y
    tt = sub4f(mul4f(rr, rr), tt); // r*r - z*z
    return lower_mask4f(uu, tt);  // x*x + y*y < r*r - z*z
}

typedef unsigned long BitField;

/**
Evaluate `cnt` pos.near(jj->pos_) using SIMD instructions
@return a bitfield representing the result of all tests
*/
BitField near_bits(vec4f const& xyzr, BigLocus const* start, int cnt)
{
    BitField res = 0;
    unsigned shift = 0;
    BigLocus const* ptr = start;
    BigLocus const* end;
#if 0
    end = start + ( cnt & ~15UL );
    /*
     Unrolling can help since all four_near() are independent,
     and the calculations can be executed out-of-order efficiently,
     but this is effective only for list size > 32
     */
    while ( ptr < end )
    {
        BitField t = four_near(xyzr, ptr);
        BitField u = four_near(xyzr, ptr+4) << 4;
        BitField v = four_near(xyzr, ptr+8) << 8;
        BitField w = four_near(xyzr, ptr+12) << 12;
        res |= (( t | u )|( v | w )) << shift;
        shift += 16;
        ptr += 16;
    }
#endif
    end = start + ( cnt & ~7UL );
    while ( ptr < end )
    {
        BitField t = four_near(xyzr, ptr);
        BitField u = four_near(xyzr, ptr+4) << 4;
        res |= ( t | u ) << shift;
        shift += 8;
        ptr += 8;
    }
    
    end = start + cnt;
    while ( ptr < end )
    {
        BitField t = four_near(xyzr, ptr);
        unsigned i = end - ptr;
        assert_true( i < 8 );
        // a mask to clear the bits past the end:
        BitField k = ~( ~0LU << i );
        res |= ( t & k ) << shift;
        shift += 4;
        ptr += 4;
    }
    
    //printf(" near_bits %2lu : %lu\n", cnt, res);
    return res;
}


/**
 Set bitL corresponding to BigLocus in first part of list,
 and bitP corresponding to BigPoints in second part of list
*/
void near_bits(BitField& bitL, BitField& bitP, vec4f const& xyzr, BigLocusList const& list, int start)
{
    int cnt = std::min((int)list.size()-start, 64);
    BitField bits = near_bits(xyzr, list.begin()+start, cnt);
    int nloc = list.num_locus();
    int shift = nloc - std::min(start, nloc);
    if ( shift < 64 )
    {
        BitField mask = ~0UL << shift;
        bitP = bits & mask;
        bitL = bits & ~mask;
    } else {
        bitL = bits;
        bitP = 0;
    }
}


/**
 This will consider once all pairs of objects from the given lists,
 assuming that the list are different and no object is repeated.

 Compared to `setSterics0()`, this performs additional tests to exclude
 objects that are too far appart to interact, based on BigVector::near()
 
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
        size_t l1 = list1.num_locus(), p1 = list1.num_points();
        size_t l2 = list2.num_locus(), p2 = list2.num_points();
        printf(" stericsU: %2lu+%2lu  :  %2lu+%2lu\n", l1, p1, l2, p2);
    }
#endif
    BigLocus const* mid1 = list1.middle();
    constexpr BitField mask(~1UL);
    BitField bitP, bitL;
    int b, u;
    
    for ( BigLocus const* ii = list1.begin(); ii < mid1; ++ii )
    {
        BigVector pos = ii->pos_;
        vec4f xyzr { pos.XX, pos.YY, pos.ZZ, -pos.RR };
        /* In most situations, the list size would be < 64 and one round would
         be sufficient, but in all generality we must handle larger list size */
        for ( int offset = 0; offset < (int)list2.size(); offset += 64 )
        {
            near_bits(bitL, bitP, xyzr, list2, offset);
            //printf(" L%lX:%lX\n", bitP, bitL);
            BigLocus const* jj = list2.begin();
#if 0
            // verify all tests:
            for ( b = 0; b < std::min(64, (int)list2.size()-offset); ++b )
            {
                BigVector vec = (jj+b+offset)->pos_;
                bool n = pos.near(vec);
                bool p = (bitL+bitP) & ( 1UL << b );
                float x = square(pos.XX-vec.XX) + square(pos.YY-vec.YY);
                float y = square(pos.ZZ-vec.ZZ) - square(pos.RR+vec.RR);
                // SIMD and scalar results can differ, near the edges (x+y ~ 0) :
                if ( n != p ) printf("!near %i+%i: %u%u %f\n", offset, b, n, p, x+y);
            }
#endif
            while ( bitL )
            {
                //printf(" LL%lX\n", bits);
                b = __builtin_ctzl(bitL);
                u = b + offset;
                if ( not_adjacentLL(*ii, *(jj+u)) )
                    checkLL(*ii, *(jj+u));
                bitL &= mask << b;
            }
            while ( bitP )
            {
                //printf(" PL%lX\n", bits);
                b = __builtin_ctzl(bitP);
                u = b + offset;
                if ( not_adjacentPL(*(jj+u), *ii) )
                    checkPL(*(jj+u), *ii);
                bitP &= mask << b;
            }
        }
    }
    
    for ( BigPoint const* ii = mid1; ii < list1.end(); ++ii )
    {
        BigVector pos = ii->pos_;
        vec4f xyzr { pos.XX, pos.YY, pos.ZZ, -pos.RR };
        /* In most situations, the list size would be < 64 and one round would
         be sufficient, but in all generality we must handle larger list size */
        for ( int offset = 0; offset < (int)list2.size(); offset += 64 )
        {
            near_bits(bitL, bitP, xyzr, list2, offset);
            //printf(" P%lX:%lX\n", bitP, bitL);
            BigLocus const* jj = list2.begin();
            while ( bitL )
            {
                //printf(" LP%lX\n", bitL);
                b = __builtin_ctzl(bitL);
                u = b + offset;
                if ( not_adjacentPL(*ii, *(jj+u)) )
                    checkPL(*ii, *(jj+u));
                bitL &= mask << b;
            }
            while ( bitP )
            {
                //printf(" PP%lX\n", bitP);
                b = __builtin_ctzl(bitP);
                u = b + offset;
                if ( not_adjacentPP(*ii, *(jj+u)) )
                    checkPP(*ii, *(jj+u));
                bitP &= mask << b;
            }
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
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int const* region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
        
        BigLocusList& base = cell_list(inx);
        setSterics0(base);
        
        for ( int reg = 1; reg < nr; ++reg )
            setSterics0(base, cell_list(inx+region[reg]));
    }
}


/** This calls setStericsT() */
void LocusGrid::setStericsT() const
{
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int const* region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
        
        BigLocusList& base = cell_list(inx);
        setStericsT(base);
        
        for ( int reg = 1; reg < nr; ++reg )
        {
#if ( DIM == 3 ) && USE_SIMD
            BigLocusList& side = cell_list(inx+region[reg]);
            if ( base.size() < side.size() )
                setStericsX(base, side);
            else
                setStericsX(side, base);
#else
            setStericsT(base, cell_list(inx+region[reg]));
#endif
        }
    }
}


void LocusGrid::setSterics() const
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
void LocusGrid::setSterics0(size_t pan) const
{
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int const* region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
       
        BigLocusList& base = cell_list(inx, pan);
        setSterics0(base);
        
        for ( int reg = 1; reg < nr; ++reg )
            setSterics0(base, cell_list(inx+region[reg], pan));
    }
}


void LocusGrid::setStericsT(size_t pan) const
{
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
         int const* region;
         int nr = pGrid.getRegion(region, inx);
         assert_true(region[0] == 0);
        
         BigLocusList& base = cell_list(inx, pan);
         setStericsT(base);
         
         for ( int reg = 1; reg < nr; ++reg )
             setStericsT(base, cell_list(inx+region[reg], pan));
    }
}


/**
 Check interactions between the FatPoints contained in Panes `pan` and `bim`,
 where ( pan1 != pan2 )
 */
void LocusGrid::setSterics0(size_t pan, size_t bim) const
{
    assert_true(pan != bim);
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int const* region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
        
        BigLocusList& base1 = cell_list(inx, pan);
        BigLocusList& base2 = cell_list(inx, bim);

        setSterics0(base1, base2);

        for ( int reg = 1; reg < nr; ++reg )
        {
            setSterics0(base1, cell_list(inx+region[reg], bim));
            setSterics0(base2, cell_list(inx+region[reg], pan));
        }
    }
}


void LocusGrid::setStericsT(size_t pan, size_t bim) const
{
    assert_true(pan != bim);
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int const* region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
        
        BigLocusList& base1 = cell_list(inx, pan);
        BigLocusList& base2 = cell_list(inx, bim);

        setStericsT(base1, base2);

        for ( int reg = 1; reg < nr; ++reg )
        {
            setStericsT(base1, cell_list(inx+region[reg], bim));
            setStericsT(base2, cell_list(inx+region[reg], pan));
        }
    }
}


void LocusGrid::setSterics(size_t pan) const
{
    if ( pGrid.isPeriodic() )
        setSterics0(pan);
    else
        setStericsT(pan);
}


void LocusGrid::setSterics(size_t pan1, size_t pan2) const
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

