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

extern Modulo const* modulo;

//------------------------------------------------------------------------------

LocusGrid::LocusGrid()
{
}


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
    pGrid.printSummary(Cytosim::log, "LocusGrid");
}


void LocusGrid::mark() const
{
    for ( size_t i = 0; i < pGrid.nbCells(); ++i )
        pGrid[i].mark();
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


/**
 This is used to check two spherical objects:
 Solid/Bead/Sphere or the terminal vertex (the tips) of a Fiber
 
 The force is applied if the objects are closer to the maximum
 of their specified range + radius.
 */
void LocusGrid::checkPP(Meca& meca, real stiff,
                        BigPoint const& aa, BigPoint const& bb)
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
    
    if ( ab2 < ran*ran )
        meca.addLongLink1(aa.vertex1(), bb.vertex1(), vab, ab2, ran, stiff);
}


/**
 This is used to check a segment of a fiber against a spherical object:
 Solid/Bead/Sphere or Fiber-tip.
 
 The force is applied if the objects are closer than the maximum of the two range + radius,
 and if the center of the sphere projects inside the segment.
 */
void LocusGrid::checkPL(Meca& meca, real stiff,
                        BigPoint const& aa, BigLocus const& bb)
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
            if ( ab2 < ran*ran )
                meca.addSideSlidingLink(bb.segment(), abs, aa.vertex1(), ran, stiff);
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
                if ( ab2 < ran*ran )
                    meca.addLongLink1(aa.vertex1(), bb.vertex2(), vab, ab2, ran, stiff);
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
        if ( ab2 < ran*ran && ( bb.isFirst() || dot(vab, bb.prevDiff()) <= 0 ))
            meca.addLongLink1(aa.vertex1(), bb.vertex1(), vab, ab2, ran, stiff);
    }
}


/**
 This is used to check a segment of a fiber against another segment of fiber,
 not including the terminal vertex of fibers.
 
 The interaction is applied only if the vertex projects 'inside' the segment.
 */
void LocusGrid::checkLL1(Meca& meca, real stiff,
                         BigLocus const& aa, BigLocus const& bb)
{
    //std::clog << "   LL1 " << aa.obj_ << " " << bb.vertex1() << '\n';
    assert_true( aa.obj_->tag() == Fiber::TAG );
    assert_true( bb.obj_->tag() == Fiber::TAG );

    const real ran = aa.rad_ + bb.rad_;
    
    // get position of bb.vertex1() with respect to segment 'aa'
    real dis2 = INFINITY;
    FiberSegment seg = aa.segment();
    real abs = seg.projectPoint0(bb.pos1(), dis2);
    
    if ((0 <= abs) & (abs <= aa.len()) & (dis2 < ran*ran))
    {
        /*
         bb.vertex1() projects inside segment 'aa'
         */
        meca.addSideSlidingLink(seg, abs, bb.vertex1(), ran, stiff);
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
                if ( ab2 < ran*ran  &&  dot(vab, bb.diff()) >= 0 )
                    meca.addLongLink1(aa.vertex1(), bb.vertex1(), vab, ab2, ran, stiff);
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
                if ( ab2 < ran*ran )
                    meca.addLongLink1(aa.vertex1(), bb.vertex1(), vab, ab2, ran, stiff);
            }
        }
    }
}


/**
 This is used to check a segment of a fiber against the terminal vertex of a fiber

 The interaction is applied only if the vertex projects 'inside' the segment.
 */
void LocusGrid::checkLL2(Meca& meca, real stiff,
                         BigLocus const& aa, BigLocus const& bb)
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
        if ( dis2 < ran*ran )
            meca.addSideSlidingLink(seg, abs, bb.vertex2(), ran, stiff);
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

            if ( ab2 < ran*ran  && dot(vab, bb.diff()) <= 0 )
                meca.addLongLink1(aa.vertex1(), bb.vertex2(), vab, ab2, ran, stiff);
        }
        else
        {
            if ( dot(vab, aa.prevDiff()) >= 0 )
            {
                real ab2 = vab.normSqr();
                if ( ab2 < ran*ran )
                    meca.addLongLink1(aa.vertex1(), bb.vertex2(), vab, ab2, ran, stiff);
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
        
        if ( ab2 < ran*ran  &&  dot(vab, bb.diff()) <= 0 )
            meca.addLongLink1(aa.vertex2(), bb.vertex2(), vab, ab2, ran, stiff);
    }
}


/**
 This is used to check two FiberSegment, that each represent a segment of a Fiber.
 The segments are tested for intersection in 3D.
 */
void LocusGrid::checkLL(Meca& meca, real stiff,
                        BigLocus const& aa, BigLocus const& bb)
{
    assert_true( aa.obj_->tag() == Fiber::TAG );
    assert_true( bb.obj_->tag() == Fiber::TAG );

#if ( DIM >= 3 )
    
    const real ran = aa.rad_ + bb.rad_;
    
    /* in 3D, check the shortest distance between two segments, and if close
     enough, use the result to build an interaction */
    real a, b, d;
    FiberSegment as = aa.segment();
    FiberSegment bs = bb.segment();
    
    if ( ! as.belowDistance(bs, ran, a, b, d) )
        return;
    
    if ( as.within(a) & bs.within(b) )
        meca.addSideSlidingLink(as, a, Interpolation(bs, b), ran, stiff);
    
#endif

    //std::clog << "   LL " << aa.obj_ << " " << bb.obj_ << '\n';
    checkLL1(meca, stiff, aa, bb);
    
    if ( aa.isLast() )
        checkLL2(meca, stiff, bb, aa);
    
    checkLL1(meca, stiff, bb, aa);
    
    if ( bb.isLast() )
        checkLL2(meca, stiff, aa, bb);
}

//------------------------------------------------------------------------------
#pragma mark - Selections of pairs excluded from Sterics

/*
 In general, these test will only exclude relatively rare pairs from interacting,
 and thus are less stringent than BigVector::near(): they should be tested after.
 */

/// excluding two spheres when they are from the same Solid
inline static bool not_adjacentPP(BigPoint const& a, BigPoint const& b)
{
    return a.obj_ != b.obj_;
}


/// excluding Fiber and Solid from the same Aster
inline static bool not_adjacentPL(BigPoint const& a, BigLocus const& b)
{
    //a->mec_->Buddy::print(std::clog);
    //b->obj_->Buddy::print(std::clog);
    return b.obj_->buddy() != a.obj_->buddy();
}


/// excluding segments that are adjacent on the same fiber, or protofilaments from Tubule
inline static bool not_adjacentLL(BigLocus const& a, BigLocus const& b)
{
#if FIBER_HAS_FAMILY
    return (( a.obj_->family_ != b.obj_->family_ )
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
void LocusGrid::setSterics0(Meca& meca, real stiff, BigLocusList const& list)
{
    BigLocus const* mid = list.middle();
    
    for ( BigLocus const* ii = list.begin(); ii < mid; ++ii )
    {
        for ( BigLocus const* jj = ii+1; jj < mid; ++jj )
            if ( not_adjacentLL(*ii, *jj) )
                checkLL(meca, stiff, *ii, *jj);
        
        for ( BigPoint const* kk = mid; kk < list.end(); ++kk )
            if ( not_adjacentPL(*kk, *ii) )
                checkPL(meca, stiff, *kk, *ii);
    }

    for ( BigPoint const* ii = mid; ii < list.end(); ++ii )
    {
        for ( BigPoint const* jj = ii+1; jj < list.end(); ++jj )
            if ( not_adjacentPP(*ii, *jj) )
                checkPP(meca, stiff, *ii, *jj);
    }
}


/**
 This will consider once all pairs of objects from the given lists,
 assuming that the list are different and no object is repeated
 */
void LocusGrid::setSterics0(Meca& meca, real stiff, BigLocusList const& list1,
                            BigLocusList const& list2)
{
    assert_true( &list1 != &list2 );
    BigLocus const* mid1 = list1.middle();
    BigLocus const* mid2 = list2.middle();

    for ( BigLocus const* ii = list1.begin(); ii < mid1; ++ii )
    {
        for ( BigLocus const* jj = list2.begin(); jj < mid2; ++jj )
            if ( not_adjacentLL(*ii, *jj)  )
                checkLL(meca, stiff, *ii, *jj);

        for ( BigPoint const* kk = mid2; kk < list2.end(); ++kk )
            if ( not_adjacentPL(*kk, *ii) )
                checkPL(meca, stiff, *kk, *ii);
    }

    for ( BigPoint const* ii = mid1; ii < list1.end(); ++ii )
    {
        for ( BigLocus const* jj = list2.begin(); jj < mid2; ++jj )
            if ( not_adjacentPL(*jj, *ii) )
                checkPL(meca, stiff, *jj, *ii);

        for ( BigPoint const* jj = mid2; jj < list2.end(); ++jj )
            if ( not_adjacentPP(*ii, *jj) )
                checkPP(meca, stiff, *ii, *jj);
    }
}


/**
 This will consider once all pairs of objects from the given lists.
 Compared to `setSterics0()`, this performs additional tests to exclude
 objects that are too far appart to interact, based on BigVector::near()
 */
void LocusGrid::setStericsT(Meca& meca, real stiff, BigLocusList const& list)
{
    BigLocus const* mid = list.middle();

    for ( BigLocus const* ii = list.begin(); ii < mid; ++ii )
    {
        const BigVector pos = ii->pos_;
        
        for ( BigLocus const* jj = ii+1; jj < mid; ++jj )
            if ( pos.near(jj->pos_) && not_adjacentLL(*ii, *jj) )
                checkLL(meca, stiff, *ii, *jj);
        
        for ( BigPoint const* jj = mid; jj < list.end(); ++jj )
            if ( pos.near(jj->pos_) && not_adjacentPL(*jj, *ii) )
                checkPL(meca, stiff, *jj, *ii);
    }

    for ( BigPoint const* ii = mid; ii < list.end(); ++ii )
    {
        const BigVector pos = ii->pos_;

        for ( BigPoint const* jj = ii+1; jj < list.end(); ++jj )
            if ( pos.near(jj->pos_) && not_adjacentPP(*ii, *jj) )
                checkPP(meca, stiff, *ii, *jj);
    }
}


/**
 This will consider once all pairs of objects from the given lists,
 assuming that the list are different and no object is repeated.

 Compared to `setSterics0()`, this performs additional tests to exclude
 objects that are too far appart to interact, based on BigVector::near()
*/
void LocusGrid::setStericsT(Meca& meca, real stiff, BigLocusList const& list1,
                            BigLocusList const& list2)
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
                checkLL(meca, stiff, *ii, *jj);

        for ( BigPoint const* jj = mid2; jj < list2.end(); ++jj )
            if ( pos.near(jj->pos_) && not_adjacentPL(*jj, *ii) )
                checkPL(meca, stiff, *jj, *ii);
    }

    for ( BigPoint const* ii = mid1; ii < list1.end(); ++ii )
    {
        const BigVector pos = ii->pos_;

        for ( BigLocus const* jj = list2.begin(); jj < mid2; ++jj )
            if ( pos.near(jj->pos_) && not_adjacentPL(*ii, *jj) )
                checkPL(meca, stiff, *ii, *jj);

        for ( BigPoint const* jj = mid2; jj < list2.end(); ++jj )
            if ( pos.near(jj->pos_) && not_adjacentPP(*ii, *jj) )
                checkPP(meca, stiff, *ii, *jj);
    }
}


#if ( DIM == 3 ) && defined(__SSE3__)

/**
 Evaluate 4 pos.near(jj->pos_) using SIMD instructions
 @return a 4-bit integer where each bit represents the result of one test
 */
inline int four_near(vec4f const& xyzr, BigLocus const* src)
{
    vec4f tt = sub4f(xyzr, loadu4f(src->pos_.data()));
    vec4f yy = sub4f(xyzr, loadu4f((src+1)->pos_.data()));
    vec4f uu = sub4f(xyzr, loadu4f((src+2)->pos_.data()));
    vec4f rr = sub4f(xyzr, loadu4f((src+3)->pos_.data()));
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
    return movemask4f(_mm_cmplt_ps(uu, tt));  // x*x + y*y < r*r - z*z
}


/**
 This will consider once all pairs of objects from the given lists,
 assuming that the list are different and no object is repeated.

 Compared to `setSterics0()`, this performs additional tests to exclude
 objects that are too far appart to interact, based on BigVector::near()
*/
void LocusGrid::setStericsU(Meca& meca, real stiff, BigLocusList const& list1,
                            BigLocusList const& list2)
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
    for ( BigLocus const* ii = list1.begin(); ii < mid1; ++ii )
    {
        BigVector pos = ii->pos_;
        vec4f xyzr { pos.XX, pos.YY, pos.ZZ, -pos.RR };
        int t;
        
        BigLocus const* jj = list2.begin();
        #pragma nounroll
        for ( ; jj < list2.pre_middle(); jj += 4 )
        {
            t = four_near(xyzr, jj);
            if ( !t ) continue;
#if 0
            bool n0 = pos.near(jj->pos_);
            bool n1 = pos.near((jj+1)->pos_);
            bool n2 = pos.near((jj+2)->pos_);
            bool n3 = pos.near((jj+3)->pos_);
            printf("%i:%u%u%u%u \n", t, n3, n2, n1, n0);
#endif
            if ( (t&1) && not_adjacentLL(*ii, *(jj  )) ) checkLL(meca, stiff, *ii, *(jj));
            if ( (t&2) && not_adjacentLL(*ii, *(jj+1)) ) checkLL(meca, stiff, *ii, *(jj+1));
            if ( (t&4) && not_adjacentLL(*ii, *(jj+2)) ) checkLL(meca, stiff, *ii, *(jj+2));
            if ( (t&8) && not_adjacentLL(*ii, *(jj+3)) ) checkLL(meca, stiff, *ii, *(jj+3));
        }
        if ( jj < list2.pre_end() )
        {
            // handle transition zone with BigLocus and BigPoint
            t = four_near(xyzr, jj);
            if ( t )
            {
                //printf(" XL%i");
                //int u = 0, k = 1;
                int u = __builtin_ctz(t), k = 1 << u;
                int m = std::min(4, int(list2.middle()-jj));
                int e = std::min(4, int(list2.end()-jj));
                for ( ; u < m; ++u, k<<=1 )
                    if ( (t&k) && not_adjacentLL(*ii, *(jj+u)) ) checkLL(meca, stiff, *ii, *(jj+u));
                for ( ; u < e; ++u, k<<=1 )
                    if ( (t&k) && not_adjacentPL(*(jj+u), *ii) ) checkPL(meca, stiff,*(jj+u), *ii);
            }
            jj += 4;
            assert_true( jj > list2.middle() );
            #pragma nounroll
            for ( ; jj < list2.pre_end(); jj += 4 )
            {
                t = four_near(xyzr, jj);
                if ( !t ) continue;
                //printf(" LP%i", t);
                if ( (t&1) && not_adjacentPL(*(jj  ), *ii) ) checkPL(meca, stiff, *(jj), *ii);
                if ( (t&2) && not_adjacentPL(*(jj+1), *ii) ) checkPL(meca, stiff, *(jj+1), *ii);
                if ( (t&4) && not_adjacentPL(*(jj+2), *ii) ) checkPL(meca, stiff, *(jj+2), *ii);
                if ( (t&8) && not_adjacentPL(*(jj+3), *ii) ) checkPL(meca, stiff, *(jj+3), *ii);
            }
        }
        else
        {
            // less than 3 BigLocus remaining
            assert_true( list2.middle() - jj < 4 );
            #pragma nounroll
            for ( ; jj < list2.middle(); ++jj )
                if ( pos.near(jj->pos_) && not_adjacentLL(*ii, *jj) ) checkLL(meca, stiff, *ii, *jj);
        }
        // less than 3 BigPoint remaining
        assert_true( list2.end() - jj < 4 );
        #pragma nounroll
        for ( ; jj < list2.end(); ++jj )
            if ( pos.near(jj->pos_) && not_adjacentPL(*jj, *ii) ) checkPL(meca, stiff, *jj, *ii);
    }
    
    for ( BigPoint const* ii = mid1; ii < list1.end(); ++ii )
    {
        BigVector pos = ii->pos_;
        vec4f xyzr { pos.XX, pos.YY, pos.ZZ, -pos.RR };
        int t;
        
        BigLocus const* jj = list2.begin();
        #pragma nounroll
        for ( ; jj < list2.pre_middle(); jj += 4 )
        {
            t = four_near(xyzr, jj);
            if ( !t ) continue;
            //printf(" PL%i", t);
            if ( (t&1) && not_adjacentPL(*ii, *(jj  )) ) checkPL(meca, stiff, *ii, *(jj));
            if ( (t&2) && not_adjacentPL(*ii, *(jj+1)) ) checkPL(meca, stiff, *ii, *(jj+1));
            if ( (t&4) && not_adjacentPL(*ii, *(jj+2)) ) checkPL(meca, stiff, *ii, *(jj+2));
            if ( (t&8) && not_adjacentPL(*ii, *(jj+3)) ) checkPL(meca, stiff, *ii, *(jj+3));
        }
        if ( jj < list2.pre_end() )
        {
            // handle transition zone with BigLocus and BigPoint
            t = four_near(xyzr, jj);
            if ( t )
            {
                //printf(" XL%i");
                //int u = 0, k = 1;
                int u = __builtin_ctz(t), k = 1 << u;
                int m = int(list2.middle()-jj);
                assert_true( m < 4 );
                int e = std::min(4, int(list2.end()-jj));
                for ( ; u < m; ++u, k<<=1 )
                    if ( (t&k) && not_adjacentPL(*ii, *(jj+u)) ) checkPL(meca, stiff, *ii, *(jj+u));
                for ( ; u < e; ++u, k<<=1 )
                    if ( (t&k) && not_adjacentPP(*ii, *(jj+u)) ) checkPP(meca, stiff, *ii, *(jj+u));
            }
            jj += 4;
            assert_true( jj > list2.middle() );
            #pragma nounroll
            for ( ; jj < list2.pre_end(); jj += 4 )
            {
                t = four_near(xyzr, jj);
                if ( !t ) continue;
                //printf(" PP%i", t);
                if ( (t&1) && not_adjacentPP(*ii, *(jj  )) ) checkPP(meca, stiff, *ii, *(jj));
                if ( (t&2) && not_adjacentPP(*ii, *(jj+1)) ) checkPP(meca, stiff, *ii, *(jj+1));
                if ( (t&4) && not_adjacentPP(*ii, *(jj+2)) ) checkPP(meca, stiff, *ii, *(jj+2));
                if ( (t&8) && not_adjacentPP(*ii, *(jj+3)) ) checkPP(meca, stiff, *ii, *(jj+3));
            }
        }
        else
        {
            // less than 3 BigLocus remaining
            assert_true( list2.middle() - jj < 4 );
            #pragma nounroll
            for ( ; jj < list2.middle(); ++jj )
                if ( pos.near(jj->pos_) && not_adjacentPL(*ii, *jj) ) checkPL(meca, stiff, *ii, *jj);
        }
        // less than 3 BigPoint remaining
        assert_true( list2.end() - jj < 4 );
        #pragma nounroll
        for ( ; jj < list2.end(); ++jj )
            if ( pos.near(jj->pos_) && not_adjacentPP(*jj, *ii) ) checkPP(meca, stiff, *jj, *ii);
    }
}

typedef unsigned long bitfield;

/**
Evaluate `cnt` pos.near(jj->pos_) using SIMD instructions
@return a bitfield representint the result of all tests
*/
bitfield near_bits(vec4f const& xyzr, BigLocus const* start, int cnt)
{
    bitfield res = 0;
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
        bitfield t = four_near(xyzr, ptr);
        bitfield u = four_near(xyzr, ptr+4) << 4;
        bitfield v = four_near(xyzr, ptr+8) << 8;
        bitfield w = four_near(xyzr, ptr+12) << 12;
        res |= (( t | u )|( v | w )) << shift;
        shift += 16;
        ptr += 16;
    }
#endif
    end = start + ( cnt & ~7UL );
    while ( ptr < end )
    {
        bitfield t = four_near(xyzr, ptr);
        bitfield u = four_near(xyzr, ptr+4) << 4;
        res |= ( t | u ) << shift;
        shift += 8;
        ptr += 8;
    }
    
    end = start + cnt;
    while ( ptr < end )
    {
        bitfield t = four_near(xyzr, ptr);
        unsigned i = end - ptr;
        assert_true( i < 4 );
        // a mask to clear the bits past the end:
        bitfield k = ~( ~0LU << i );
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
void near_bits(bitfield& bitL, bitfield& bitP, vec4f const& xyzr, BigLocusList const& list, int start)
{
    int cnt = std::min((int)list.size()-start, 64);
    bitfield bits = near_bits(xyzr, list.begin()+start, cnt);
    int nloc = list.num_locus();
    int shift = nloc - std::min(start, nloc);
    if ( shift < 64 )
    {
        bitfield mask = ~0UL << shift;
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
*/
void LocusGrid::setStericsX(Meca& meca, real stiff, BigLocusList const& list1,
                            BigLocusList const& list2)
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
    constexpr bitfield mask(~1UL);
    bitfield bitP, bitL;
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
                if ( n != p ) printf("near %i+%i: %u%u %f\n", offset, b, n, p, x+y);
            }
#endif
            while ( bitL )
            {
                //printf(" LL%lX\n", bits);
                b = __builtin_ctzl(bitL);
                u = b + offset;
                if ( not_adjacentLL(*ii, *(jj+u)) )
                    checkLL(meca, stiff, *ii, *(jj+u));
                bitL &= mask << b;
            }
            while ( bitP )
            {
                //printf(" PL%lX\n", bits);
                b = __builtin_ctzl(bitP);
                u = b + offset;
                if ( not_adjacentPL(*(jj+u), *ii) )
                    checkPL(meca, stiff, *(jj+u), *ii);
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
                    checkPL(meca, stiff, *ii, *(jj+u));
                bitL &= mask << b;
            }
            while ( bitP )
            {
                //printf(" PP%lX\n", bitP);
                b = __builtin_ctzl(bitP);
                u = b + offset;
                if ( not_adjacentPP(*ii, *(jj+u)) )
                    checkPP(meca, stiff, *ii, *(jj+u));
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
 Scan all cells to examine each pair each pair of objects (ii, jj) only once.
 */
void LocusGrid::setSterics0(Meca& meca, real stiff) const
{
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int * region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
        
        BigLocusList& base = cell_list(inx);
        setSterics0(meca, stiff, base);
        
        for ( int reg = 1; reg < nr; ++reg )
            setSterics0(meca, stiff, base, cell_list(inx+region[reg]));
    }
}


/** This calls setStericsT() */
void LocusGrid::setStericsT(Meca& meca, real stiff) const
{
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int * region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
        
        BigLocusList& base = cell_list(inx);
        setStericsT(meca, stiff, base);
        
        for ( int reg = 1; reg < nr; ++reg )
        {
#if ( DIM == 3 ) && defined(__SSE3__)
            BigLocusList& side = cell_list(inx+region[reg]);
            if ( base.size() < side.size() )
                setStericsX(meca, stiff, base, side);
            else
                setStericsX(meca, stiff, side, base);
#else
            setStericsT(meca, stiff, base, cell_list(inx+region[reg]));
#endif
        }
    }
}


void LocusGrid::setInteractions(Meca& meca, real stiff) const
{
    //std::clog << "----" << '\n';
    if ( pGrid.isPeriodic() )
        setSterics0(meca, stiff);
    else
        setStericsT(meca, stiff);
}

#else

/**
 Check interactions between objects contained in the pane `pan`:
 Scan all cells to examine each pair each pair of objects (ii, jj) only once.
 */
void LocusGrid::setSterics0(Meca& meca, real stiff, size_t pan) const
{
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int * region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
       
        BigLocusList& base = cell_list(inx, pan);
        setSterics0(meca, stiff, base);
        
        for ( int reg = 1; reg < nr; ++reg )
            setSterics0(meca, stiff, base, cell_list(inx+region[reg], pan));
    }
}


void LocusGrid::setStericsT(Meca& meca, real stiff, size_t pan) const
{
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
         int * region;
         int nr = pGrid.getRegion(region, inx);
         assert_true(region[0] == 0);
        
         BigLocusList& base = cell_list(inx, pan);
         setStericsT(meca, stiff, base);
         
         for ( int reg = 1; reg < nr; ++reg )
             setStericsT(meca, stiff, base, cell_list(inx+region[reg], pan));
    }
}


/**
 Check interactions between the FatPoints contained in Panes `pan` and `bim`,
 where ( pan1 != pan2 )
 */
void LocusGrid::setSterics0(Meca& meca, real stiff, size_t pan, size_t bim) const
{
    assert_true(pan != bim);
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int * region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
        
        BigLocusList& base1 = cell_list(inx, pan);
        BigLocusList& base2 = cell_list(inx, bim);

        setSterics0(meca, stiff, base1, base2);

        for ( int reg = 1; reg < nr; ++reg )
        {
            setSterics0(meca, stiff, base1, cell_list(inx+region[reg], bim));
            setSterics0(meca, stiff, base2, cell_list(inx+region[reg], pan));
        }
    }
}


void LocusGrid::setStericsT(Meca& meca, real stiff, size_t pan, size_t bim) const
{
    assert_true(pan != bim);
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int * region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
        
        BigLocusList& base1 = cell_list(inx, pan);
        BigLocusList& base2 = cell_list(inx, bim);

        setStericsT(meca, stiff, base1, base2);

        for ( int reg = 1; reg < nr; ++reg )
        {
            setStericsT(meca, stiff, base1, cell_list(inx+region[reg], bim));
            setStericsT(meca, stiff, base2, cell_list(inx+region[reg], pan));
        }
    }
}


void LocusGrid::setInteractions(Meca& meca, real stiff, size_t pan) const
{
    if ( pGrid.isPeriodic() )
        setSterics0(meca, stiff, pan);
    else
        setStericsT(meca, stiff, pan);
}


void LocusGrid::setInteractions(Meca& meca, real stiff, size_t pan1, size_t pan2) const
{
    if ( pGrid.isPeriodic() )
        setSterics0(meca, stiff, pan1, pan2);
    else
        setStericsT(meca, stiff, pan1, pan2);
}

#endif

//------------------------------------------------------------------------------
#pragma mark - Display

#ifdef DISPLAY

#include "opengl.h"
void drawBoundaries(Map<DIM> const&);

void LocusGrid::drawGrid() const
{
#if ( DIM <= 3 )
    glDisable(GL_LIGHTING);
    glColor3f(1, 0, 0);
    glLineWidth(0.25);
    drawBoundaries(pGrid);
#endif
}
#endif

