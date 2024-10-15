// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "assert_macro.h"
#include "point_grid.h"
#include "exceptions.h"
#include "messages.h"
#include "modulo.h"
#include "space.h"
#include "solid.h"
#include "meca.h"

//------------------------------------------------------------------------------

index_t PointGrid::setGrid(Vector inf, Vector sup, real min_width)
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


void PointGrid::createCells()
{
    if ( pGrid.createCells() )
        pGrid.printSummary(Cytosim::log, "   PointGrid");

    // create side regions suitable for pairwise interactions:
    if ( !pGrid.hasRegions() )
        pGrid.createSideRegions(1);
}


size_t PointGrid::capacity() const
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
 
 The force is applied if the objects are closer than the
 sum of their radiuses.
 */
void PointGrid::checkPP(FatPoint const& aa, FatPoint const& bb) const
{
    Vector vab = bb.cen() - aa.cen();
    const real len = aa.rad_ + bb.rad_;
    const real ran = std::max(aa.rge_+bb.rad_, aa.rad_+bb.rge_);
    if ( modulo )
        modulo->fold(vab);
    real ab2 = vab.normSqr();
    if ( below(ab2, ran) )
    {
        real stiff = sign_select(ab2-len*len, push_, pull_);
        meca_.addLongLink(aa.pnt_, bb.pnt_, len, stiff);
        //meca_.addLongLink1(aa.pnt_, bb.pnt_, vab, ab2, len, push_);
    }
    //std::clog << "   PP- " << bb.pnt_ << " " << aa.pnt_ << "  " << std::sqrt(ab2) << '\n';
}


/**
 This is used to check a segment of a fiber against a spherical object:
 Solid/Bead/Sphere or Fiber-tip.
 
 The force is applied if the objects are closer than the sum of their radiuses.
 */
void PointGrid::checkPL(FatPoint const& aa, FatLocus const& bb) const
{
    //std::clog << "   PL- " << bb.seg_ << " " << aa.pnt_ << '\n';
    const real ran = aa.rad_ + bb.rad_;
    
    // get position of point with respect to segment:
    real ab2 = INFINITY;
    real abs = bb.seg_.projectPoint0(aa.cen(), ab2);
    
    if ( 0 <= abs )
    {
        if ( abs <= bb.len() )
        {
            // the point projects inside the segment
            if ( below(ab2, ran) )
                meca_.addSideSlidingLink(bb.seg_, abs, aa.pnt_, ran, push_);
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
                    meca_.addLongLink1(aa.pnt_, bb.vertex2(), vab, ab2, ran, push_);
            }
        }
    }
    else
    {
        
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
            meca_.addLongLink1(aa.pnt_, bb.vertex1(), vab, ab2, ran, push_);
    }
}


/**
 This is used to check a segment of a fiber against another segment of fiber,
 not including the terminal vertex of fibers.
 
 The interaction is applied only if the vertex projects 'inside' the segment.
 */
void PointGrid::checkLL1(FatLocus const& aa, FatLocus const& bb) const
{
    //std::clog << "   LL1 " << aa.seg_ << " " << bb.vertex1() << '\n';
    const real ran = std::max(aa.rge_+bb.rad_, aa.rad_+bb.rge_);
    
    // get position of bb.vertex1() with respect to segment 'aa'
    real dis2 = INFINITY;
    real abs = aa.seg_.projectPoint0(bb.pos1(), dis2);
    
    if ( below(dis2, ran) && ((0 <= abs) & (abs <= aa.len())))
    {
        /*
         bb.vertex1() projects inside segment 'aa'
         */
        const real len = aa.rad_ + bb.rad_;
        real stiff = sign_select(dis2-len*len, push_,pull_);
        meca_.addSideSlidingLink(aa.seg_, abs, bb.vertex1(), len, stiff);
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
                real len = aa.rad_ + bb.rad_;
                if ( below(ab2, len)  &&  dot(vab, bb.diff()) >= 0 )
                    meca_.addLongLink1(aa.vertex1(), bb.vertex1(), vab, ab2, len, push_);
            }
        }
        else
        {
            /*
             Check the projection to the segment located before 'aa',
             and interact if 'bb.vertex1()' falls on the right side of it
             */
            Vector vab = bb.pos1() - aa.pos1();
            if ( modulo )
                modulo->fold(vab);
            if ( dot(vab, aa.prevDiff()) >= 0 )
            {
                real ab2 = vab.normSqr();
                if ( below(ab2, ran) )
                {
                    real len = aa.rad_ + bb.rad_;
                    real stiff = sign_select(ab2-len*len, push_,pull_);
                    meca_.addLongLink2(aa.vertex1(), bb.vertex1(), vab, ab2, len, stiff);
                }
            }
        }
    }
}


/**
 This is used to check a segment of a fiber against the terminal vertex of a fiber
 
 The interaction is applied only if the vertex projects 'inside' the segment.
 */
void PointGrid::checkLL2(FatLocus const& aa, FatLocus const& bb) const
{
    //std::clog << "   LL2 " << aa.seg_ << " " << bb.vertex2() << '\n';
    const real ran = std::max(aa.rge_+bb.rad_, aa.rad_+bb.rge_);
    
    // get position of bb.vertex2() with respect to segment 'aa'
    real dis2 = INFINITY;
    real abs = aa.seg_.projectPoint0(bb.pos2(), dis2);
    
    if ((0 <= abs) & (abs <= aa.len()))
    {
        /*
         bb.vertex2() projects inside segment 'aa'
         */
        if ( below(dis2, ran) )
        {
            const real len = aa.rad_ + bb.rad_;
            real stiff = sign_select(dis2-len*len, push_,pull_);
            meca_.addSideSlidingLink(aa.seg_, abs, bb.vertex2(), len, stiff);
        }
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
            real len = aa.rad_ + bb.rad_;
            if ( below(ab2, len)  &&  dot(vab, bb.diff()) <= 0 )
                meca_.addLongLink1(aa.vertex1(), bb.vertex2(), vab, ab2, len, push_);
        }
        else
        {
            if ( dot(vab, aa.prevDiff()) >= 0 )
            {
                real ab2 = vab.normSqr();
                if ( below(ab2, ran) )
                {
                    real len = aa.rad_ + bb.rad_;
                    real stiff = sign_select(ab2-len*len, push_,pull_);
                    meca_.addLongLink2(aa.vertex1(), bb.vertex2(), vab, ab2, len, stiff);
                }
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
        if ( modulo )
            modulo->fold(vab);
        real ab2 = vab.normSqr();
        real len = aa.rad_ + bb.rad_;
        if ( below(ab2, len)  &&  dot(vab, bb.diff()) <= 0 )
            meca_.addLongLink1(aa.vertex2(), bb.vertex2(), vab, ab2, len, push_);
    }
}


/**
 This is used to check two FiberSegment, each representing the segment of a Fiber.
 in 3D, the segments are tested for getting within the requested distance.
 in 2D, only the vertices are checked.
 */
void PointGrid::checkLL(FatLocus const& aa, FatLocus const& bb) const
{
#if ( DIM >= 3 )
    
    const real ran = std::max(bb.rad_+aa.rge_, aa.rad_+bb.rge_);

    /* in 3D, check the shortest distance between two segments, and if close
     enough, use the result to build an interaction */
    real a, b;
    real dis2 = aa.seg_.shortestDistanceSqr(bb.seg_, a, b);
    
    if ( dis2 < ran * ran )
    {
        if ( aa.seg_.within(a) && bb.seg_.within(b) )
        {
            const real len = aa.rad_ + bb.rad_;
            real stiff = sign_select(dis2-len*len, push_,pull_);
            meca_.addSideSlidingLink(aa.seg_, a, Interpolation(bb.seg_, b), len, stiff);
        }
    }
    else
    {
        /* If the shortest distance between the lines is greater than 'ran', then
         the vertices associated with this segment will also be too far to interact */
        return;
    }
#endif
    
    //std::clog << "LL " << aa.seg_ << " " << bb.seg_ << '\n';
    checkLL1(aa, bb);
    
    if ( aa.isLast() )
        checkLL2(bb, aa);
    
    checkLL1(bb, aa);
    
    if ( bb.isLast() )
        checkLL2(aa, bb);
}


//------------------------------------------------------------------------------
#pragma mark - Functions to exclude certain pairs from Sterics

/*
 In general, these test will only exclude relatively rare pairs from interacting,
 and thus are less stringent than FatVector::near(): they should be tested after.
 */

/// excluding two spheres when they are from the same Solid
static inline bool not_adjacentPP(FatPoint const* a, FatPoint const* b)
{
    return a->pnt_.mecable() != b->pnt_.mecable();
}


/// excluding Fiber and Solid from the same Aster
static inline bool not_adjacentPL(FatPoint const* a, FatLocus const* b)
{
    //a->pnt_.mecable()->Buddy::print(std::clog);
    //b->seg_.fiber()->Buddy::print(std::clog);
    return ! b->seg_.fiber()->isBuddy(a->pnt_.mecable());
}


/// excluding segments that are adjacent on the same fiber, or protofilaments from Tubule
static inline bool not_adjacentLL(FatLocus const* a, FatLocus const* b)
{
#if FIBER_HAS_FAMILY
    return (( a->seg_.fiber()->family_ != b->seg_.fiber()->family_ )
            || (( a->seg_.point() > 1 + b->seg_.point() ) || ( b->seg_.point() > 1 + a->seg_.point() )));
#else
    return (( a->seg_.fiber() != b->seg_.fiber() )
            || (( a->seg_.point() > 1 + b->seg_.point() ) || ( b->seg_.point() > 1 + a->seg_.point() )));
#endif
    // we cannot use abs() above because `FatLocus::point()` is unsigned
}

//------------------------------------------------------------------------------
#pragma mark - Check all possible object pairs from two Cells

/**
 This will consider once all pairs of objects from the given lists
 */
void PointGrid::setSterics0(FatPointList & pots, FatLocusList & locs) const
{
    for ( FatPoint* ii = pots.begin(); ii < pots.end(); ++ii )
    {
        for ( FatPoint* jj = ii+1; jj < pots.end(); ++jj )
            if ( not_adjacentPP(ii, jj) )
                checkPP(*ii, *jj);
        
        for ( FatLocus* kk = locs.begin(); kk < locs.end(); ++kk )
            if ( not_adjacentPL(ii, kk) )
                checkPL(*ii, *kk);
    }

    for ( FatLocus* ii = locs.begin(); ii < locs.end(); ++ii )
    {
        for ( FatLocus* jj = ii+1; jj < locs.end(); ++jj )
            if ( not_adjacentLL(ii, jj) )
                checkLL(*ii, *jj);
    }
}


/**
 This will consider once all pairs of objects from the given lists,
 assuming that the list are different and no object is repeated
 */
void PointGrid::setSterics0(FatPointList & pots1, FatLocusList & locs1,
                            FatPointList & pots2, FatLocusList & locs2) const
{
    assert_true( &pots1 != &pots2 );
    assert_true( &locs1 != &locs2 );
    
    for ( FatPoint* ii = pots1.begin(); ii < pots1.end(); ++ii )
    {
        for ( FatPoint* jj = pots2.begin(); jj < pots2.end(); ++jj )
            if ( not_adjacentPP(ii, jj) )
                checkPP(*ii, *jj);
        
        for ( FatLocus* kk = locs2.begin(); kk < locs2.end(); ++kk )
            if ( not_adjacentPL(ii, kk) )
                checkPL(*ii, *kk);
    }
    
    for ( FatLocus* ii = locs1.begin(); ii < locs1.end(); ++ii )
    {
        for ( FatPoint* jj = pots2.begin(); jj < pots2.end(); ++jj )
            if ( not_adjacentPL(jj, ii) )
                checkPL(*jj, *ii);
        
        for ( FatLocus* kk = locs2.begin(); kk < locs2.end(); ++kk )
        {
            if ( not_adjacentLL(ii, kk)  )
                checkLL(*ii, *kk);
        }
    }
}


/**
 This will consider once all pairs of objects from the given lists.
 Compared to `setSterics0()`, this performs additional tests to exclude
 objects that are too far apart to interact, based on FatVector::near()
 */
void PointGrid::setStericsT(FatPointList & pots, FatLocusList & locs) const
{
    for ( FatPoint* ii = pots.begin(); ii < pots.end(); ++ii )
    {
        const FatVector pos = ii->pos_;
        
        for ( FatPoint* jj = ii+1; jj < pots.end(); ++jj )
            if ( pos.near(jj->pos_) && not_adjacentPP(ii, jj) )
                checkPP(*ii, *jj);
        
        for ( FatLocus* kk = locs.begin(); kk < locs.end(); ++kk )
            if ( pos.near(kk->pos_) && not_adjacentPL(ii, kk) )
                checkPL(*ii, *kk);
    }

    for ( FatLocus* ii = locs.begin(); ii < locs.end(); ++ii )
    {
        const FatVector pos = ii->pos_;
        
        for ( FatLocus* jj = ii+1; jj < locs.end(); ++jj )
            if ( pos.near(jj->pos_) && not_adjacentLL(ii, jj) )
                checkLL(*ii, *jj);
    }
}


/**
 This will consider once all pairs of objects from the given lists,
 assuming that the list are different and no object is repeated.

 Compared to `setSterics0()`, this performs additional tests to exclude
 objects that are too far apart to interact, based on FatVector::near()
 */
void PointGrid::setStericsT(FatPointList & pots1, FatLocusList & locs1,
                            FatPointList & pots2, FatLocusList & locs2) const
{
    assert_true( &pots1 != &pots2 );
    assert_true( &locs1 != &locs2 );

    for ( FatPoint* ii = pots1.begin(); ii < pots1.end(); ++ii )
    {
        const FatVector pos = ii->pos_;

        for ( FatPoint* jj = pots2.begin(); jj < pots2.end(); ++jj )
            if ( pos.near(jj->pos_) && not_adjacentPP(ii, jj) )
                checkPP(*ii, *jj);
        
        for ( FatLocus* kk = locs2.begin(); kk < locs2.end(); ++kk )
            if ( pos.near(kk->pos_) && not_adjacentPL(ii, kk) )
                checkPL(*ii, *kk);
    }
    
    for ( FatLocus* ii = locs1.begin(); ii < locs1.end(); ++ii )
    {
        const FatVector pos = ii->pos_;

        for ( FatPoint* jj = pots2.begin(); jj < pots2.end(); ++jj )
            if ( pos.near(jj->pos_) && not_adjacentPL(jj, ii) )
                checkPL(*jj, *ii);
        
        for ( FatLocus* kk = locs2.begin(); kk < locs2.end(); ++kk )
            if ( pos.near(kk->pos_) && not_adjacentLL(ii, kk) )
                checkLL(*ii, *kk);
    }
}

//------------------------------------------------------------------------------
#pragma mark - Check all pairs of Cells


#if ( NUM_STERIC_PANES == 1 )

/**
 Check interactions between objects contained in the grid.
 */
void PointGrid::setSterics0() const
{
    // scan all cells to examine each pair of particles:
    for ( index_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int const* region;
        int nR = pGrid.getRegion(region, inx);
        
        FatPointList & baseP = point_list(inx);
        FatLocusList & baseL = locus_list(inx);
        setSterics0(baseP, baseL);
        
        for ( int reg = 0; reg < nR; ++reg )
        {
            FatPointList & sideP = point_list(inx+region[reg]);
            FatLocusList & sideL = locus_list(inx+region[reg]);
            setSterics0(baseP, baseL, sideP, sideL);
        }
    }
}


void PointGrid::setStericsT() const
{
    assert_false( pGrid.isPeriodic() );
    
    // scan all cells to examine each pair of particles:
    for ( index_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int const* region;
        int nR = pGrid.getRegion(region, inx);
        
        FatPointList & baseP = point_list(inx);
        FatLocusList & baseL = locus_list(inx);
        
        setStericsT(baseP, baseL);
        
        for ( int reg = 0; reg < nR; ++reg )
        {
            FatPointList & sideP = point_list(inx+region[reg]);
            FatLocusList & sideL = locus_list(inx+region[reg]);
            setStericsT(baseP, baseL, sideP, sideL);
        }
    }
}


void PointGrid::setSterics() const
{
#if 0
    // printouts for debugging session, 15/04/2025
    pGrid.printRegions(std::clog, "PointGrid");
    for ( index_t i = 0; i < pGrid.nbCells(); ++i )
    {
        size_t nP = point_list(i).size();
        size_t nL = locus_list(i).size();
        std::clog << i << "  P " << nP << " L " << nL << "\n";
    }
#endif

    //std::clog << "----" << '\n';
    if ( pGrid.isPeriodic() )
        setSterics0();
    else
        setStericsT();
}

#else

/**
 Check interactions between the FatPoints contained in Pane `pan`.
 */
void PointGrid::setSterics0(const index_t pan) const
{
    // scan all cells to examine each pair of particles:
    for ( index_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int const* region;
        int nR = pGrid.getRegion(region, inx);
        
        FatPointList & baseP = point_list(inx, pan);
        FatLocusList & baseL = locus_list(inx, pan);
        
        setSterics0(baseP, baseL);
        
        for ( int reg = 0; reg < nR; ++reg )
        {
            FatPointList & sideP = point_list(inx+region[reg], pan);
            FatLocusList & sideL = locus_list(inx+region[reg], pan);
            setSterics0(baseP, baseL, sideP, sideL);
        }
    }
}


void PointGrid::setStericsT(const index_t pan) const
{
    assert_false( pGrid.isPeriodic() );
    // scan all cells to examine each pair of particles:
    for ( index_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int const* region;
        int nR = pGrid.getRegion(region, inx);
        
        FatPointList & baseP = point_list(inx, pan);
        FatLocusList & baseL = locus_list(inx, pan);
        
        setStericsT(baseP, baseL);
        
        for ( int reg = 0; reg < nR; ++reg )
        {
            FatPointList & sideP = point_list(inx+region[reg], pan);
            FatLocusList & sideL = locus_list(inx+region[reg], pan);
            setStericsT(baseP, baseL, sideP, sideL);
        }
    }
}


/**
 Check interactions between the FatPoints contained in Panes `pan1` and `pan2`,
 where ( pan1 != pan2 )
 */
void PointGrid::setSterics0(const index_t pan1, const index_t pan2) const
{
    assert_true(pan1 != pan2);
    
    // scan all cells to examine each pair of particles:
    for ( index_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int const* region;
        int nR = pGrid.getRegion(region, inx);
        
        FatPointList & baseP = point_list(inx, pan1);
        FatLocusList & baseL = locus_list(inx, pan1);
        
        FatPointList & baseP2 = point_list(inx, pan2);
        FatLocusList & baseL2 = locus_list(inx, pan2);

        setSterics0(baseP, baseL, baseP2, baseL2);

        for ( int reg = 0; reg < nR; ++reg )
        {
            FatPointList & sideP = point_list(inx+region[reg], pan2);
            FatLocusList & sideL = locus_list(inx+region[reg], pan2);
            setSterics0(baseP, baseL, sideP, sideL);
        }
        
        for ( int reg = 0; reg < nR; ++reg )
        {
            FatPointList & sideP = point_list(inx+region[reg], pan1);
            FatLocusList & sideL = locus_list(inx+region[reg], pan1);
            setSterics0(baseP2, baseL2, sideP, sideL);
        }
    }
}

void PointGrid::setStericsT(const index_t pan1, const index_t pan2) const
{
    assert_false( pGrid.isPeriodic() );
    assert_true(pan1 != pan2);
    
    // scan all cells to examine each pair of particles:
    for ( index_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int const* region;
        int nR = pGrid.getRegion(region, inx);
        
        FatPointList & baseP = point_list(inx, pan1);
        FatLocusList & baseL = locus_list(inx, pan1);
        
        FatPointList & baseP2 = point_list(inx, pan2);
        FatLocusList & baseL2 = locus_list(inx, pan2);

        setStericsT(baseP, baseL, baseP2, baseL2);

        for ( int reg = 0; reg < nR; ++reg )
        {
            FatPointList & sideP = point_list(inx+region[reg], pan2);
            FatLocusList & sideL = locus_list(inx+region[reg], pan2);
            setStericsT(baseP, baseL, sideP, sideL);
        }
        
        for ( int reg = 0; reg < nR; ++reg )
        {
            FatPointList & sideP = point_list(inx+region[reg], pan1);
            FatLocusList & sideL = locus_list(inx+region[reg], pan1);
            setStericsT(baseP2, baseL2, sideP, sideL);
        }
    }
}


void PointGrid::setSterics(index_t pan) const
{
    if ( pGrid.isPeriodic() )
        setSterics0(pan);
    else
        setStericsT(pan);
}


void PointGrid::setSterics(index_t pan1, index_t pan2) const
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

void PointGrid::drawGrid() const
{
#if ( DIM <= 3 )
    gym::ref_view();
    gym::disableLighting();
    gym::color(0,1,0);
    drawBoundaries(pGrid, 0.5f);
#endif
}
#endif

