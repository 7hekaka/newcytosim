// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "assert_macro.h"
#include "locus_grid.h"
#include "mecapoint.h"
#include "fiber_segment.h"
#include "exceptions.h"
#include "messages.h"
#include "modulo.h"
#include "space.h"
#include "meca.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------

inline Mecapoint BigPoint::point() const { return Mecapoint(mec_, pti_); }

inline FiberSegment BigLocus::segment() const { return FiberSegment(fib_, pti_); }

inline Mecapoint BigLocus::point1() const { return Mecapoint(fib_, pti_); }

inline Mecapoint BigLocus::point2() const { return Mecapoint(fib_, pti_+1); }

//------------------------------------------------------------------------------

LocusGrid::LocusGrid()
{
}


size_t LocusGrid::setGrid(Space const* spc, real min_step)
{
    assert_true( modulo == spc->getModulo() );
    assert_true( min_step > 0 );
    Vector inf, sup;
    spc->boundaries(inf, sup);
    
    size_t cnt[3] = { 1, 1, 1 };
    for ( size_t d = 0; d < DIM; ++d )
    {
        real n = ( sup[d] - inf[d] ) / min_step;
        
        if ( n < 0 )
            throw InvalidParameter("invalid space:boundaries");
        
#if GRID_HAS_PERIODIC
        if ( modulo  &&  modulo->isPeriodic(d) )
        {
            assert_small( modulo->period_[d] - sup[d] + inf[d] );
            //adjust the grid to match the edges
            cnt[d] = std::max((size_t)1, (size_t)std::floor(n));
            pGrid.setPeriodic(d, true);
        }
        else
#endif
        {
            //add a border in any dimension which is not periodic
            cnt[d] = (size_t)std::ceil(n) + 2;
            n = cnt[d] * 0.5 * min_step;
            real mid = inf[d] + sup[d];
            inf[d] = mid - n;
            sup[d] = mid + n;
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


size_t LocusGrid::capacity() const
{
    size_t res = 0;
    for ( size_t i = 0; i < pGrid.nbCells(); ++i )
        res += pGrid[i].capacity();
    return res;
}


//------------------------------------------------------------------------------
#pragma mark -

#if ( MAX_STERIC_PANES != 1 )

void LocusGrid::add(size_t pan, Mecable const* mec, size_t inx, real rad)
{
    if ( pan == 0 || pan > MAX_STERIC_PANES )
        throw InvalidParameter("point:steric is out-of-range");
    
    Vector w = mec->posPoint(inx);
    point_list(w, pan).emplace(mec, inx, rad, w);
}


void LocusGrid::add(size_t pan, Fiber const* fib, size_t inx, real rad, real rge)
{
    if ( pan == 0 || pan > MAX_STERIC_PANES )
        throw InvalidParameter("line:steric is out-of-range");
    
    // link in the cell containing the middle of the segment:
    Vector w = fib->posPoint(inx, 0.5);
    locus_list(w, pan).emplace(fib, inx, rad, rge, w);
}


#endif


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
    //std::clog << "   PP- " << bb.pnt << " " << aa.pnt << '\n';
    Vector vab = bb.cen() - aa.cen();
    const real ran = aa.rad() + bb.rad();

#if GRID_HAS_PERIODIC
    if ( modulo )
        modulo->fold(vab);
#endif
    real ab2 = vab.normSqr();
    
    if ( ab2 < ran*ran )
        meca.addLongLink1(aa.point(), bb.point(), vab, ab2, ran, stiff);
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
    //std::clog << "   PL- " << bb.seg << " " << aa.pnt << '\n';
    const real ran = aa.rad() + bb.rad_;
    
    // get position of point with respect to segment:
    real dis2 = INFINITY;
    real abs = bb.segment().projectPoint0(aa.cen(), dis2);
    
    if ( 0 <= abs )
    {
        if ( abs <= bb.len() )
        {
            if ( dis2 < ran*ran )
                meca.addSideSlidingLink(bb.segment(), abs, aa.point(), ran, stiff);
        }
        else
        {
            if ( bb.isLast() )
                checkPP(meca, stiff, aa, bb.bigPoint2());
        }
    }
    else
    {
        if ( bb.isFirst() )
            checkPP(meca, stiff, aa, bb.bigPoint1());
        else
        {
            /* we check the projection to the previous segment,
             and if it falls on the right of it, then we interact with the node */
            Vector vab = bb.pos1() - aa.cen();
            
#if GRID_HAS_PERIODIC
            if ( modulo )
                modulo->fold(vab);
#endif
            if ( dot(vab, bb.prevDiff()) <= 0 )
            {
                real ab2 = vab.normSqr();
                if ( ab2 < ran*ran )
                    meca.addLongLink1(aa.point(), bb.point1(), vab, ab2, ran, stiff);
            }
        }
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
    //std::clog << "   LL1 " << aa.seg << " " << bb.point1() << '\n';
    const real ran = aa.rad_ + bb.rad_;
    
    // get position of bb.point1() with respect to segment 'aa'
    real dis2 = INFINITY;
    FiberSegment seg = aa.segment();
    real abs = seg.projectPoint0(bb.pos1(), dis2);
    
    if ((0 <= abs) & (abs <= aa.len()) & (dis2 < ran*ran))
    {
        /*
         bb.point1() projects inside segment 'aa'
         */
        meca.addSideSlidingLink(seg, abs, bb.point1(), ran, stiff);
    }
    else if ( abs < 0 )
    {
        if ( aa.isFirst() )
        {
            /*
             Check the projection of aa.point1(),
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
                    meca.addLongLink1(aa.point1(), bb.point1(), vab, ab2, ran, stiff);
            }
        }
        else
        {
            /*
             Check the projection to the segment located before 'aa',
             and interact if 'bb.point1()' falls on the right side of it
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
                    meca.addLongLink1(aa.point1(), bb.point1(), vab, ab2, ran, stiff);
            }
        }
    }
}


/**
 This is used to check a segment of a fiber against the terminal vertex of fiber

 The interaction is applied only if the vertex projects 'inside' the segment.
 */
void LocusGrid::checkLL2(Meca& meca, real stiff,
                         BigLocus const& aa, BigLocus const& bb)
{
    //std::clog << "   LL2 " << aa.seg << " " << bb.point2() << '\n';
    const real ran = aa.rad_ + bb.rad_;
    
    // get position of bb.point2() with respect to segment 'aa'
    real dis2 = INFINITY;
    FiberSegment seg = aa.segment();
    real abs = seg.projectPoint0(bb.pos2(), dis2);
    
    if ((0 <= abs) & (dis2 < ran*ran) & (abs <= aa.len()))
    {
        /*
         bb.point2() projects inside segment 'aa'
         */
        meca.addSideSlidingLink(seg, abs, bb.point2(), ran, stiff);
    }
    else if ( abs < 0 )
    {
        /*
         Check the projection to the segment located before 'aa',
         and interact if 'bb.point1()' falls on the right side of it
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
                meca.addLongLink1(aa.point1(), bb.point2(), vab, ab2, ran, stiff);
        }
        else
        {
            if ( dot(vab, aa.prevDiff()) >= 0 )
            {
                real ab2 = vab.normSqr();
                if ( ab2 < ran*ran )
                    meca.addLongLink1(aa.point1(), bb.point2(), vab, ab2, ran, stiff);
            }
        }
    }
    else if ( &bb < &aa  &&  aa.isLast()  &&  abs > aa.len() )
    {
        /*
         Check the projection of aa.point2(),
         on the segment represented by 'bb'
         */
        assert_true(bb.isLast());
        
        Vector vab = bb.pos2() - aa.pos2();
#if GRID_HAS_PERIODIC
        if ( modulo )
            modulo->fold(vab);
#endif
        real ab2 = vab.normSqr();
        
        if ( ab2 < ran*ran  &&  dot(vab, bb.diff()) <= 0 )
            meca.addLongLink1(aa.point2(), bb.point2(), vab, ab2, ran, stiff);
    }
}


/**
 This is used to check two FiberSegment, that each represent a segment of a Fiber.
 The segments are tested for intersection in 3D.
 */
void LocusGrid::checkLL(Meca& meca, real stiff,
                        BigLocus const& aa, BigLocus const& bb)
{
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

    //std::clog << "LL " << aa.seg << " " << bb.seg << '\n';
    checkLL1(meca, stiff, aa, bb);
    
    if ( aa.isLast() )
        checkLL2(meca, stiff, bb, aa);
    
    checkLL1(meca, stiff, bb, aa);
    
    if ( bb.isLast() )
        checkLL2(meca, stiff, aa, bb);
}

//------------------------------------------------------------------------------
#pragma mark - Selections of pairs excluded from Sterics


/// excluding two spheres when they are from the same Solid
inline bool not_adjacent(BigPoint const* a, BigPoint const* b)
{
    return a->mec_ != b->mec_;
}


/// excluding Fiber and Solid from the same Aster
inline bool not_adjacent(BigPoint const* a, BigLocus const* b)
{
    //a->mec_->Buddy::print(std::clog);
    //b->fib_->Buddy::print(std::clog);
    return b->fib_->buddy() != a->mec_->buddy();
}


/// excluding segments that are adjacent on the same fiber, or protofilaments from Tubule
inline bool not_adjacent(BigLocus const* a, BigLocus const* b)
{
#if FIBER_HAS_FAMILY
    return (( a->fib_->family_ != b->fib_->family_ )
#else
    return (( a->fib_ != b->fib_ )
#endif
            || (( a->pti_ >= 2 + b->pti_ ) | ( b->pti_ >= 2 + a->pti_ )));
            // we cannot use abs() above because `pti_` is unsigned
}

//------------------------------------------------------------------------------
#pragma mark - Check all possible object pairs from two Cells

/**
 This will consider once all pairs of objects from the given lists
 */
void LocusGrid::setSterics0(Meca& meca, real stiff,
                            BigPointList & pots, BigLocusList & locs)
{
    for ( BigPoint* ii = pots.begin(); ii < pots.end(); ++ii )
    {
        for ( BigPoint* jj = ii+1; jj < pots.end(); ++jj )
            if ( not_adjacent(ii, jj) )
                checkPP(meca, stiff, *ii, *jj);
        
        for ( BigLocus* kk = locs.begin(); kk < locs.end(); ++kk )
            if ( not_adjacent(ii, kk) )
                checkPL(meca, stiff, *ii, *kk);
    }

    for ( BigLocus* ii = locs.begin(); ii < locs.end(); ++ii )
    {
        for ( BigLocus* jj = ii+1; jj < locs.end(); ++jj )
            if ( not_adjacent(ii, jj) )
                checkLL(meca, stiff, *ii, *jj);
    }
}


/**
 This will consider once all pairs of objects from the given lists,
 assuming that the list are different and no object is repeated
 */
void LocusGrid::setSterics0(Meca& meca, real stiff,
                            BigPointList & pots1, BigLocusList & locs1,
                            BigPointList & pots2, BigLocusList & locs2)
{
    assert_true( &pots1 != &pots2 );
    assert_true( &locs1 != &locs2 );
    
    for ( BigPoint* ii = pots1.begin(); ii < pots1.end(); ++ii )
    {
        for ( BigPoint* jj = pots2.begin(); jj < pots2.end(); ++jj )
            if ( not_adjacent(ii, jj) )
                checkPP(meca, stiff, *ii, *jj);
        
        for ( BigLocus* kk = locs2.begin(); kk < locs2.end(); ++kk )
            if ( not_adjacent(ii, kk) )
                checkPL(meca, stiff, *ii, *kk);
    }
    
    for ( BigLocus* ii = locs1.begin(); ii < locs1.end(); ++ii )
    {
        for ( BigPoint* jj = pots2.begin(); jj < pots2.end(); ++jj )
            if ( not_adjacent(jj, ii) )
                checkPL(meca, stiff, *jj, *ii);
        
        for ( BigLocus* kk = locs2.begin(); kk < locs2.end(); ++kk )
        {
            if ( not_adjacent(ii, kk)  )
                checkLL(meca, stiff, *ii, *kk);
        }
    }
}


/**
 This will consider once all pairs of objects from the given lists.
 Compared to `setSterics0()`, this performs additional tests to exclude
 objects that are too far appart to interact, based on BigVector::near()
 
 `rad` is the maximum radius of Points, and the `seg` the maximum range of Locuses
 */
void LocusGrid::setStericsT(Meca& meca, real stiff,
                            BigPointList & pots, BigLocusList & locs)
{
    for ( BigPoint* ii = pots.begin(); ii < pots.end(); ++ii )
    {
        BigVector pos = ii->pos_;
        
        for ( BigPoint* jj = ii+1; jj < pots.end(); ++jj )
            if ( not_adjacent(ii, jj) && pos.near(jj->pos_) )
                checkPP(meca, stiff, *ii, *jj);
        
        for ( BigLocus* kk = locs.begin(); kk < locs.end(); ++kk )
            if ( not_adjacent(ii, kk) && pos.near(kk->pos_) )
                checkPL(meca, stiff, *ii, *kk);
    }

    for ( BigLocus* ii = locs.begin(); ii < locs.end(); ++ii )
    {
        BigVector pos = ii->pos_;
        
        for ( BigLocus* jj = ii+1; jj < locs.end(); ++jj )
            if ( not_adjacent(ii, jj) && pos.near(jj->pos_) )
                checkLL(meca, stiff, *ii, *jj);
    }
}


/**
 This will consider once all pairs of objects from the given lists,
 assuming that the list are different and no object is repeated.

 Compared to `setSterics0()`, this performs additional tests to exclude
 objects that are too far appart to interact, based on BigVector::near()

 `rad` is the maximum radius of Points, and the `seg` the maximum range of Locuses
*/
void LocusGrid::setStericsT(Meca& meca, real stiff,
                            BigPointList & pots1, BigLocusList & locs1,
                            BigPointList & pots2, BigLocusList & locs2)
{
    assert_true( &pots1 != &pots2 );
    assert_true( &locs1 != &locs2 );

    for ( BigPoint* ii = pots1.begin(); ii < pots1.end(); ++ii )
    {
        const BigVector pos = ii->pos_;

        for ( BigPoint* jj = pots2.begin(); jj < pots2.end(); ++jj )
            if ( not_adjacent(ii, jj) && pos.near(jj->pos_) )
                checkPP(meca, stiff, *ii, *jj);
        
        for ( BigLocus* kk = locs2.begin(); kk < locs2.end(); ++kk )
            if ( not_adjacent(ii, kk) && pos.near(kk->pos_) )
                checkPL(meca, stiff, *ii, *kk);
    }
    
    for ( BigLocus* ii = locs1.begin(); ii < locs1.end(); ++ii )
    {
        const BigVector pos = ii->pos_;

        for ( BigPoint* jj = pots2.begin(); jj < pots2.end(); ++jj )
            if ( not_adjacent(jj, ii) && pos.near(jj->pos_) )
                checkPL(meca, stiff, *jj, *ii);
        
        for ( BigLocus* kk = locs2.begin(); kk < locs2.end(); ++kk )
        {
            if ( not_adjacent(ii, kk) && pos.near(kk->pos_) )
                checkLL(meca, stiff, *ii, *kk);
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - Check all pairs of Cells

#if ( MAX_STERIC_PANES == 1 )

/**
 Check interactions between objects contained in the grid.
 */
void LocusGrid::setInteractions(Meca& meca, real stiff) const
{
    //std::clog << "----" << '\n';
    
    // scan all cells to examine each pair of particles:
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int * region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
        
        // We consider each pair of objects (ii, jj) only once:
        
        BigPointList & baseP = point_list(inx);
        BigLocusList & baseL = locus_list(inx);
        
        if ( isPeriodic() )
        {
            setSterics0(meca, stiff, baseP, baseL);
            
            for ( int reg = 1; reg < nr; ++reg )
            {
                BigPointList & sideP = point_list(inx+region[reg]);
                BigLocusList & sideL = locus_list(inx+region[reg]);
                setSterics0(meca, stiff, baseP, baseL, sideP, sideL);
            }
        }
        else
        {
            setStericsT(meca, stiff, baseP, baseL);
            
            for ( int reg = 1; reg < nr; ++reg )
            {
                BigPointList & sideP = point_list(inx+region[reg]);
                BigLocusList & sideL = locus_list(inx+region[reg]);
                setStericsT(meca, stiff, baseP, baseL, sideP, sideL);
            }
        }
    }
}


#else

/**
 Check interactions between the FatPoints contained in Pane `pan`.
 */
void LocusGrid::setInteractions(Meca& meca, real stiff,
                                const size_t pan) const
{
    // scan all cells to examine each pair of particles:
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int * region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
        
        // We consider each pair of objects (ii, jj) only once:
        
        BigPointList & baseP = point_list(inx, pan);
        BigLocusList & baseL = locus_list(inx, pan);
        
        if ( isPeriodic() )
        {
            setSterics0(meca, stiff, baseP, baseL);
            
            for ( int reg = 1; reg < nr; ++reg )
            {
                BigPointList & sideP = point_list(inx+region[reg], pan);
                BigLocusList & sideL = locus_list(inx+region[reg], pan);
                setSterics0(meca, stiff, baseP, baseL, sideP, sideL);
            }
        }
        else
        {
            setStericsT(meca, stiff, baseP, baseL);
            
            for ( int reg = 1; reg < nr; ++reg )
            {
                BigPointList & sideP = point_list(inx+region[reg], pan);
                BigLocusList & sideL = locus_list(inx+region[reg], pan);
                setStericsT(meca, stiff, baseP, baseL, sideP, sideL);
            }
        }
    }
}


/**
 Check interactions between the FatPoints contained in Panes `pan1` and `pan2`,
 where ( pan1 != pan2 )
 */
void LocusGrid::setInteractions(Meca& meca, real stiff,
                                const size_t pan1, const size_t pan2) const
{
    assert_true(pan1 != pan2);
    
    // scan all cells to examine each pair of particles:
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int * region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
        
        // We consider each pair of objects (ii, jj) only once:
        
        BigPointList & baseP = point_list(inx, pan1);
        BigLocusList & baseL = locus_list(inx, pan1);
        
        if ( isPeriodic() )
        {
            for ( int reg = 0; reg < nr; ++reg )
            {
                BigPointList & sideP = point_list(inx+region[reg], pan2);
                BigLocusList & sideL = locus_list(inx+region[reg], pan2);
                setSterics0(meca, stiff, baseP, baseL, sideP, sideL);
            }

            BigPointList & baseP2 = point_list(inx, pan2);
            BigLocusList & baseL2 = locus_list(inx, pan2);
            
            for ( int reg = 1; reg < nr; ++reg )
            {
                BigPointList & sideP = point_list(inx+region[reg], pan1);
                BigLocusList & sideL = locus_list(inx+region[reg], pan1);
                setSterics0(meca, stiff, baseP2, baseL2, sideP, sideL);
            }
        }
        else
        {
            for ( int reg = 0; reg < nr; ++reg )
            {
                BigPointList & sideP = point_list(inx+region[reg], pan2);
                BigLocusList & sideL = locus_list(inx+region[reg], pan2);
                setStericsT(meca, stiff, baseP, baseL, sideP, sideL);
            }

            BigPointList & baseP2 = point_list(inx, pan2);
            BigLocusList & baseL2 = locus_list(inx, pan2);
            
            for ( int reg = 1; reg < nr; ++reg )
            {
                BigPointList & sideP = point_list(inx+region[reg], pan1);
                BigLocusList & sideL = locus_list(inx+region[reg], pan1);
                setStericsT(meca, stiff, baseP2, baseL2, sideP, sideL);
            }
        }
    }
}


#endif

//------------------------------------------------------------------------------
#pragma mark - Display

#ifdef DISPLAY

#include "opengl.h"
void drawBoundaries(Map<DIM> const&);

void LocusGrid::draw() const
{
#if ( DIM <= 3 )
    glDisable(GL_LIGHTING);
    glColor3f(1, 0, 0);
    glLineWidth(0.25);
    drawBoundaries(pGrid);
#endif
}
#endif

