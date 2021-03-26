// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "assert_macro.h"
#include "point_grid.h"
#include "exceptions.h"
#include "messages.h"
#include "modulo.h"
#include "space.h"
#include "meca.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------

PointGrid::PointGrid()
: max_diameter(0)
{
}


size_t PointGrid::setGrid(Space const* spc, real min_step)
{
    assert_true( min_step > 0 );
    Vector inf, sup;
    spc->boundaries(inf, sup);
    
    size_t cnt[3] = { 1, 1, 1 };
    for ( size_t d = 0; d < DIM; ++d )
    {
        real n = ( sup[d] - inf[d] ) / min_step;
        
        if ( n < 0 )
            throw InvalidParameter("invalid space:boundaries");
        
        if ( modulo  &&  modulo->isPeriodic(d) )
        {
            //adjust the grid to match the edges
            cnt[d] = std::max((size_t)1, (size_t)std::floor(n));
            pGrid.setPeriodic(d, true);
        }
        else
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


void PointGrid::createCells()
{
    pGrid.createCells();
    
    //Create side regions suitable for pairwise interactions:
    pGrid.createSideRegions(1);
    
    //The maximum allowed diameter of particles is half the minimum cell width
    max_diameter = pGrid.minimumWidth(1);
    
    //report the grid size used
    pGrid.printSummary(Cytosim::log, "StericGrid");
}


//------------------------------------------------------------------------------
#pragma mark -

/// include verifications that the grid range is appropriate
#define CHECK_STERIC_RANGE 0


#if ( N_STERIC_PANES != 1 )

void PointGrid::add(size_t pan, Mecable const* mec, size_t inx, real rad, real rg) const
{
    if ( pan == 0 || pan > N_STERIC_PANES )
        throw InvalidParameter("point:steric is out-of-range");
    
    Vector w = mec->posPoint(inx);
    point_list(w, pan).emplace(Mecapoint(mec, inx), rad, rg, w);
    
#if ( CHECK_STERIC_RANGE )
    //we check that the grid would correctly detect collision of two particles
    if ( max_diameter < 1.999 * rg )
    {
        InvalidParameter e("simul:steric_max_range is too short");
        e << PREF << "steric_max_range should be greater than 2 * ( particle_radius + extra_range )\n";
        e << PREF << "= " << 2 * rg << " for some particles\n";
        throw e;
    }
#endif
}


void PointGrid::add(size_t pan, Fiber const* fib, size_t inx, real rd, real rg) const
{
    if ( pan == 0 || pan > N_STERIC_PANES )
        throw InvalidParameter("line:steric is out-of-range");
    
    // link in the cell containing the middle of the segment:
    Vector w = fib->posPoint(inx, 0.5);
    locus_list(w, pan).emplace(FiberSegment(fib, inx), rd, rg, w);
    
#if ( CHECK_STERIC_RANGE )
    //we check that the grid would correctly detect collision of two segments
    //along the diagonal, corresponding to the worst-case scenario
    real diag = square(fib->segmentation()) + square(2*rg);
    if ( square(max_diameter) * 1.001 < diag )
    {
        InvalidParameter e("simul:steric_max_range is too short");
        e << PREF << "steric_max_range should be greater than std::sqrt( sqr(segment_length) + 4*sqr(range) )\n";
        e << PREF << "with segment_length ~ 4/3 segmentation\n";
        e << PREF << "= " << diag << " for some fibers\n";
        throw e;
    }
#endif
}


#endif


//------------------------------------------------------------------------------
#pragma mark - Check two Objects: P = Point; L = Line segment


/**
 This is used to check two spherical objects:
 Solid/Bead/Sphere or the terminal vertex (the tips) of a Fiber
 
 The force is applied if the objects are closer than the
 sum of their radiuses.
 */
void PointGrid::checkPP(Meca& meca, StericParam const& pam,
                        FatPoint const& aa, FatPoint const& bb)
{
    //std::clog << "   PP- " << bb.pnt << " " << aa.pnt << '\n';
    const real len = aa.radius + bb.radius;
    Vector vab = bb.pos - aa.pos;
    
    if ( modulo )
        modulo->fold(vab);
    
    if ( vab.normSqr() < len*len )
        meca.addLongLink(aa.pnt, bb.pnt, len, pam.stiff_push);
}


/**
 This is used to check a segment of a fiber against a spherical object:
 Solid/Bead/Sphere or Fiber-tip.
 
 The force is applied if the objects are closer than the sum of their radiuses.
 */
void PointGrid::checkPL(Meca& meca, StericParam const& pam,
                        FatPoint const& aa, FatLocus const& bb)
{
    //std::clog << "   PL- " << bb.seg << " " << aa.pnt << '\n';
    const real len = aa.radius + bb.radius;
    
    // get position of point with respect to segment:
    real dis2 = INFINITY;
    real abs = bb.seg.projectPoint0(aa.pos, dis2);
    
    if ( 0 <= abs )
    {
        if ( abs <= bb.seg.len() )
        {
            if ( dis2 < len*len )
                meca.addSideSlidingLink(Interpolation(bb.seg, abs), aa.pnt, len, pam.stiff_push);
        }
        else
        {
            if ( bb.isLast() )
                checkPP(meca, pam, aa, bb.fatPoint2());
        }
    }
    else
    {
        if ( bb.isFirst() )
            checkPP(meca, pam, aa, bb.fatPoint1());
        else
        {
            /* we check the projection to the previous segment,
             and if it falls on the right of it, then we interact with the node */
            Vector vab = aa.pos - bb.seg.pos1();
            
            if ( modulo )
                modulo->fold(vab);
            
            if ( dot(vab, bb.seg.fiber()->diffPoints(bb.seg.point()-1)) >= 0 )
            {
                if ( vab.normSqr() < len*len )
                    meca.addLongLink(aa.pnt, bb.seg.exact1(), len, pam.stiff_push);
            }
        }
    }
}


/**
 This is used to check a segment of a fiber against another segment of fiber,
 not including the terminal vertex of fibers.
 
 The interaction is applied only if the vertex projects 'inside' the segment.
 */
void PointGrid::checkLL1(Meca& meca, StericParam const& pam,
                         FatLocus const& aa, FatLocus const& bb)
{
    //std::clog << "   LL1 " << aa.seg << " " << bb.point1() << '\n';
    const real ran = aa.range + bb.radius;
    
    // get position of bb.point1() with respect to segment 'aa'
    real dis2 = INFINITY;
    real abs = aa.seg.projectPoint0(bb.seg.pos1(), dis2);
    
    if ((0 <= abs) & (abs <= aa.seg.len()) & (dis2 < ran*ran))
    {
        /*
         bb.point1() projects inside segment 'aa'
         */
        const real len = aa.radius + bb.radius;
        real stiff = sign_select(dis2-len*len, pam.stiff_push, pam.stiff_pull);
        meca.addSideSlidingLink(Interpolation(aa.seg, abs), bb.seg.exact1(), len, stiff);
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
                Vector vab = bb.seg.pos1() - aa.seg.pos1();
                
                if ( modulo )
                    modulo->fold(vab);
                
                const real len = aa.radius + bb.radius;
                if ( vab.normSqr() < len*len  &&  dot(vab, bb.seg.diff()) >= 0 )
                    meca.addLongLink(aa.seg.exact1(), bb.seg.exact1(), len, pam.stiff_push);
            }
        }
        else
        {
            /*
             Check the projection to the segment located before 'aa',
             and interact if 'bb.point1()' falls on the right side of it
             */
            Vector vab = bb.seg.pos1() - aa.seg.pos1();
            
            if ( modulo )
                modulo->fold(vab);
            
            if ( dot(vab, aa.seg.fiber()->diffPoints(aa.seg.point()-1)) >= 0 )
            {
                const real d = vab.normSqr();
                if ( d < ran*ran )
                {
                    const real len = aa.radius + bb.radius;
                    real stiff = sign_select(d-len*len, pam.stiff_push, pam.stiff_pull);
                    meca.addLongLink(aa.seg.exact1(), bb.seg.exact1(), len, stiff);
                }
            }
        }
    }
}


/**
 This is used to check a segment of a fiber against the terminal vertex of fiber
 
 The interaction is applied only if the vertex projects 'inside' the segment.
 */
void PointGrid::checkLL2(Meca& meca, StericParam const& pam,
                         FatLocus const& aa, FatLocus const& bb)
{
    //std::clog << "   LL2 " << aa.seg << " " << bb.point2() << '\n';
    const real ran = aa.range + bb.radius;
    
    // get position of bb.point2() with respect to segment 'aa'
    real dis2 = INFINITY;
    real abs = aa.seg.projectPoint0(bb.seg.pos2(), dis2);
    
    if ((0 <= abs) & (dis2 < ran*ran) & (abs <= aa.seg.len()))
    {
        /*
         bb.point2() projects inside segment 'aa'
         */
        const real len = aa.radius + bb.radius;
        real stiff = sign_select(dis2-len*len, pam.stiff_push, pam.stiff_pull);
        meca.addSideSlidingLink(Interpolation(aa.seg, abs), bb.seg.exact2(), len, stiff);
    }
    else if ( abs < 0 )
    {
        /*
         Check the projection to the segment located before 'aa',
         and interact if 'bb.point1()' falls on the right side of it
         */
        Vector vab = bb.seg.pos2() - aa.seg.pos1();
        
        if ( modulo )
            modulo->fold(vab);
        
        if ( aa.isFirst() )
        {
            assert_true(bb.isLast());
            const real len = aa.radius + bb.radius;
            if ( vab.normSqr() < len*len  && dot(vab, bb.seg.diff()) <= 0 )
                meca.addLongLink(aa.seg.exact1(), bb.seg.exact2(), len, pam.stiff_push);
        }
        else
        {
            if ( dot(vab, aa.seg.fiber()->diffPoints(aa.seg.point()-1)) >= 0 )
            {
                const real d = vab.normSqr();
                if ( d < ran*ran )
                {
                    const real len = aa.radius + bb.radius;
                    real stiff = sign_select(d-len*len, pam.stiff_push, pam.stiff_pull);
                    meca.addLongLink(aa.seg.exact1(), bb.seg.exact2(), len, stiff);
                }
            }
        }
    }
    else if ( &bb < &aa  &&  aa.isLast()  &&  abs > aa.seg.len() )
    {
        /*
         Check the projection of aa.point2(),
         on the segment represented by 'bb'
         */
        assert_true(bb.isLast());
        
        Vector vab = bb.seg.pos2() - aa.seg.pos2();
        
        if ( modulo )
            modulo->fold(vab);
        
        const real len = aa.radius + bb.radius;
        if ( vab.normSqr() < len*len  &&  dot(vab, bb.seg.diff()) <= 0 )
            meca.addLongLink(aa.seg.exact2(), bb.seg.exact2(), len, pam.stiff_push);
    }
}


/**
 This is used to check two FiberSegment, each representing the segment of a Fiber.
 in 3D, the segments are tested for getting within the requested distance.
 in 2D, only the vertices are checked.
 */
void PointGrid::checkLL(Meca& meca, StericParam const& pam,
                        FatLocus const& aa, FatLocus const& bb)
{
#if ( DIM == 3 )
    
    const real ran = std::max(aa.range+bb.radius, aa.radius+bb.range);

    /* in 3D, we use shortestDistance() to calculate the closest distance
     between two segments, and use the result to build an interaction */
    real a, b;
    real d = aa.seg.shortestDistance(bb.seg, a, b);
    if ( d >= ran*ran )
        return;
    
    if ( aa.seg.within(a) & bb.seg.within(b) )
    {
        const real len = aa.radius + bb.radius;
        real stiff = sign_select(d-len*len, pam.stiff_push, pam.stiff_pull);
        meca.addSideSlidingLink(Interpolation(aa.seg, a), Interpolation(bb.seg, b), ran, stiff);
    }
#endif
    
    //std::clog << "LL " << aa.seg << " " << bb.seg << '\n';
    checkLL1(meca, pam, aa, bb);
    
    if ( aa.isLast() )
        checkLL2(meca, pam, bb, aa);
    
    checkLL1(meca, pam, bb, aa);
    
    if ( bb.isLast() )
        checkLL2(meca, pam, aa, bb);
}


//------------------------------------------------------------------------------
#pragma mark - Selections of pairs excluded from Sterics


/// excluding two spheres when they are from the same Solid
inline bool adjacent(FatPoint const* a, FatPoint const* b)
{
    return ( a->pnt.mecable() == b->pnt.mecable() );
}


/// excluding Fiber and Solid from the same Aster
inline bool adjacent(FatPoint const* a, FatLocus const* b)
{
    //a->pnt.mecable()->Buddy::print(std::clog);
    //b->seg.fiber()->Buddy::print(std::clog);
    return b->seg.fiber()->buddy() == a->pnt.mecable()->buddy();
}


/// excluding segments that are adjacent on the same fiber, or protofilaments from Tubule
inline bool adjacent(FatLocus const* a, FatLocus const* b)
{
#if FIBER_HAS_FAMILY
    return (( a->seg.fiber()->family_ == b->seg.fiber()->family_ )
#else
    return (( a->seg.fiber() == b->seg.fiber() )
#endif
            & ( a->seg.point() < 2 + b->seg.point() ) & ( b->seg.point() < 2 + a->seg.point() ));
}

//------------------------------------------------------------------------------
#pragma mark - Check all possible object pairs from two Cells

/**
 This will consider once all pairs of objects from the given lists
 */
void PointGrid::setInteractions(Meca& meca, StericParam const& stiff,
                                FatPointList & pots, FatLocusList & locs)
{
    for ( FatPoint* ii = pots.begin(); ii < pots.end(); ++ii )
    {
        for ( FatPoint* jj = ii+1; jj < pots.end(); ++jj )
            if ( !adjacent(ii, jj) )
                checkPP(meca, stiff, *ii, *jj);
        
        for ( FatLocus* kk = locs.begin(); kk < locs.end(); ++kk )
            if ( !adjacent(ii, kk) )
                checkPL(meca, stiff, *ii, *kk);
    }

    for ( FatLocus* ii = locs.begin(); ii < locs.end(); ++ii )
    {
        for ( FatLocus* jj = ii+1; jj < locs.end(); ++jj )
            if ( !adjacent(ii, jj) )
                checkLL(meca, stiff, *ii, *jj);
    }
}


/**
 This will consider once all pairs of objects from the given lists,
 assuming that the list are different and no object is repeated
 */
void PointGrid::setInteractions(Meca& meca, StericParam const& pam,
                                FatPointList & pots1, FatLocusList & locs1,
                                FatPointList & pots2, FatLocusList & locs2)
{
    assert_true( &pots1 != &pots2 );
    assert_true( &locs1 != &locs2 );
    
    for ( FatPoint* ii = pots1.begin(); ii < pots1.end(); ++ii )
    {
        for ( FatPoint* jj = pots2.begin(); jj < pots2.end(); ++jj )
            if ( !adjacent(ii, jj) )
                checkPP(meca, pam, *ii, *jj);
        
        for ( FatLocus* kk = locs2.begin(); kk < locs2.end(); ++kk )
            if ( !adjacent(ii, kk) )
                checkPL(meca, pam, *ii, *kk);
    }
    
    for ( FatLocus* ii = locs1.begin(); ii < locs1.end(); ++ii )
    {
        for ( FatPoint* jj = pots2.begin(); jj < pots2.end(); ++jj )
            if ( !adjacent(jj, ii) )
                checkPL(meca, pam, *jj, *ii);
        
        for ( FatLocus* kk = locs2.begin(); kk < locs2.end(); ++kk )
        {
            if ( !adjacent(ii, kk)  )
                checkLL(meca, pam, *ii, *kk);
        }
    }
}


/**
 This will consider once all pairs of objects from the given lists
 Compared to `setInteractions()`, this performs an additional test
 for the distance between the object's `pos` is below `max_diameter`
 */
void PointGrid::setInteractions(Meca& meca, StericParam const& pam, real sup,
                                FatPointList & pots, FatLocusList & locs)
{
    for ( FatPoint* ii = pots.begin(); ii < pots.end(); ++ii )
    {
        Vector pos = ii->pos;
        for ( FatPoint* jj = ii+1; jj < pots.end(); ++jj )
            if ( !adjacent(ii, jj) && distanceSqr(pos, jj->pos) <= sup )
                checkPP(meca, pam, *ii, *jj);
        
        for ( FatLocus* kk = locs.begin(); kk < locs.end(); ++kk )
            if ( !adjacent(ii, kk) && distanceSqr(pos, kk->pos) <= sup )
                checkPL(meca, pam, *ii, *kk);
    }

    for ( FatLocus* ii = locs.begin(); ii < locs.end(); ++ii )
    {
        Vector pos = ii->pos;
        for ( FatLocus* jj = ii+1; jj < locs.end(); ++jj )
            if ( !adjacent(ii, jj) && distanceSqr(pos, jj->pos) <= sup )
                checkLL(meca, pam, *ii, *jj);
    }
}


/**
 This will consider once all pairs of objects from the given lists,
 assuming that the list are different and no object is repeated.

 Compared to `setInteractions()`, this performs an additional test
 for the distance between the object's `pos` is below `max_diameter`
 */
void PointGrid::setInteractions(Meca& meca, StericParam const& pam, real sup,
                                FatPointList & pots1, FatLocusList & locs1,
                                FatPointList & pots2, FatLocusList & locs2)
{
    assert_true( &pots1 != &pots2 );
    assert_true( &locs1 != &locs2 );

    for ( FatPoint* ii = pots1.begin(); ii < pots1.end(); ++ii )
    {
        const Vector pos = ii->pos;

        for ( FatPoint* jj = pots2.begin(); jj < pots2.end(); ++jj )
            if ( !adjacent(ii, jj) && distanceSqr(pos, jj->pos) <= sup )
                checkPP(meca, pam, *ii, *jj);
        
        for ( FatLocus* kk = locs2.begin(); kk < locs2.end(); ++kk )
            if ( !adjacent(ii, kk) && distanceSqr(pos, kk->pos) <= sup )
                checkPL(meca, pam, *ii, *kk);
    }
    
    for ( FatLocus* ii = locs1.begin(); ii < locs1.end(); ++ii )
    {
        const Vector pos = ii->pos;

        for ( FatPoint* jj = pots2.begin(); jj < pots2.end(); ++jj )
            if ( !adjacent(jj, ii) && distanceSqr(pos, jj->pos) <= sup )
                checkPL(meca, pam, *jj, *ii);
        
        for ( FatLocus* kk = locs2.begin(); kk < locs2.end(); ++kk )
        {
            if ( !adjacent(ii, kk) && distanceSqr(pos, kk->pos) <= sup )
                checkLL(meca, pam, *ii, *kk);
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - Check all pairs of Cells


#if ( N_STERIC_PANES == 1 )

/**
 Check interactions between objects contained in the grid.
 */
void PointGrid::setInteractions(Meca& meca, StericParam const& pam) const
{
    //std::clog << "----" << '\n';
    
    // scan all cells to examine each pair of particles:
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int * region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
        
        // We consider each pair of objects (ii, jj) only once:
        
        FatPointList & baseP = point_list(inx);
        FatLocusList & baseL = locus_list(inx);
        
        if ( isPeriodic() )
        {
            setInteractions(meca, pam, baseP, baseL);
            
            for ( int reg = 1; reg < nr; ++reg )
            {
                FatPointList & sideP = point_list(inx+region[reg]);
                FatLocusList & sideL = locus_list(inx+region[reg]);
                setInteractions(meca, pam, baseP, baseL, sideP, sideL);
            }
        }
        else
        {
            const real sup = square(max_diameter);
            setInteractions(meca, pam, sup, baseP, baseL);
            
            for ( int reg = 1; reg < nr; ++reg )
            {
                FatPointList & sideP = point_list(inx+region[reg]);
                FatLocusList & sideL = locus_list(inx+region[reg]);
                setInteractions(meca, pam, sup, baseP, baseL, sideP, sideL);
            }
        }
    }
}


#else

/**
 Check interactions between the FatPoints contained in Pane `pan`.
 */
void PointGrid::setInteractions(Meca& meca, StericParam const& pam,
                                const size_t pan) const
{
    // scan all cells to examine each pair of particles:
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int * region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
        
        // We consider each pair of objects (ii, jj) only once:
        
        FatPointList & baseP = point_list(inx, pan);
        FatLocusList & baseL = locus_list(inx, pan);
        
        if ( isPeriodic() )
        {
            setInteractions(meca, pam, baseP, baseL);
            
            for ( int reg = 1; reg < nr; ++reg )
            {
                FatPointList & sideP = point_list(inx+region[reg], pan);
                FatLocusList & sideL = locus_list(inx+region[reg], pan);
                setInteractions(meca, pam, baseP, baseL, sideP, sideL);
            }
        }
        else
        {
            const real sup = square(max_diameter);
            setInteractions(meca, pam, sup, baseP, baseL);
            
            for ( int reg = 1; reg < nr; ++reg )
            {
                FatPointList & sideP = point_list(inx+region[reg], pan);
                FatLocusList & sideL = locus_list(inx+region[reg], pan);
                setInteractions(meca, pam, sup, baseP, baseL, sideP, sideL);
            }
        }
    }
}


/**
 Check interactions between the FatPoints contained in Panes `pan1` and `pan2`,
 where ( pan1 != pan2 )
 */
void PointGrid::setInteractions(Meca& meca, StericParam const& pam,
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
        
        FatPointList & baseP = point_list(inx, pan1);
        FatLocusList & baseL = locus_list(inx, pan1);
        
        if ( isPeriodic() )
        {
            for ( int reg = 0; reg < nr; ++reg )
            {
                FatPointList & sideP = point_list(inx+region[reg], pan2);
                FatLocusList & sideL = locus_list(inx+region[reg], pan2);
                setInteractions(meca, pam, baseP, baseL, sideP, sideL);
            }
            
            FatPointList & baseP2 = point_list(inx, pan2);
            FatLocusList & baseL2 = locus_list(inx, pan2);
            
            for ( int reg = 1; reg < nr; ++reg )
            {
                FatPointList & sideP = point_list(inx+region[reg], pan1);
                FatLocusList & sideL = locus_list(inx+region[reg], pan1);
                setInteractions(meca, pam, baseP2, baseL2, sideP, sideL);
            }
        }
        else
        {
            const real sup = square(max_diameter);
            for ( int reg = 0; reg < nr; ++reg )
            {
                BigPointList & sideP = point_list(inx+region[reg], pan2);
                BigLocusList & sideL = locus_list(inx+region[reg], pan2);
                setInteractions(meca, pam, sup, baseP, baseL, sideP, sideL);
            }

            BigPointList & baseP2 = point_list(inx, pan2);
            BigLocusList & baseL2 = locus_list(inx, pan2);
            
            for ( int reg = 1; reg < nr; ++reg )
            {
                BigPointList & sideP = point_list(inx+region[reg], pan1);
                BigLocusList & sideL = locus_list(inx+region[reg], pan1);
                setInteractions(meca, pam, sup, baseP2, baseL2, sideP, sideL);
            }
        }

    }
}


#endif


//------------------------------------------------------------------------------
#pragma mark - Display

#ifdef DISPLAY

#include "opengl.h"
void drawEdges(Map<DIM> const&);

void PointGrid::draw() const
{
    glDisable(GL_LIGHTING);
    glColor3f(0, 1, 0);
    glLineWidth(0.25);
    drawEdges(pGrid);
}
#endif

