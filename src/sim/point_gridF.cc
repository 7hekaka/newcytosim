// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "assert_macro.h"
#include "point_gridF.h"
#include "exceptions.h"
#include "messages.h"
#include "modulo.h"
#include "space.h"
#include "meca.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------

PointGridF::PointGridF()
: max_diameter(0)
{
}


size_t PointGridF::setGrid(Space const* spc, real min_step)
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
        
#if GRID_HAS_PERIODIC
        if ( modulo  &&  modulo->isPeriodic(d) )
        {
            //adjust the grid to match the edges
            cnt[d] = std::max((size_t)1, (size_t)floor(n));
            pGrid.setPeriodic(d, true);
        }
        else
#endif
        {
            //add a border in any dimension which is not periodic
            cnt[d] = (size_t)ceil(n) + 2;
            n = cnt[d] * 0.5 * min_step;
            real mid = inf[d] + sup[d];
            inf[d] = mid - n;
            sup[d] = mid + n;
        }
    }
    
    pGrid.setDimensions(inf, sup, cnt);
    return pGrid.nbCells();
}

void PointGridF::createCells()
{
    pGrid.createCells();
    
    //Create side regions suitable for pairwise interactions:
    pGrid.createSideRegions(1);
    
    //The maximum allowed diameter of particles is half the minimum cell width
    max_diameter = pGrid.minimumWidth(1);
    
    //report the grid size used
    if ( pGrid.nbCells() > 4096 )
        pGrid.printSummary(std::clog, "StericGrid");
}


//------------------------------------------------------------------------------
#pragma mark -

/// include verifications that the grid range is appropriate
#define CHECK_STERIC_RANGE 0


#if ( MAX_STERIC_PANES != 1 )

void PointGridF::add(size_t pan, Mecapoint const& pe, real rd) const
{
    if ( pan == 0 || pan > MAX_STERIC_PANES )
        throw InvalidParameter("point:steric is out-of-range");
    
    Vector w = pe.pos();
    point_list(w, pan).emplace_back(pe, rd, w);
    
#if ( CHECK_STERIC_RANGE )
    //we check that the grid would correctly detect collision of two particles
    if ( max_diameter < 1.999 * rd )
    {
        InvalidParameter e("simul:steric_max_range is too short");
        e << PREF << "steric_max_range should be greater than 2 * ( particle_radius + extra_range )\n";
        e << PREF << "= " << 2 * rd << " for some particles\n";
        throw e;
    }
#endif
}


void PointGridF::add(size_t pan, FiberSegment const& fl, real rd) const
{
    if ( pan == 0 || pan > MAX_STERIC_PANES )
        throw InvalidParameter("line:steric is out-of-range");
    
    // link in the cell containing the middle of the segment:
    Vector w = fl.center();
    locus_list(w, pan).emplace_back(fl, rd, w);
    
#if ( CHECK_STERIC_RANGE )
    //we check that the grid would correctly detect collision of two segments
    //along the diagonal, corresponding to the worst-case scenario
    real diag = square(fl.len()) + square(2*rd);
    if ( square(max_diameter) * 1.001 < diag )
    {
        InvalidParameter e("simul:steric_max_range is too short");
        e << PREF << "steric_max_range should be greater than sqrt( sqr(segment_length) + 4*sqr(range) )\n";
        e << PREF << "with segment_length ~ 4/3 segmentation\n";
        e << PREF << "= " << diag << " for some fibers\n";
        throw e;
    }
#endif
}


#endif


//------------------------------------------------------------------------------
#pragma mark - Steric functions


/**
 This is used to check two spherical objects:
 Solid/Bead/Sphere or the terminal vertex (the tips) of a Fiber
 
 The force is applied if the objects are closer to the maximum
 of their specified range + radius.
 */
void PointGridF::checkPP(Meca& meca, real stiff,
                         FatPointF const& aa, FatPointF const& bb) const
{
    //std::clog << "   PP- " << bb.pnt << " " << aa.pnt << std::endl;
    const real ran = aa.radius + bb.radius;
    Vector vab = bb.pos - aa.pos;
    
#if GRID_HAS_PERIODIC
    if ( modulo )
        modulo->fold(vab);
#endif
    if ( vab.normSqr() < ran*ran )
        meca.addLongLink(aa.pnt, bb.pnt, ran, stiff);
}


/**
 This is used to check a segment of a fiber against a spherical object:
 Solid/Bead/Sphere or Fiber-tip.
 
 The force is applied if the objects are closer than the maximum of the two range + radius,
 and if the center of the sphere projects inside the segment.
 */
void PointGridF::checkPL(Meca& meca, real stiff,
                         FatPointF const& aa, FatLocusF const& bb) const
{
    //std::clog << "   PL- " << bb.seg << " " << aa.pnt << std::endl;
    const real ran = aa.radius + bb.radius;
    
    // get position of point with respect to segment:
    real dis2 = INFINITY;
    real abs = bb.seg.projectPoint1(aa.pos, dis2);
    
    if ( 0 <= abs )
    {
        if ( abs <= bb.seg.len() )
        {
            if ( dis2 < ran*ran )
                meca.addSideSlidingLink(Interpolation(bb.seg, abs), aa.pnt, ran, stiff);
        }
        else
        {
            if ( bb.isLast() )
                checkPP(meca, stiff, aa, bb.point2());
        }
    }
    else
    {
        if ( bb.isFirst() )
            checkPP(meca, stiff, aa, bb.point1());
        else
        {
            /* we check the projection to the previous segment,
             and if it falls on the right of it, then we interact with the node */
            Vector vab = aa.pos - bb.seg.pos1();
            
#if GRID_HAS_PERIODIC
            if ( modulo )
                modulo->fold(vab);
#endif
            if ( dot(vab, bb.seg.fiber()->diffPoints(bb.seg.point()-1)) >= 0 )
            {
                if ( vab.normSqr() < ran*ran )
                    meca.addLongLink(aa.pnt, bb.seg.exact1(), ran, stiff);
            }
        }
    }
}


/**
 This is used to check a segment of a fiber against another segment of fiber,
 not including the terminal vertex of fibers.
 
 The interaction is applied only if the vertex projects 'inside' the segment.
 */
void PointGridF::checkLL1(Meca& meca, real stiff,
                          FatLocusF const& aa, FatLocusF const& bb) const
{
    //std::clog << "   LL1 " << aa.seg << " " << bb.point1() << std::endl;
    const real ran = aa.radius + bb.radius;
    
    // get position of bb.point1() with respect to segment 'aa'
    real dis2 = INFINITY;
    real abs = aa.seg.projectPoint1(bb.seg.pos1(), dis2);
    
    if ((0 <= abs)  & (abs <= aa.seg.len())  & (dis2 < ran*ran))
    {
        /*
         bb.point1() projects inside segment 'aa'
         */
        meca.addSideSlidingLink(Interpolation(aa.seg, abs), bb.seg.exact1(), ran, stiff);
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
                
#if GRID_HAS_PERIODIC
                if ( modulo )
                    modulo->fold(vab);
#endif
                if ( vab.normSqr() < ran*ran  &&  dot(vab, bb.seg.diff()) >= 0 )
                    meca.addLongLink(aa.seg.exact1(), bb.seg.exact1(), ran, stiff);
            }
        }
        else
        {
            /*
             Check the projection to the segment located before 'aa',
             and interact if 'bb.point1()' falls on the right side of it
             */
            Vector vab = bb.seg.pos1() - aa.seg.pos1();
            
#if GRID_HAS_PERIODIC
            if ( modulo )
                modulo->fold(vab);
#endif
            if ( dot(vab, aa.seg.fiber()->diffPoints(aa.seg.point()-1)) >= 0 )
            {
                if ( vab.normSqr() < ran*ran )
                    meca.addLongLink(aa.seg.exact1(), bb.seg.exact1(), ran, stiff);
            }
        }
    }
}


/**
 This is used to check a segment of a fiber against another segment of fiber,
 the non-terminal vertex of a fiber.
 
 The interaction is applied only if the vertex projects 'inside' the segment.
 */
void PointGridF::checkLL2(Meca& meca, real stiff,
                          FatLocusF const& aa, FatLocusF const& bb) const
{
    //std::clog << "   LL2 " << aa.seg << " " << bb.point2() << std::endl;
    const real ran = aa.radius + bb.radius;
    
    // get position of bb.point2() with respect to segment 'aa'
    real dis2 = INFINITY;
    real abs = aa.seg.projectPoint1(bb.seg.pos2(), dis2);
    
    if ((0 <= abs) & (dis2 < ran*ran) & (abs <= aa.seg.len()))
    {
        /*
         bb.point2() projects inside segment 'aa'
         */
        const real len = aa.radius + bb.radius;
        meca.addSideSlidingLink(Interpolation(aa.seg, abs), bb.seg.exact2(), len, stiff);
    }
    else if ( abs < 0 )
    {
        /*
         Check the projection to the segment located before 'aa',
         and interact if 'bb.point1()' falls on the right side of it
         */
        Vector vab = bb.seg.pos2() - aa.seg.pos1();
        
#if GRID_HAS_PERIODIC
        if ( modulo )
            modulo->fold(vab);
#endif
        if ( aa.isFirst() )
        {
            assert_true(bb.isLast());
            if ( vab.normSqr() < ran*ran  && dot(vab, bb.seg.diff()) <= 0 )
                meca.addLongLink(aa.seg.exact1(), bb.seg.exact2(), ran, stiff);
        }
        else
        {
            if ( dot(vab, aa.seg.fiber()->diffPoints(aa.seg.point()-1)) >= 0 )
            {
                if ( vab.normSqr() < ran*ran )
                    meca.addLongLink(aa.seg.exact1(), bb.seg.exact2(), ran, stiff);
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
        
#if GRID_HAS_PERIODIC
        if ( modulo )
            modulo->fold(vab);
#endif
        if ( vab.normSqr() < ran*ran  &&  dot(vab, bb.seg.diff()) <= 0 )
            meca.addLongLink(aa.seg.exact2(), bb.seg.exact2(), ran, stiff);
    }
}


/**
 This is used to check two FiberSegment, that each represent a segment of a Fiber.
 The segments are tested for intersection in 3D.
 */
void PointGridF::checkLL(Meca& meca, real stiff,
                         FatLocusF const& aa, FatLocusF const& bb) const
{
#if ( DIM == 3 )
    
    const real ran = aa.radius + bb.radius;
    
    /* in 3D, we use shortestDistance() to calculate the closest distance
     between two segments, and use the result to build an interaction */
    real a, b;
    real d = aa.seg.shortestDistance(bb.seg, a, b);
    if ( d >= ran*ran )
        return;
    
    if ( aa.seg.within(a) & bb.seg.within(b) )
        meca.addSideSlidingLink(Interpolation(aa.seg, a), Interpolation(bb.seg, b), ran, stiff);
    
#endif

    //std::clog << "LL " << aa.seg << " " << bb.seg << std::endl;
    checkLL1(meca, stiff, aa, bb);
    
    if ( aa.isLast() )
        checkLL2(meca, stiff, bb, aa);
    
    checkLL1(meca, stiff, bb, aa);
    
    if ( bb.isLast() )
        checkLL2(meca, stiff, aa, bb);
}


//------------------------------------------------------------------------------
#pragma mark -


/**
 This will consider once all pairs of objects from the given lists
 */
void PointGridF::setInteractions(Meca& meca, real stiff,
                                 FatPointListF & fpl, FatLocusListF & fll) const
{
    for ( FatPointF* ii = fpl.begin(); ii < fpl.end(); ++ii )
    {
        for ( FatPointF* jj = ii+1; jj < fpl.end(); ++jj )
            checkPP(meca, stiff, *ii, *jj);
        
        for ( FatLocusF* kk = fll.begin(); kk < fll.end(); ++kk )
            checkPL(meca, stiff, *ii, *kk);
    }

    if ( isPeriodic() )
    {
        for ( FatLocusF* ii = fll.begin(); ii < fll.end(); ++ii )
        {
            for ( FatLocusF* jj = ii+1; jj < fll.end(); ++jj )
                if ( !adjacent(ii->seg, jj->seg) )
                    checkLL(meca, stiff, *ii, *jj);
        }
    }
    else
    {
        const real sup = square(max_diameter);
        for ( FatLocusF* ii = fll.begin(); ii < fll.end(); ++ii )
        {
            Vector pos = ii->pos;
            for ( FatLocusF* jj = ii+1; jj < fll.end(); ++jj )
                if ( !adjacent(ii->seg, jj->seg) && distanceSqr(pos, jj->pos) <= sup )
                    checkLL(meca, stiff, *ii, *jj);
        }
    }
}


/**
 This will consider once all pairs of objects from the given lists,
 assuming that the list are different and no object is repeated
 */
void PointGridF::setInteractions(Meca& meca, real stiff,
                                 FatPointListF & fpl1, FatLocusListF & fll1,
                                 FatPointListF & fpl2, FatLocusListF & fll2) const
{
    assert_true( &fpl1 != &fpl2 );
    assert_true( &fll1 != &fll2 );
    
    for ( FatPointF* ii = fpl1.begin(); ii < fpl1.end(); ++ii )
    {
        for ( FatPointF* jj = fpl2.begin(); jj < fpl2.end(); ++jj )
            checkPP(meca, stiff, *ii, *jj);
        
        for ( FatLocusF* kk = fll2.begin(); kk < fll2.end(); ++kk )
            checkPL(meca, stiff, *ii, *kk);
    }
    
    if ( isPeriodic() )
    {
        for ( FatLocusF* ii = fll1.begin(); ii < fll1.end(); ++ii )
        {
            for ( FatPointF* jj = fpl2.begin(); jj < fpl2.end(); ++jj )
                checkPL(meca, stiff, *jj, *ii);
            
            for ( FatLocusF* kk = fll2.begin(); kk < fll2.end(); ++kk )
            {
                if ( !adjacent(ii->seg, kk->seg)  )
                    checkLL(meca, stiff, *ii, *kk);
            }
        }
    }
    else
    {
        const real sup = square(max_diameter);
        for ( FatLocusF* ii = fll1.begin(); ii < fll1.end(); ++ii )
        {
            for ( FatPointF* jj = fpl2.begin(); jj < fpl2.end(); ++jj )
                checkPL(meca, stiff, *jj, *ii);
            
            Vector pos = ii->pos;
            for ( FatLocusF* kk = fll2.begin(); kk < fll2.end(); ++kk )
            {
                if ( !adjacent(ii->seg, kk->seg) && distanceSqr(pos, kk->pos) <= sup )
                    checkLL(meca, stiff, *ii, *kk);
            }
        }
    }
}


#if ( MAX_STERIC_PANES == 1 )

/**
 Check interactions between objects contained in the grid.
 */
void PointGridF::setInteractions(Meca& meca, real stiff) const
{
    //std::clog << "----" << std::endl;
    
    // scan all cells to examine each pair of particles:
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int * region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
        
        // We consider each pair of objects (ii, jj) only once:
        
        FatPointListF & baseP = point_list(inx);
        FatLocusListF & baseL = locus_list(inx);
        
        setInteractions(meca, stiff, baseP, baseL);
        
        for ( int reg = 1; reg < nr; ++reg )
        {
            FatPointListF & sideP = point_list(inx+region[reg]);
            FatLocusListF & sideL = locus_list(inx+region[reg]);
            
            setInteractions(meca, stiff, baseP, baseL, sideP, sideL);
        }
    }
}


#else

/**
 Check interactions between the FatPoints contained in Pane `pan`.
 */
void PointGridF::setInteractions(Meca& meca, real stiff,
                                 const size_t pan) const
{
    // scan all cells to examine each pair of particles:
    for ( size_t inx = 0; inx < pGrid.nbCells(); ++inx )
    {
        int * region;
        int nr = pGrid.getRegion(region, inx);
        assert_true(region[0] == 0);
        
        // We consider each pair of objects (ii, jj) only once:
        
        FatPointListF & baseP = point_list(inx, pan);
        FatLocusListF & baseL = locus_list(inx, pan);
        
        setInteractions(meca, stiff, baseP, baseL);
        
        for ( int reg = 1; reg < nr; ++reg )
        {
            FatPointListF & sideP = point_list(inx+region[reg], pan);
            FatLocusListF & sideL = locus_list(inx+region[reg], pan);
            
            setInteractions(meca, stiff, baseP, baseL, sideP, sideL);
        }
    }
}


/**
 Check interactions between the FatPoints contained in Panes `pan1` and `pan2`,
 where ( pan1 != pan2 )
 */
void PointGridF::setInteractions(Meca& meca, real stiff,
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
        
        FatPointListF & baseP = point_list(inx, pan1);
        FatLocusListF & baseL = locus_list(inx, pan1);
        
        for ( int reg = 0; reg < nr; ++reg )
        {
            FatPointListF & sideP = point_list(inx+region[reg], pan2);
            FatLocusListF & sideL = locus_list(inx+region[reg], pan2);
            
            setInteractions(meca, stiff, baseP, baseL, sideP, sideL);
        }
        
        FatPointListF & baseP2 = point_list(inx, pan2);
        FatLocusListF & baseL2 = locus_list(inx, pan2);
        
        for ( int reg = 1; reg < nr; ++reg )
        {
            FatPointListF & sideP = point_list(inx+region[reg], pan1);
            FatLocusListF & sideL = locus_list(inx+region[reg], pan1);
            
            setInteractions(meca, stiff, baseP2, baseL2, sideP, sideL);
        }
    }
}


#endif

