// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "assert_macro.h"
#include "rasterizer.h"
#include "fiber_grid.h"
#include "exceptions.h"
#include "fiber_segment.h"
#include "fiber_site.h"
#include "messages.h"
#include "space.h"
#include "modulo.h"
#include "hand.h"
#include "hand_prop.h"
#include "simul.h"
#include "sim.h"

extern Modulo const* modulo;


#if ( 0 )
// use a naive implementation, which is slow but helpful for debugging
#   include "fiber_grid2.cc"
#else

/**
 Creates a grid where the dimensions of the cells are `max_step` at most.
 @returns the numbers of cells to be created with this cell size.
 
 The results will be correct for any value of `max_step`, but efficiency of the
 algorithm is affected by `max_step`:
     -if `max_step` is too small, paintGrid() will be slow,
     -if `max_step` is too large, tryToAttach() will be slow.
 A compromise is to adjust `max_step` to the length of the segments.
 */
size_t FiberGrid::setGrid(Space const* space, real max_step)
{
    assert_true(max_step > 0);
    Vector inf, sup;
    space->boundaries(inf, sup);
    
    size_t cnt[3] = { 1, 1, 1 };
    
    for ( size_t d = 0; d < DIM; ++d )
    {
        real n = ( sup[d] - inf[d] ) / max_step;
        
        if ( n < 0 )
            throw InvalidParameter("invalid space:boundaries");

        if ( modulo  &&  modulo->isPeriodic(d) )
        {
            //adjust the grid to match the edges exactly
            cnt[d] = std::max((size_t)1, (size_t)ceil(n));
            fGrid.setPeriodic(d, true);
        }
        else
        {
            //extend the grid by one cell on each side
            cnt[d]  = (size_t)ceil(n) + 2;
            inf[d] -= max_step;
            sup[d] += max_step;
        }
    }
    
    fGrid.setDimensions(inf, sup, cnt);
    return fGrid.nbCells();
}


void FiberGrid::createCells()
{
    fGrid.createCells();
    fGrid.printSummary(Cytosim::log, "BindingGrid");
}


size_t FiberGrid::nbCells() const
{
    return fGrid.nbCells();
}


size_t FiberGrid::hasGrid() const
{
    return fGrid.hasCells();
}


size_t FiberGrid::nbTargets() const
{
    size_t res = 0;
    SegmentList * ptr = fGrid.data();
    const SegmentList * end = ptr + fGrid.nbCells();
    
    for ( ; ptr < end; ++ptr )
        res += ptr->size();
    
    return res;
}

//------------------------------------------------------------------------------
#pragma mark - Paint

/// Structure used by FiberGrid::paintGrid to find Hand's attachement
struct PaintJob
{
    FiberGrid::grid_type * grid;
    FiberSegment segment;
};


/**
 paintCell(x,y,z) adds a Segment to the SegmentList associated with
 the grid point (x,y,z).
 It is called by the rasterizer function paintFatLine().
 
 This version uses the fact that cells with consecutive
 X-coordinates should be consecutive also in the Grid
 */
void paintCell(const int x_inf, const int x_sup, const int y, const int z, void * arg)
{
    auto* grid = static_cast<PaintJob*>(arg)->grid;
    const auto& seg = static_cast<PaintJob*>(arg)->segment;
    //printf("paint %p in (%i to %i, %i, %i)\n", seg, x_inf, x_sup, y, z);

#if   ( DIM == 1 )
    FiberGrid::SegmentList * inf = & grid->icell1D( x_inf );
    FiberGrid::SegmentList * sup = & grid->icell1D( x_sup );
#elif ( DIM == 2 )
    FiberGrid::SegmentList * inf = & grid->icell2D( x_inf, y );
    FiberGrid::SegmentList * sup = & grid->icell2D( x_sup, y );
#else
    FiberGrid::SegmentList * inf = & grid->icell3D( x_inf, y, z );
    FiberGrid::SegmentList * sup = & grid->icell3D( x_sup, y, z );
#endif
    
    // Since all the lists are independent, they can be updated in parallel
    # pragma ivdep
    for ( FiberGrid::SegmentList * list = inf; list <= sup; ++list )
        list->push_back(seg);
}


/** 
 paintCellPeriodic(x,y,z) adds a Segment in the SegmentList associated with
 the grid point (x,y,z). 
 It is called by the rasterizer function paintFatLine()
 */

void paintCellPeriodic(const int x_inf, const int x_sup, const int y, const int z, void * arg)
{
    auto* grid = static_cast<PaintJob*>(arg)->grid;
    const auto& seg = static_cast<PaintJob*>(arg)->segment;
    //printf("paint %p in (%i to %i, %i, %i)\n", seg, x_inf, x_sup, y, z);
    
    # pragma ivdep
    for ( int x = x_inf; x <= x_sup; ++x )
    {
        /*
        Here icell3() will call pack3D() to calculate the index, switching
        3 times to always select the periodic image() function, which is wastefull.
        @todo write/call a specialized function for periodic: icellP1D
         */
#if   ( DIM == 1 )
        grid->icell1D( x ).push_back(seg);
#elif ( DIM == 2 )
        grid->icell2D( x, y ).push_back(seg);
#elif ( DIM == 3 )
        grid->icell3D( x, y, z ).push_back(seg);
#endif
    }
}


/**
paintGrid(first_fiber, last_fiber) links all segments found in 'fiber' and its
 descendant, in the point-list GP that match distance(GP, segment) < H.
 
 'H' is calculated such that tryToAttach() finds any segment closer than 'grid:range':
 
 To determine H, we start from a relation on the sides of a triangle:
 (A) distance( GP, segment ) < distance( GP, X ) + distance( X, segment )
 where GP (grid-point) is the closest point on the grid to X.
 
 Since GP in tryToAttach() is the closest point on fGrid to X, we have:
 (B) distance( GP, X ) < 0.5 * fGrid.diagonalLength()
 
 Thus to find all rods for which:
 (B) distance( X, segment ) < grid::range
 we simply use:
 
     H = grid::range + 0.5 * fGrid.diagonalLength();
 
 Note: H is calculated by paintGrid() and grid::range by setFiberGrid().
 
 Linking all segments is done in an inverse way:
 for each segment, we cover all points of the grid inside a volume obtained
 by inflating the segment by the length H. We use for that the raterizer which
 calls the function paint() above.
 */

void FiberGrid::paintGrid(const Fiber * first, const Fiber * last, real range)
{
    assert_true(hasGrid());
    assert_true(range >= 0);
    
    fGrid.clear();
    const Vector offset(fGrid.inf());
    const Vector deltas(fGrid.delta());
    const real width = range + 0.5 * fGrid.diagonalLength();
    
    //define the painting function used:
    void (*paint)(int, int, int, int, void*) = modulo ? paintCellPeriodic : paintCell;
    
    for ( const Fiber * fib = first; fib != last ; fib=fib->next() )
    {
        // do not paint fibers on which binding is disabled:
        if ( 0 == fib->prop->binding_key )
            continue;
        PaintJob job;
        job.grid = &fGrid;
        Vector P, Q = fib->posP(0);
        const real iPQ = 1.0 / fib->segmentation();

        for ( size_t n = 1; n < fib->nbPoints(); ++n )
        {
            P = Q;
            Q = fib->posP(n);
            job.segment.set(fib, n-1);

#if   ( DIM == 1 )
            Rasterizer::paintFatLine1D(paint, &job, P, Q, width, offset, deltas);
#elif ( DIM == 2 )
            Rasterizer::paintFatLine2D(paint, &job, P, Q, iPQ, width, offset, deltas);
#else
            //Rasterizer::paintHexLine3D(paint, &job, P, Q, iPQ, width, offset, deltas);
            Rasterizer::paintFatLine3D(paint, &job, P, Q, iPQ, width, offset, deltas);
            //Rasterizer::paintBox3D(paint, &job, P, Q, width, offset, deltas);
#endif
        }
    }
    
    //Cytosim::log("distributed %lu segments over %lu cells\n", nbTargets(), nbCells());
}


//------------------------------------------------------------------------------
#pragma mark - Access

#if ATTACH_CLOSEST_FIBER

class HeavySegment
{
public:
    FiberSegment seg_;
    real         dis_;
    real         abs_;
    HeavySegment() { }
    HeavySegment(FiberSegment const& s, real d, real a) { seg_ = s; dis_ = d; abs_ = a; }
};


typedef Array<HeavySegment> HeavySegmentList;

/// used to qsort segments according to distance
int compareSegments(const void * a, const void * b)
{
    real ad = ((HeavySegment const*)(a))->dis_;
    real bd = ((HeavySegment const*)(b))->dis_;
    
    return ( ad > bd ) - ( bd > ad );
}


/**
 This will bind the given Hand to any Fiber found within `binding_range`, with a
 probability that is encoded in `prob`.
 The test is `RNG.pint() < prob`, and with 'prob = 1<<30', the chance is 1/4.
 The result is thus stochastic, and will depend on the number of Fiber
 within the range, but it will saturate if there are more than '4' possible targets.
 
 NOTE:
 The distance at which Fibers are detected is limited to the range given in paintGrid()
 */
void FiberGrid::tryToAttach(Vector const& place, Hand& ha) const
{
    assert_true( hasGrid() );
    /*
     find the cell whose center is closest to the position, and
     get the list of segments associated with this cell in FiberGrid::paintGrid
     */
    SegmentList & segments = fGrid.icell(fGrid.index(place, 0.5));
    HeavySegmentList targets;

    // calculate distance to all targets
    for ( FiberSegment const& seg : segments )
    {
        real dis_ = INFINITY;
        real abs_ = seg.projectPoint(place, dis_);
        if ( dis_ < ha.prop->binding_range_sqr )
            targets.emplace_back(seg, dis_, abs_);
    }
    
    // sort targets within range from close to distant:
    if ( targets.size() > 1 )
        targets.sort(compareSegments);
    else if ( targets.empty() )
        return;
    
    //std::clog << "tryToAttachClosest has " << targets.size() << " targets\n";

    /**
     Instead of flipping a coin for each target, we could use a single random
     number to get the index of the next target that will bind, using a Poisson
     distribution */
    const real prob = ha.prop->binding_prob;
    for ( HeavySegment const& target : targets )
    {
        //printf("    trying segment F%u:%lu dis %.6f\n", seg.fiber()->identity(), seg.point(), seg.dis_);
        if ( RNG.test(prob) )
        {
            FiberSegment const& seg = target.seg_;
            Fiber * fib = const_cast<Fiber*>(seg.fiber());
            FiberSite pos(fib, seg.abscissa1()+target.abs_);
            if ( ha.attachmentAllowed(pos) )
            {
                ha.attach(pos);
                return;
            }
        }
    }
}

#else


void FiberGrid::tryToAttach(Vector const& place, Hand& ha) const
{
    assert_true( hasGrid() );
    /*
     find the cell whose center is closest to the position, and
     get the list of segments associated with this cell in FiberGrid::paintGrid
     */
    SegmentList & segments = fGrid.icell(fGrid.index(place, 0.5));
    
    // randomize the list, to make attachments more fair:
    if ( segments.size() > 1 )
        segments.shuffle();
    else if ( segments.empty() )
        return;
    
    //std::clog << "tryToAttach has " << segments.size() << " targets\n";

    for ( FiberSegment const& seg : segments )
    {
#if !TRICKY_HAND_ATTACHMENT
        if ( RNG.test(ha.prop->binding_prob) )
#else
        if ( RNG.flip_8th() )
#endif
        {
            real dis = INFINITY;
            // Compute the distance from the hand to the rod, and abscissa of projection:
            real abs = seg.projectPoint(place, dis);      // always works
            //real abs = seg->projectPointF(place, dis);    // faster, but not compatible with periodic boundaries
            
            /*
             Compare to the maximum attachment range of the hand,
             and compare a newly tossed random number with 'prob'
             */
            if ( dis < ha.prop->binding_range_sqr )
            {
                Fiber * fib = const_cast<Fiber*>(seg.fiber());
                FiberSite pos(fib, seg.abscissa1()+abs);
                
                if ( ha.attachmentAllowed(pos) )
                {
                    ha.attach(pos);
                    return;
                }
            }
        }
    }
}

#endif

/**
 This function is limited to the range given in paintGrid();
 */
FiberGrid::SegmentList FiberGrid::nearbySegments(Vector const& place, const real DD, Fiber * exclude) const
{
    SegmentList res;
    
    //get the grid node list index closest to the position in space:
    const auto indx = fGrid.index(place, 0.5);
    
    //get the list of rods associated with this cell:
    for ( FiberSegment const& seg : fGrid.icell(indx) )
    {
        if ( seg.fiber() != exclude )
        {
            real dis = INFINITY;
            seg.projectPoint(place, dis);
            
            if ( dis < DD )
                res.push_back(seg);
        }
    }
    
    return res;
}


FiberSegment FiberGrid::closestSegment(Vector const& place) const
{
    //get the cell index from the position in space:
    const auto indx = fGrid.index(place, 0.5);
    
    FiberSegment res(nullptr, 0);
    real hit = INFINITY;
    
    //get the list of rods associated with this cell:
    for ( FiberSegment const& seg : fGrid.icell(indx) )
    {
        //we compute the distance from the hand to the candidate rod,
        //and compare it to the best we have so far.
        real dis = INFINITY;
        seg.projectPoint(place, dis);
        
        if ( dis < hit )
        {
            hit = dis;
            res = seg;
        }
    }
    return res;
}


#endif


//==============================================================================
//===                        TEST  ATTACHMENT                               ====
//==============================================================================
#pragma mark - Test

#include <map>
#include "simul.h"


/// used for debugging
unsigned mingle(FiberSegment const& seg)
{
    return ( seg.fiber()->identity() << 16 ) | seg.point();
}

/**
Function testAttach() is given a position in space,
 it calls tryToAttach() from this position to check that:
 - attachement has equal probability to all targets,
 - no target is missed,
 - attachment are not made to targets that are beyond binding_range
 */
void FiberGrid::testAttach(FILE* out, const Vector pos, FiberSet const& set, HandProp const* hp) const
{
    typedef std::map < unsigned, int > map_type;
    map_type hits;

    // create a test Hand with a dummy HandMonitor:
    HandMonitor hm;
    Hand ha(hp, &hm);
    real sup = square(hp->binding_range);
    
    //check all the segments to find those close enough from pos:
    for ( Fiber const* fib=set.first(); fib; fib=fib->next() )
    {
        for ( size_t p = 0; p < fib->nbSegments(); ++p )
        {
            FiberSegment seg(fib, p);
            real dis = INFINITY;
            seg.projectPoint(pos, dis);
            if ( dis < sup )
                hits[mingle(seg)] = 0;
        }
    }
    
    const size_t n_targets = hits.size();
    // call tryTyAttach 100 times per target:
    for ( size_t n = 0; n < n_targets; ++n )
    for ( size_t i = 0; i < 100; ++i )
    {
        tryToAttach(pos, ha);
        if ( ha.attached() )
        {
            Interpolation inter = ha.fiber()->interpolate(ha.abscissa());
            FiberSegment seg(ha.fiber(), inter.point1());
            
            if ( hits.find(mingle(seg)) != hits.end() )
                ++hits[mingle(seg)];
            else
                hits[mingle(seg)] = -2;
            //fprintf(out, "   attached to f%04d abscissa %7.3f\n", ha.fiber()->identity(), ha.abscissa());
            ha.detach();
        }
    }
    
    if ( hits.empty() )
        return;

    //detect segments that have been missed or mistargeted:
    bool doprint = false;
    for ( auto const& i : hits )
        doprint |= ( i.second < 50 );
    
    if ( doprint )
    {
        // print a summary of all targets:
        fprintf(out, "FiberGrid::testAttach %lu target(s) within %.3f um of", n_targets, hp->binding_range);
        pos.println(out);
#if ( 0 )
        //report content of grid's list
        const auto indx = fGrid.index(pos, 0.5);
        for ( FiberSegment const& seg : fGrid.icell(indx) )
            fprintf(out, "    target f%04d:%02i\n", seg.fiber()->identity(), seg.point());
#endif
        //report for all the segments that were targeted:
        for ( auto const& hit : hits )
        {
            ObjectID id = hit.first >> 16;    // opposite of mingle()
            unsigned pt = hit.first & 65535;  // opposite of mingle()
            Fiber const* fib = set.findID(id);
            FiberSegment seg(fib, pt);
            real dis = INFINITY;
            real abs = seg.projectPoint(pos, dis);
            
            fprintf(out, "    rod f%04d:%02i at %12.7f um, abs %+.2f : ", id, pt, dis, abs);
            if ( hit.second == 0 )
                fprintf(out, "missed");
            else if ( hit.second < 0 )
                fprintf(out, "found, although out of range");
            else if ( hit.second > 0 )
                fprintf(out, "%-3i hits", hit.second);
            fprintf(out, "\n");
        }
    }
}

//==============================================================================
#pragma mark - Display

#ifdef DISPLAY

#  include "grid_display.h"

void FiberGrid::draw() const
{
    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glColor4f(0, 1, 1, 1);
    glLineWidth(0.5);
    drawEdges(fGrid);
    glPopAttrib();
}
#endif
