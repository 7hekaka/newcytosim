// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020
#include "assert_macro.h"
#include "rasterizer.h"
#include "fiber_grid.h"
#include "exceptions.h"
#include "fiber_segment.h"
#include "fiber_site.h"
#include "messages.h"
#include "space.h"
#include "hand.h"
#include "hand_prop.h"
#include "cymdef.h"


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
index_t FiberGrid::setGrid(Vector inf, Vector sup, real max_step)
{
    assert_true(max_step > 0);
    modulo_ = modulo;
    
    index_t cnt[3] = { 1, 1, 1 };
    
    for ( int d = 0; d < DIM; ++d )
    {
        real n = ( sup[d] - inf[d] ) / max_step;
        
        if ( n < 0 )
            throw InvalidParameter("invalid space:boundaries");

        if ( modulo_ && modulo_->isPeriodic(d) )
        {
            //adjust the grid to match the edges exactly
            cnt[d] = std::max((index_t)1, (index_t)std::ceil(n));
            fGrid.setPeriodic(d, true);
        }
        else
        {
            //extend the grid by one cell on each side
            cnt[d]  = (index_t)std::ceil(n) + 2;
            inf[d] -= max_step;
            sup[d] += max_step;
        }
    }
    
    fGrid.setDimensions(inf, sup, cnt);
    return fGrid.nbCells();
}


void FiberGrid::createCells()
{
    if ( fGrid.createCells() )
        fGrid.printSummary(Cytosim::log, "   FiberGrid");
}


index_t FiberGrid::nbCells() const
{
    return fGrid.nbCells();
}


size_t FiberGrid::hasGrid() const
{
    return fGrid.hasCells();
}


index_t FiberGrid::nbTargets() const
{
    index_t res = 0;
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
    FiberGrid::SegmentGrid * grid;
    FiberSegment segment;
    PaintJob(FiberGrid::SegmentGrid* g, Fiber const* f) { grid=g; segment.set(f, 0); }
};


/**
 paintCell(x,y,z) adds a Segment to the SegmentList associated with
 the grid point (x,y,z).
 It is called by the rasterizer function paintThickLine().
 
 This version uses the fact that cells with consecutive
 X-coordinates should be consecutive also in the Grid
 */
void paintCell(const int x_inf, const int x_sup, const int y, const int z, void * arg)
{
    const auto* grid = static_cast<PaintJob*>(arg)->grid;
    const auto& seg = static_cast<PaintJob*>(arg)->segment;
    //printf("paint %i:%i (%i --> %i, %i, %i)\n", seg.identity(), seg.point(), x_inf, x_sup, y, z);

#if   ( DIM == 1 )
    FiberGrid::SegmentList * list = & grid->icell1D_clamped(x_inf);
    FiberGrid::SegmentList * last = & grid->icell1D_clamped(x_sup);
#elif ( DIM == 2 )
    FiberGrid::SegmentList * list = & grid->icell2D_clamped(x_inf, y);
    FiberGrid::SegmentList * last = & grid->icell2D_clamped(x_sup, y);
#else
    FiberGrid::SegmentList * list = & grid->icell3D_clamped(x_inf, y, z);
    FiberGrid::SegmentList * last = & grid->icell3D_clamped(x_sup, y, z);
#endif
    
    // Since all the lists are independent, they could be updated in parallel
    #pragma ivdep
    for ( ; list <= last; ++list )
        list->push_back(seg);
}


/** 
 paintCellPeriodic(x,y,z) adds a Segment in the SegmentList associated with
 the grid point (x,y,z). 
 It is called by the rasterizer function paintThickLine()
 */

void paintCellPeriodic(const int x_inf, const int x_sup, const int y, const int z, void * arg)
{
    auto* grid = static_cast<PaintJob*>(arg)->grid;
    const auto& seg = static_cast<PaintJob*>(arg)->segment;
    //printf("paintPeriodic %i:%i (%i --> %i, %i, %i)\n", seg.identity(), seg.point(), x_inf, x_sup, y, z);

    // safety, to avoid instances wehere x_inf was a very large negative number
    size_t w = grid->breadth(0);
    
    #pragma ivdep
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
#else
        grid->icell3D( x, y, z ).push_back(seg);
#endif
        if ( --w == 0 )
        {
            printf("aborted paintCellPeriodic %i:%i (%i --> %i, %i, %i)\n", seg.identity(), seg.point(), x_inf, x_sup, y, z);
            break;
        }
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
 (B) distance( GP, X ) < fGrid.cellRadius()
 
 Thus to find all rods for which:
 (B) distance( X, segment ) < grid::range
 we simply use:
 
     H = grid::range + fGrid.cellRadius();
 
 Note: H is calculated by paintGrid() and grid::range by estimateFiberGridStep().
 
 Linking all segments is done in an inverse way:
 for each segment, we cover all points of the grid inside a volume obtained
 by inflating the segment by the length H. We use for that the raterizer which
 calls the function paint() above.
 */

void FiberGrid::paintGrid(const Fiber * first, const Fiber * last, real range)
{
    assert_true(range >= 0);
    
    fGrid.clearCells();
    const Vector offset(fGrid.inf());
    const Vector deltas(fGrid.delta());
    const real width = range + fGrid.cellRadius();
    
    //define the painting function used:
    void (*paint)(int, int, int, int, void*) = modulo_ ? paintCellPeriodic : paintCell;
    
    for ( const Fiber * fib = first; fib != last ; fib=fib->next() )
    {
        // skip any Fiber on which binding is disabled:
        if ( 0 == fib->prop->binding_key )
            continue;
        
        PaintJob job(&fGrid, fib);
        Vector P, Q = fib->posP(0);
        const real iPQ = 1.0 / fib->segmentation();

        for ( unsigned n = 1; n < fib->nbPoints(); ++n )
        {
            P = Q;
            Q = fib->posP(n);
            job.segment.point(n-1);

#if ( DIM == 1 )
            Rasterizer::paintThickLine1D(paint, &job, P, Q, width, offset, deltas);
#elif ( DIM == 2 )
            Rasterizer::paintRectangle(paint, &job, P, Q, iPQ, width, offset, deltas);
#else
            //Rasterizer::paintHexagonalPrism(paint, &job, P, Q, iPQ, width, offset, deltas);
            Rasterizer::paintCuboid(paint, &job, P, Q, iPQ, width, offset, deltas);
            //Rasterizer::paintBox3D(paint, &job, P, Q, width, offset, deltas);
#endif
        }
    }
    
    //if ( first ) printf("FiberGrid has %lu segments in %lu cells\n", nbTargets(), nbCells());
}


//------------------------------------------------------------------------------
#pragma mark - Access

#if BIND_CLOSEST_FIBER

std::string FiberGrid::BindingTarget::to_string() const
{
    std::ostringstream oss;
    oss << " " << sit_;
    // convert distance squared to nm
    oss << " " << std::setprecision(1) << std::fixed << 1000 * sqrt(dis_);
    return oss.str();
}
    
/// qsort function comparing target segments from near to far
static int compareTargets(const void * A, const void * B)
{
    real a = static_cast<FiberGrid::BindingTarget const*>(A)->dis_;
    real b = static_cast<FiberGrid::BindingTarget const*>(B)->dis_;
    
    return ( a > b ) - ( b > a );
}


/**
 This will bind the given Hand to any Fiber found within `HandProp::binding_range`,
 with a probability `HandProp::binding_prob`
 
 NOTE:
 The distance at which Fibers are detected is limited to the range given in paintGrid()
 
 The operation can be done in parallel, if a thread_local target list is used (targets_),
 since the grid is read-only during the operation
 */
void FiberGrid::tryToAttach(Vector const& place, Hand& ha) const
{
    assert_true( hasGrid() );
    /*
     find the cell whose center is closest to the position, and
     get the list of segments associated with this cell in FiberGrid::paintGrid
     */
    SegmentList const& segments = fGrid.icell(fGrid.index(place, 0.5));
    // using here the member variable 'targets' as a temporary list:
    targets_.clear();
    
    // calculate distance to all targets
    const real sup = square(ha.property()->binding_range);
    for ( FiberSegment const& seg : segments )
    {
        if ( ha.keyMatch(seg.fiber()) )
        {
            real dis = INFINITY;
            real abs = seg.projectPoint(place, dis);
            if ( dis < sup )
            {
                // ATTENTION: convert `abs` relative to the segment to Fiber's abscissa
                FiberSite sit(seg.fiber(), seg.abscissa1()+abs);
                if ( ha.attachmentAllowed(sit) )
                {
                    /* Might be better to insert the new target at the right position,
                     to directly produce an ordered list, avoiding the sorting stage below */
                    targets_.emplace(sit, dis);
                    //std::clog << "   target " << sit << " at " << 1000*std::sqrt(dis) << " nm\n";
                }
                //else std::clog << "      xxx " << sit << " at " << 1000*std::sqrt(dis) << " nm\n";
            }
            //else std::clog << "      far " << seg << " at " << 1000*std::sqrt(dis) << " nm\n";
        }
    }
    
    // sort targets within range from close to distant:
    if ( targets_.size() > 1 )
        targets_.quick_sort(compareTargets);
    else if ( targets_.empty() )
        return;
    
    /**
     Instead of flipping a coin for each target, we could use a single random
     number to get the index of the next target that will bind, using a Poisson
     distribution */
    const uint64_t prob = 0x1p+32 * ha.property()->binding_prob;
    //std::clog << &ha << " target";
    for ( BindingTarget const& hit : targets_ )
    {
        //std::clog << " " << hit.to_string();
        if ( RNG.pint32() < prob )
        {
            ha.attach(hit.site());
            break;
        }
    }
    //std::clog << "\n";
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

    const uint64_t prob = 0x1p+32 * ha.property()->binding_prob;
    const real sup = square(ha.property()->binding_range);
    for ( FiberSegment const& seg : segments )
    {
        if ( RNG.pint32() < prob )
        {
            real dis = INFINITY;
            // Compute the distance from the hand to the rod, and abscissa of projection:
            real abs = seg.projectPoint(place, dis);      // always works
            //real abs = seg->projectPointF(place, dis);    // faster, but not compatible with periodic boundaries
            
            /*
             Compare to the maximum attachment range of the hand,
             and compare a newly tossed random number with 'prob'
             */
            if ( dis < sup )
            {
                FiberSite sit(seg.fiber(), seg.abscissa1()+abs);
                
                if ( ha.keyMatch(seg.fiber()) && ha.attachmentAllowed(sit) )
                {
                    ha.attach(sit);
                    return;
                }
            }
        }
    }
}

#endif


FiberGrid::SegmentList FiberGrid::nearbySegments(Vector const& place) const
{
    const auto inx = fGrid.index(place, 0.5);
    return fGrid[inx];
}

/**
 This function is limited to the range given in paintGrid();
 */
FiberGrid::SegmentList FiberGrid::nearbySegments(Vector const& place, const real DD, Fiber const* exclude) const
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
#include "hand_monitor.h"
#include "fiber_set.h"


/// used for debugging
index_t mingle(FiberSegment const& seg)
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
    const real sup = square(hp->binding_range);
    
    //check all the segments to find those close enough from pos:
    for ( Fiber const* fib=set.first(); fib; fib=fib->next() )
    {
        for ( index_t p = 0; p < fib->nbSegments(); ++p )
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
            Interpolation inter = ha.fiber()->interpolateAbs(ha.abscissa());
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
            ObjectID id = hit.first >> 16;   // undo mingle()
            index_t pt = hit.first & 65535;  // undo mingle()
            Fiber const* fib = set.identifyObject(id);
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


//------------------------------------------------------------------------------
#pragma mark - Display

#ifdef DISPLAY

#include "gym_view.h"
#include "gym_draw.h"
#include "gym_cap.h"

void drawBoundaries(Map<DIM> const&, float);

void FiberGrid::drawGrid() const
{
#if ( DIM <= 3 )
    gym::ref_view();
    gym::disableLighting();
    gym::color(0,0,1);
    drawBoundaries(fGrid, 0.5f);
#endif
}
#endif
