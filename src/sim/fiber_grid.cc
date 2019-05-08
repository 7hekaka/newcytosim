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
// this includes a naive implementation, which is slow but helpful for debugging
#   include "fiber_grid2.cc"
#else

/**
 Creates a grid where the dimensions of the cells are `max_step` at most.
 If the numbers of cells that need to be created is greater than `max_nb_cells`,
 the function returns 1 without building the grid.
 The return value is zero in case of success.
 
 The algorithm works with any value of `max_step` (the results are always correct),
 but `max_step` affects the efficiency (speed) of the algorithm:
 -if `max_step` is too small, paintGrid() will be slow,
 -if `max_step` is too large, tryToAttach() will be slow.
 A good compromise is to set `max_step` equivalent to the attachment distance,
 or at least to the size of the segments of the Fibers.
 */
unsigned FiberGrid::setGrid(Space const* space, real max_step)
{
    if ( max_step <= 0 )
        throw InvalidParameter("simul:binding_grid_step should be > 0");
    
    Vector inf, sup;
    space->boundaries(inf, sup);
    
    int n_cell[3] = { 1, 1, 1 };
    
    for ( int d = 0; d < DIM; ++d )
    {
        n_cell[d] = (int) ceil( ( sup[d] - inf[d] ) / max_step );
        
        if ( n_cell[d] < 0 )
            throw InvalidParameter("invalid space:boundaries");
        
        if ( modulo  &&  modulo->isPeriodic(d) )
        {
            //adjust the grid to match the edges exactly
            fGrid.setPeriodic(d, true);
        }
        else
        {
            //extend the grid by one cell on each side
            inf[d]    -= max_step;
            sup[d]    += max_step;
            n_cell[d] += 2;
        }
        
        if ( n_cell[d] <= 0 )
            n_cell[d] = 1;
    }

    //create the grid using the calculated dimensions:
    fGrid.setDimensions(inf, sup, n_cell);
    return fGrid.nbCells();
}


void FiberGrid::createCells()
{
    fGrid.createCells();
    //fGrid.printSummary(std::cerr, "FiberGrid");
}


size_t FiberGrid::hasGrid() const
{
    return fGrid.hasCells();
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
        //@todo write/call a specialized function for periodic: icellP1D
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
 
 'H' is calculated such that tryToAttach() finds any segment closer than 'gridRange':
 
 To determine H, we start from a relation on the sides of a triangle:
 (A) distance( GP, segment ) < distance( GP, X ) + distance( X, segment )
 where GP (grid-point) is the closest point on the grid to X.
 
 Since GP in tryToAttach() is the closest point on fGrid to X, we have:
 (B) distance( GP, X ) < 0.5 * fGrid.diagonalLength()
 
 Thus to find all rods for which:
 (B) distance( X, segment ) < gridRange
 we simply use:
 
     H = gridRange + 0.5 * fGrid.diagonalLength();
 
 Note: H is calculated by paintGrid() and gridRange by setFiberGrid().
 
 Linking all segments is done in an inverse way:
 for each segment, we cover all points of the grid inside a volume obtained
 by inflating the segment by the length H. We use for that the raterizer which
 calls the function paint() above.
 */

void FiberGrid::paintGrid(const Fiber * first, const Fiber * last)
{
    assert_true(hasGrid());
    assert_true(gridRange >= 0);
    
    fGrid.clear();
    const Vector offset(fGrid.inf());
    const Vector deltas(fGrid.delta());
    const real width = gridRange + 0.5 * fGrid.diagonalLength();
    
    //define the painting function used:
    void (*paint)(int, int, int, int, void*) = modulo ? paintCellPeriodic : paintCell;
    
    for ( const Fiber * fib = first; fib != last ; fib=fib->next() )
    {
        PaintJob job;
        job.grid = &fGrid;
        Vector P, Q = fib->posP(0);
        const real S = fib->segmentation();

        for ( unsigned n = 1; n < fib->nbPoints(); ++n )
        {
            P = Q;
            Q = fib->posP(n);
            job.segment.set(fib, n-1);

#if   ( DIM == 1 )
            Rasterizer::paintFatLine1D(paint, &job, P, Q, width, offset, deltas);
#elif ( DIM == 2 )
            Rasterizer::paintFatLine2D(paint, &job, P, Q, S, width, offset, deltas);
#else
            //Rasterizer::paintHexLine3D(paint, &job, P, Q, S, width, offset, deltas);
            Rasterizer::paintFatLine3D(paint, &job, P, Q, S, width, offset, deltas);
            //Rasterizer::paintBox3D(paint, &job, P, Q, width, offset, deltas);
#endif
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark - Access

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
    
    //get the cell index closest to the position in space:
    const auto indx = fGrid.index(place, 0.5);
    
    //get the list of rods associated with this cell:
    SegmentList & segments = fGrid.icell(indx);

    //randomize the list, to make attachments more fair:
    if ( segments.size() > 1 )
    {
        // randomize the list order
        //std::random_shuffle(segments.begin(), segments.end());
        segments.shuffle();
    }
    else if ( segments.empty() )
        return;
    
    //std::clog << "tryToAttach has " << segments.size() << " segments\n";
    
    for ( FiberSegment const& seg : segments )
    {
#if !TRICKY_HAND_ATTACHMENT
        if ( RNG.test(ha.prop->binding_rate_dt) )
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


/**
 This function is limited to the range given in paintGrid();
 */
FiberGrid::SegmentList FiberGrid::nearbySegments(Vector const& place, const real D, Fiber * exclude)
{
    if ( gridRange < 0 )
        throw InvalidParameter("the Grid was not initialized");
    if ( gridRange < D )
    {
        printf("gridRange = %.4f < range = %.4f\n", gridRange, D);
        throw InvalidParameter("the Grid maximum distance was exceeded");
    }
    
    SegmentList res;
    
    //get the grid node list index closest to the position in space:
    const auto indx = fGrid.index(place, 0.5);
    
    //get the list of rods associated with this cell:
    SegmentList & segments = fGrid.icell(indx);
    
    const real DD = D*D;
    for ( FiberSegment const& seg : segments )
    {
        if ( seg.fiber() == exclude )
            continue;
        
        real dis = INFINITY;
        seg.projectPoint(place, dis);
        
        if ( dis < DD )
            res.push_back(seg);
    }
    
    return res;
}


FiberSegment FiberGrid::closestSegment(Vector const& place)
{
    //get the cell index from the position in space:
    const auto indx = fGrid.index(place, 0.5);
    
    //get the list of rods associated with this cell:
    SegmentList & segments =  fGrid.icell(indx);
    
    FiberSegment res(nullptr, 0);
    real closest = 4 * gridRange * gridRange;
    
    for ( FiberSegment const& seg : segments )
    {
        //we compute the distance from the hand to the candidate rod,
        //and compare it to the best we have so far.
        real dis = INFINITY;
        seg.projectPoint(place, dis);
        
        if ( dis < closest )
        {
            closest = dis;
            res = seg;
        }
    }
    return res;
}


#endif


//============================================================================
//===                        TEST  ATTACHMENT                             ====
//============================================================================
#pragma mark - Test

#include <map>
#include "simul.h"

/**
Function testAttach() is given a position in space,
 it calls tryToAttach() from this position to check that:
 - attachement has equal probability to all targets,
 - no target is missed,
 - attachment are not made to targets that are beyond binding_range
 */
void FiberGrid::testAttach(FILE * out, const Vector pos, Fiber * start, HandProp const* hp)
{
    //create a test motor with a dummy HandMonitor:
    HandMonitor hm;
    Hand ha(hp, &hm);
    real dsq = hp->binding_range_sqr;
    
    typedef std::map < FiberSegment const*, int > map_type;
    map_type hits;
    
    //go through all the segments to find those close enough from pos:
    for ( Fiber * fib=start; fib; fib=fib->next() )
    {
        for ( unsigned p = 0; p < fib->nbSegments(); ++p )
        {
            FiberSegment seg(fib, p);
            real dis = INFINITY;
            seg.projectPoint(pos, dis);
            
            if ( dis < dsq )
                hits[&seg] = 0;
        }
    }
    
    const size_t targets = hits.size();
    
    if ( targets == 0 )
    {
        //fprintf(out, "no target here\n");
        return;
    }
    
    //call tryTyAttach NB times to check to which rods the Hand binds:
    const size_t NB = 100 * targets;
    for ( size_t n = 0; n < NB; ++n )
    {
        // we set 'prob' to bind immediately
        tryToAttach(pos, ha);
        if ( ha.attached() )
        {
            Interpolation inter = ha.fiber()->interpolate(ha.abscissa());
            FiberSegment seg(ha.fiber(), inter.point1());
            
            if ( hits.find(&seg) != hits.end() )
                ++hits[&seg];
            else
                hits[&seg] = -2;
            
            ha.detach();
        }
    }
    
    //detect segments that have been missed or mistargeted:
    int verbose = 1;
    for ( auto const& i : hits )
    {
        if ( i.second <= 50 )
            verbose = 1;
        if ( i.second < 0 )
            verbose = 2;
    }
    
    if ( verbose )
    {
        // print a summary of all targets:
        fprintf(out, "FiberGrid::testAttach %lu target(s) within %.3f um of", targets, hp->binding_range);
        pos.print(out);

#if ( 0 )
        const auto indx = fGrid.index(pos, 0.5);
        for ( FiberSegment const& seg : fGrid.icell(indx) )
            fprintf(out, "\n    target f%04d:%02i", seg->fiber()->identity(), seg->point());
#endif
        
        //report for all the segments that were targeted:
        for ( auto const& it : hits )
        {
            FiberSegment const* seg = it.first;
            Fiber const* fib = seg->fiber();
            real dis = INFINITY;
            real abs = seg->projectPoint(pos, dis);
            
            fprintf(out, "\n    rod f%04d:%02i at %5.3f um, abs %+.2f : ", fib->identity(), seg->point(), dis, abs);
            if ( hits[seg] == 0 )
                fprintf(out, "missed");
            else if ( hits[seg] < 0 )
                fprintf(out, "found, although out of range");
            else if ( hits[seg] > 0 )
                fprintf(out, "%-3i hits", hits[seg]);
        }
        fprintf(out, "\n");
    }
}

