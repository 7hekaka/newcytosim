// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef FIBER_GRID_H
#define FIBER_GRID_H

#include "dim.h"
#include "vector.h"
#include "modulo.h"
#include "array.h"
#include "grid.h"
#include "fiber.h"
#include "fiber_segment.h"
#include "fiber_site.h"
//#include <vector>

class Simul;
class PropertyList;
class FiberSet;
class Modulo;
class Space;
class HandProp;
class Hand;


/// Divide-and-Conquer method to find all FiberSegment located near a given point in space
/**
A divide-and-conquer algorithm is used to find all segments of fibers close to a given point:
 
 -# It uses a grid 'fGrid' covering the space, initialized by setGrid().
    After initialization, each cell of the grid has an empty SegmentList (a list of FiberSegment*).
 -# clear() resets all lists on the grid
 -# paintGrid() distributes the segments specified in the arguments to the cell-associated SegmentList.
    One of the argument specifies a maximum distance to be queried (`max_range`).
    After the distribution, tryToAttach() is able to find any segment
    located at a distance `max_range` or less from any given point, in linear time.
 -# The function tryToAttach(X, ...) finds the cell on fGrid that contain `X`. 
    The associated SegmentList will then contains all the segments located at distance `max_range` or less from `X`. 
    tryToAttach() calls a function distanceSqr() sequentially for all the segments in this list,
    to calculate the exact Euclidian distance. 
    Finally, using a random number it tests the probability of attachment for the Hand given as argument.
 .
 
 @todo we could call paintGrid() only if the objects have moved by a certain threshold.
 This would work if we also extend the painted area around the rod, by the same threshold.
 we must also redo the paintGrid() when MT points are added or removed.
 
 Such algortihm should lead to large CPU gain, if calling clear() or paintGrid() is limiting,
 which is the case in particular in 3D, because the number of grid-cells is large.
*/

class FiberGrid 
{
public:
    
    /// variable-size list of FiberSegment
    typedef Array<FiberSegment> SegmentList;
    //typedef std::vector<FiberSegment> SegmentList;

    /// a grid of SegmentList
    typedef Grid<SegmentList, DIM> SegmentGrid;

#if BIND_CLOSEST_FIBER

    /// BindingTarget is used to calculate distance of segments to a point
    class BindingTarget
    {
    public:
        FiberSite sit_;   ///< The segment
        real      dis_;   ///< shortest distance squared from point to segment
        
        BindingTarget() {}
        BindingTarget(FiberSite const& s, real d) { sit_ = s; dis_ = d; }
        
        FiberSite site() const { return sit_; }
        
        friend bool operator <(const BindingTarget& a, const BindingTarget& b) { return a.dis_ < b.dis_; }

        std::string to_string() const;
    };

    /// list of SegmentHits used in `tryToAttach`
    mutable Array<BindingTarget> targets_;
    
#endif
    
private:
    
    /// grid for divide-and-conquer strategies:
    SegmentGrid fGrid;
    
    /// Object for periodic boundary conditions
    Modulo const* modulo_;
    
public:
    
    /// constructor
    FiberGrid() : modulo_(nullptr) { }

    /// set a grid to cover the specified Volume with cells of width `max_step` at most
    index_t setGrid(Vector inf, Vector sup, real max_step);
    
    /// true if the grid was initialized by setGrid(); return allocated size
    size_t hasGrid() const;

    /// allocate memory for the grid, with the dimensions set by setGrid()
    void createCells();
    
    /// number of cells in grid
    index_t nbCells() const;

    /// distribute the Fiber segments over the grid cells
    void paintGrid(const Fiber * first, const Fiber * last, real range);
    
    /// given a position, find nearby Fiber segments and test attachement of the provided Hand
    void tryToAttach(Vector const&, Hand&) const;

    /// return all Fiber segments located in cell grid at specified position
    SegmentList nearbySegments(Vector const& vec) const;

    /// return all Fiber segments located near P, within distance squared, except those belonging to `exclude`
    SegmentList nearbySegments(Vector const&, real disSqr, Fiber const* exclude = nullptr) const;
    
    /// Among the segments closer than grid:range, return the closest one
    FiberSegment closestSegment(Vector const&) const;
    
    /// total number of segments in grid
    index_t nbTargets() const;

    /// return segment list associated with cell containing 'pos'
    SegmentList& cellTargets(Vector const& pos) const
    {
        // get the cell index from the position in space:
        const index_t indx = fGrid.index(pos, 0.5);
        // get the list of rods associated with this cell:
        return fGrid.icell(indx);
    }
    
    /// test the results of tryToAttach(), at a particular position
    void testAttach(FILE *, Vector place, FiberSet const&, HandProp const*) const;
    
    /// underlying spatial grid
    Map<DIM> const& map() const { return fGrid; }

    /// OpenGL display function
    void drawGrid() const;
};


#endif
