// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.
#ifndef POINT_GRID_H
#define POINT_GRID_H

#include "grid.h"
#include "dim.h"
#include "vector.h"
#include "mecapoint.h"
#include "fiber_segment.h"
#include "array.h"

class Space;
class Modulo;
class Simul;
class Fiber;


/// number of panes in the steric engine
/** This should normally be set equal to 1, for optimal performance */
#define NUM_STERIC_PANES 1


/// Used for early exclusion of potential pairs, representing { position, interaction radius }
/** This uses single precision arithmetics, hopefully sufficient for exclusion tests */
class FatVector
{
public:

    float XX, YY, ZZ;
    float RR;
    
    FatVector() { XX = 0; YY = 0; ZZ = 0; RR = 0; }

    FatVector(Vector1 v, real r) { XX = v.XX; YY = 0; ZZ = 0; RR = r; }
    FatVector(Vector2 v, real r) { XX = v.XX; YY = v.YY; ZZ = 0; RR = r; }
    FatVector(Vector3 v, real r) { XX = v.XX; YY = v.YY; ZZ = v.ZZ; RR = r; }

    /// @return result of test `distance(this, arg) < sum_of_ranges`
    bool near(FatVector const& arg) const
    {
        float x = XX - arg.XX;
#if ( DIM == 1 )
        float r = RR + arg.RR;
        return ( std::abs(x) <= r );
#elif ( DIM == 2 )
        float y = YY - arg.YY;
        float r = RR + arg.RR;
        return ( x*x + y*y <= r*r );
#else
        float y = YY - arg.YY;
        float z = ZZ - arg.ZZ;
        float r = RR + arg.RR;
        return ( z*z + x*x <= r*r - y*y );
#endif
    }
};

//------------------------------------------------------------------------------

/// represents a Mecapoint for steric interactions
class FatPoint
{
    friend class PointGrid;
    
public:
    
    /// position of center
    FatVector pos_;
    
    /// indicates one vertex in a Mecable
    Mecapoint pnt_;

    /// equilibrium radius of the interaction (distance where force is zero)
    real      rad_;
    
    /// interaction range (maximum distance at which the force can operate)
    real      rge_;
    
public:
    
    FatPoint() {}
    
    FatPoint(Mecapoint const& p, real r, real e, Vector const& w)
    : pos_(w, e)
    {
        rad_ = r;
        rge_ = e;
        pnt_ = p;
    }
    
    /// position of center
    Vector cen() const { return pnt_.pos(); }
};

//------------------------------------------------------------------------------

/// represents the Segment of a Fiber for steric interactions
class FatLocus
{
    friend class PointGrid;
    
public:
    
    /// position of center
    FatVector pos_;
    
    /// indicates one segment of a Fiber
    FiberSegment seg_;

    /// equilibrium radius of the interaction (distance where force is zero)
    real rad_;
    
    /// interaction range (maximum distance at which the force can operate)
    real rge_;
    
public:
    
    FatLocus() {}
    
    FatLocus(FiberSegment const& p, real r, real e, real u, Vector const& w)
    : pos_(w, u)
    {
        rad_ = r;
        rge_ = e;
        seg_ = p;
    }
    
    /// true if the segment is the first of the Fiber
    bool isFirst() const { return seg_.isFirst(); }
    
    /// true if the segment is the last of the Fiber
    bool isLast() const { return seg_.isLast(); }
    
    /// position of point 1
    Vector pos1() const { return seg_.pos1(); }
    
    /// position of point 2
    Vector pos2() const { return seg_.pos2(); }
    
    /// offset = point2 - point1
    Vector diff() const { return seg_.diff(); }

    /// offset = point1 - point0
    Vector prevDiff() const { return seg_.fiber()->diffPoints(seg_.point()-1); }
    
    /// length of segment
    real len() const { return seg_.len(); }

    /// Mecapoint to point 1
    Mecapoint vertex1() const { return seg_.vertex1(); }
    
    /// Mecapoint to point 2
    Mecapoint vertex2() const { return seg_.vertex2(); }

};

//------------------------------------------------------------------------------

/// type for a list of FatPoint
typedef Array<FatPoint> FatPointList;

/// type for a list of FatLocus
typedef Array<FatLocus> FatLocusList;


/// a set of lists associated with the same location
class PointGridCell
{
    friend class PointGrid;
    
#if ( NUM_STERIC_PANES == 1 )
    
    /// unique steric pane
    FatPointList point_pane;
    
    /// unique steric pane
    FatLocusList locus_pane;
    
#else
    
    /// different steric panes
    FatPointList point_panes_0[NUM_STERIC_PANES];
    
    /// different steric panes
    FatLocusList locus_panes_0[NUM_STERIC_PANES];
    
    /// alias to the array of panes, with index 1 refering to point_panes_0[0]
    FatPointList * point_panes;
    
    /// alias to the array of panes, with index 1 refering to locus_panes_0[0]
    FatLocusList * locus_panes;
    
#endif
    
public:
    
#if ( NUM_STERIC_PANES == 1 )
    
    PointGridCell()
    {
    }
    
    /// clear all panes
    void clear()
    {
        point_pane.clear();
        locus_pane.clear();
    }
    
    size_t capacity() const
    {
        return point_pane.capacity() + locus_pane.capacity();
    }

#else
    
    PointGridCell() : point_panes(point_panes_0), locus_panes(locus_panes_0)
    {
        --point_panes;
        --locus_panes;
    }
    
    /// clear all panes
    void clear()
    {
        for ( size_t p = 1; p <= NUM_STERIC_PANES; ++p )
        {
            point_panes[p].clear();
            locus_panes[p].clear();
        }
    }
    
    FatPointList& point_list(size_t p)
    {
        assert_true( 0 < p && p <= NUM_STERIC_PANES );
        return point_panes[p];
    }
    
    FatLocusList& locus_list(size_t p)
    {
        assert_true( 0 < p && p <= NUM_STERIC_PANES );
        return locus_panes[p];
    }
    
    size_t capacity() const
    {
        size_t res = 0;
        for ( int i = 0; i < NUM_STERIC_PANES; ++i )
            res += point_panes[i].capacity() + locus_panes[i].capacity();
        return res;
    }

#endif
};

//------------------------------------------------------------------------------

/// PointGrid implements a Cell Lists approach to steric interactions
/**
 This implements a divide-and-conquer method to find particles that are within
 a certain cutoff distance from each other. In brief:
 - It covers the space with a Grid `pGrid`, initialized by `setGrid()`
 - A list of class `PointGridCell` is associated with each cell of this grid.
 - `PointGrid::add()` links a `FatPoint` or a `FatLocus` to the appropriate cell of the grid.
 - `PointGrid::setSterics()` checks all pairs of Point/Locus that may overlap,
    calculating their distance, and calling Meca::addLink() if they are interacting
 .
 The related class `LocusGrid`, is a simpler, streamline version of this class.
 For periodic boundary conditions, this follows the [Periodic wrapping] method.
 
 Check the [general introduction on Cell Lists](https://en.wikipedia.org/wiki/Cell_lists)
 */
class PointGrid
{
private:
    
    /// grid for divide-and-conquer strategies:
    Grid<PointGridCell, DIM> pGrid;
    
    /// Meca
    Meca& meca;
    
    /// stiffness
    real push;
    
    ///
    real pull;

private:
    
    /// check two Spheres
    void checkPP(FatPoint const&, FatPoint const&) const;
    
    /// check Sphere against Line segment
    void checkPL(FatPoint const&, FatLocus const&) const;
    
    /// check Line segment against Sphere
    void checkLL1(FatLocus const&, FatLocus const&) const;
    
    /// check Line segment against the terminal Sphere of a Fiber
    void checkLL2(FatLocus const&, FatLocus const&) const;
    
    /// check two Line segments
    void checkLL(FatLocus const&, FatLocus const&) const;
                 
    /// check all pairs between the two lists
    void setSterics0(FatPointList &, FatLocusList &) const;
    
    /// check all pairs between the two lists
    void setSterics0(FatPointList &, FatLocusList &,
                     FatPointList &, FatLocusList &) const;
    
    /// check all pairs between two lists, checking center-to-center distance
    void setStericsT(FatPointList &, FatLocusList &) const;
    
    /// check all pairs between two lists, checking center-to-center distance
    void setStericsT(FatPointList &, FatLocusList &,
                     FatPointList &, FatLocusList &) const;

#if ( NUM_STERIC_PANES == 1 )
    
    /// cell corresponding to position `w`
    FatPointList& point_list(Vector const& w) const
    {
        return pGrid.cell(w).point_pane;
    }
    
    /// cell corresponding to position `w`
    FatLocusList& locus_list(Vector const& w) const
    {
        return pGrid.cell(w).locus_pane;
    }
    
    /// cell corresponding to index `w`
    FatPointList& point_list(const size_t w) const
    {
        return pGrid.icell(w).point_pane;
    }
    
    /// cell corresponding to index `w`
    FatLocusList& locus_list(const size_t w) const
    {
        return pGrid.icell(w).locus_pane;
    }
    
    /// enter interactions into Meca
    void setSterics0() const;
    
    /// enter interactions into Meca
    void setStericsT() const;

#else
    
    /// cell corresponding to position `w`, and pane `p`
    FatPointList& point_list(Vector const& w, const size_t p) const
    {
        assert_true( 0 < p && p <= NUM_STERIC_PANES );
        return pGrid.cell(w).point_panes[p];
    }
    
    /// cell corresponding to position `w`, and pane `p`
    FatLocusList& locus_list(Vector const& w, const size_t p) const
    {
        assert_true( 0 < p && p <= NUM_STERIC_PANES );
        return pGrid.cell(w).locus_panes[p];
    }
    
    /// cell corresponding to index `c`, and pane `p`
    FatPointList& point_list(const size_t c, const size_t p) const
    {
        assert_true( 0 < p && p <= NUM_STERIC_PANES );
        return pGrid.icell(c).point_panes[p];
    }
    
    /// cell corresponding to index `c`, and pane `p`
    FatLocusList& locus_list(const size_t c, const size_t p) const
    {
        assert_true( 0 < p && p <= NUM_STERIC_PANES );
        return pGrid.icell(c).locus_panes[p];
    }
    
    /// enter interactions into Meca in one pane
    void setSterics0(size_t pan) const;
    
    /// enter interactions into Meca in one pane
    void setStericsT(size_t pan) const;
    
    /// enter interactions into Meca between two panes
    void setSterics0(size_t pan1, size_t pan2) const;
    
    /// enter interactions into Meca between two panes
    void setStericsT(size_t pan1, size_t pan2) const;

#endif
    
public:
    
    /// creator
    PointGrid(Meca& m) : meca(m) { push = 0; pull = 0; }
    
    /// number of panes
    static int nbPanes() { return NUM_STERIC_PANES; }

    /// set stiffness
    void stiffness(real h, real l) { push = h; pull = l; }

    /// define grid covering specified Space, given a minimal cell size requirement
    size_t setGrid(Space const*, real min_width);
    
    /// allocate memory for grid
    void createCells();
    
    /// true if the grid was initialized by calling setGrid()
    size_t hasGrid() const { return pGrid.hasCells(); }
    
    /// sum of allocated size of lists for all cells
    size_t capacity() const;

    /// clear the grid
    void clear() { pGrid.clear(); }
        
    /// place FiberSegment on the grid
    void add(Fiber const* fib, size_t inx, real rad, real rge, real sup) const
    {
        // link in cell containing the middle of the segment
        Vector w = fib->midPoint(inx, 0.5);
#if ( NUM_STERIC_PANES <= 1 )
        locus_list(w).emplace(FiberSegment(fib, inx), rad, rge, sup, w);
#else
        locus_list(w, fib->prop->steric).emplace(FiberSegment(fib, inx), rad, rge, sup, w);
#endif
    }
    
    /// place Mecapoint on the grid
    template <typename MECABLE>
    void add(MECABLE const* mec, size_t inx, real rad, real rge) const
    {
        Vector w = mec->posPoint(inx);
#if ( NUM_STERIC_PANES <= 1 )
        point_list(w).emplace(Mecapoint(mec, inx), rad, rge, w);
#else
        point_list(w, mec->prop->steric).emplace(Mecapoint(mec, inx), rad, rge, w);
#endif
    }

    /// enter interactions into Meca
    void setSterics() const;
    
#if ( NUM_STERIC_PANES > 1 )
    
    /// enter interactions into Meca in one pane
    void setSterics(size_t pan) const;
    
    /// enter interactions into Meca between two panes
    void setSterics(size_t pan1, size_t pan2) const;
    
#endif
    
    /// underlying spatial grid
    Map<DIM> const& map() const { return pGrid; }
    
    /// OpenGL display function
    void drawGrid() const;
};


#endif
