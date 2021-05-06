// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef LOCUS_GRID_H
#define LOCUS_GRID_H

#include "grid.h"
#include "dim.h"
#include "vector.h"
#include "mecapoint.h"
#include "fiber.h"
#include "fiber_segment.h"
#include "array.h"

class Space;
class Modulo;
class Simul;
class Mecable;
class FiberSegment;


/// number of panes in the steric engine
/** This should normally be set equal to 1, for optimal performance */
#define MAX_STERIC_PANES 1


/// Used for early exclusing of potential pairs, representing { position, interaction radius }
/** This uses single precision arithmetics, hopefully sufficient for exclusion tests */
class BigVector
{
public:

    float XX, YY, ZZ;
    float RR;
    
    BigVector() { XX = 0; YY = 0; ZZ = 0; RR = 0; }

    BigVector(Vector1 v, real r) { XX = v.XX; YY = 0; ZZ = 0; RR = r; }
    BigVector(Vector2 v, real r) { XX = v.XX; YY = v.YY; ZZ = 0; RR = r; }
    BigVector(Vector3 v, real r) { XX = v.XX; YY = v.YY; ZZ = v.ZZ; RR = r; }

    float const* data() const { return &XX; }
    
    /// @return result of test `distance(this, arg) < sum_of_ranges`
    bool near(BigVector const& arg) const
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
        return ( x*x + y*y < r*r - z*z );
#endif
    }
};


/// represents the Segment of a Fiber for steric interactions
class BigLocus
{
    friend class LocusGrid;
    
public:
    
    /// position of center, and radius of interaction
    BigVector pos_;

    /// Fiber to which the segment belongs to
    Mecable const* obj_;
    
    /// equilibrium radius of the interaction (distance where force is zero)
    real     rad_;

    /// index of segment's first point
    unsigned vix_;
    
public:
    
    BigLocus() {}
    
    BigLocus(Mecable const* f, size_t i, real r, real e, Vector const& w)
    : pos_(w, e)
    {
        obj_ = f;
        rad_ = r;
        vix_ = static_cast<unsigned>(i);
        assert_true( i == vix_ );
    }
    
    /// position of center
    Vector cen() const { return obj_->posPoint(vix_); }
    //Vector cen() const { return Vector(&pos_.XX); }

    /// position of point 1
    Vector pos1() const { return obj_->posPoint(vix_); }
    
    /// position of point 2
    Vector pos2() const { return obj_->posPoint(vix_+1); }
    
    /// offset = point2 - point1
    Vector diff() const { return obj_->diffPoints(vix_); }
    
    /// offset = point1 - point0
    Vector prevDiff() const { return obj_->diffPoints(vix_-1); }
    
    /// length of segment
    real len() const { return static_cast<Fiber const*>(obj_)->segmentation(); }

    /// true if the segment is the first of the Fiber
    bool isFirst() const { return ( vix_ == 0 ); }
    
    /// true if the segment is the last of the Fiber
    bool isLast() const { return ( vix_+2 == obj_->nbPoints() ); }

    /// Mecapoint to point 1
    Mecapoint vertex1() const { return Mecapoint(obj_, vix_); }
    
    /// Mecapoint to point 2
    Mecapoint vertex2() const { return Mecapoint(obj_, vix_+1); }
    
    /// FiberSegment
    FiberSegment segment() const { return FiberSegment(static_cast<Fiber const*>(obj_), vix_); }

};

/// we used this alias for clarity and backward compatibility
typedef BigLocus BigPoint;

//------------------------------------------------------------------------------

/// a list containing BigLocus and BigPoint
/**
 The Fiber segments are contained in the first part of the list, index in [0, border[
 The other elements are in the second part, index in [border, end[
 */
class BigLocusList
{
    friend class LocusGrid;

    /// the list containing objects
    Array<BigLocus> pane;
    
    /// index of first non Fiber element in list
    size_t border;
    
public:
    
    /// constructor
    BigLocusList()
    {
        border = 0;
    }
    
    /// clear all panes
    void clear()
    {
        pane.clear();
        border = 0;
    }
    
    /// mark the edge where non-Fiber elements start
    void mark() { border = pane.size(); }
    
    /// number of elements in list
    size_t size() const { return pane.size(); }
    
    /// number of BigLocus in list
    size_t num_locus() const { return border; }
    
    /// number of BigPoints in list
    size_t num_points() const { return pane.size() - border; }

    /// first element in list
    BigLocus const* begin() const { return pane.begin(); }
    
    /// first BigPoint in list
    BigLocus const* pre_middle() const { return pane.begin() + (border & (~3UL)); }

    /// first BigPoint in list
    BigLocus const* middle() const { return pane.data() + border; }
    
    /// one past last element in list
    BigLocus const* pre_end() const { return pane.begin() + (pane.size() & (~3UL)); }
    
    /// one past last element in list
    BigLocus const* end() const { return pane.end(); }
    
    /// reference to Object at index ii (val_[ii])
    BigLocus const& operator[](const size_t i) const { return pane.at(i); }

    /// allocated size
    size_t capacity() const { return pane.capacity(); }

};


#if ( MAX_STERIC_PANES > 1 )

/// a set of lists associated with one cell of the grid
class LocusGridCell
{
    friend class LocusGrid;
    
    /// different steric panes
    BigLocusList panes_0[MAX_STERIC_PANES];
    
    /// alias to the array of panes, to use indices starting from 1
    BigLocusList * panes;
    
public:
    
    LocusGridCell() : panes(panes_0)
    {
        --panes;
    }
    
    /// clear all panes
    void clear()
    {
        for ( size_t p = 1; p <= MAX_STERIC_PANES; ++p )
            panes[p].clear();
    }
    
    BigLocusList& cell_list(size_t p)
    {
        assert_true( 0 < p && p <= MAX_STERIC_PANES );
        return panes[p];
    }
    
    void mark() const
    {
        for ( size_t p = 1; p <= MAX_STERIC_PANES; ++p )
            panes[p].mark();
    }

    size_t capacity() const
    {
        size_t res = 0;
        for ( int i = 0; i < MAX_STERIC_PANES; ++i )
            res += panes[i].capacity();
        return res;
    }
};

#endif

//------------------------------------------------------------------------------

/// LocusGrid implements a *Cell Lists* approach to steric interactions
/**
 This implements a divide-and-conquer method to find particles that are within a
 certain cutoff distance from each other. In brief:
 - It covers the space with a Grid `pGrid`, initialized by `setGrid()`
 - A list of class `LocusGridCell` is associated with each cell of `pGrid`.
 - `LocusGrid::add()` links `BigLocus` or `BigLocus` to the appropriate cell of the grid.
 - `LocusGrid::setInteractions()` checks all pairs of particles that may overlap,
    calculating their actual distance, and calling Meca::addLink() as necessary
 .
 Compared to PointGrid, LocusGrid only supports repulsive interactions.
 For periodic boundary conditions, this follows the [Periodic wrapping] method.
 
 Check the [general introduction on Cell Lists](https://en.wikipedia.org/wiki/Cell_lists)
 */
class LocusGrid
{
private:
    
#if ( MAX_STERIC_PANES == 1 )
    /// grid for divide-and-conquer strategies:
    Grid<BigLocusList, DIM> pGrid;
#else
    /// grid for divide-and-conquer strategies:
    Grid<LocusGridCell, DIM> pGrid;
#endif

private:
    
    /// check two Spheres
    static void checkPP(Meca&, real stiff, BigPoint const&, BigPoint const&);
    
    /// check Sphere against Line segment
    static void checkPL(Meca&, real stiff, BigPoint const&, BigLocus const&);
    
    /// check Line segment against Sphere
    static void checkLL1(Meca&, real stiff, BigLocus const&, BigLocus const&);
    
    /// check Line segment against the terminal Sphere of a Fiber
    static void checkLL2(Meca&, real stiff, BigLocus const&, BigLocus const&);
    
    /// check two Line segments
    static void checkLL(Meca&, real stiff, BigLocus const&, BigLocus const&);
    
    
    /// check all pairs between the two lists
    static void setSterics0(Meca&, real stiff, BigLocusList const&);
    
    /// check all pairs between the two lists
    static void setSterics0(Meca&, real stiff, BigLocusList const&, BigLocusList const&);
    
    /// check all pairs between the two lists, checking center-to-center distance
    static void setStericsT(Meca&, real stiff, BigLocusList const&);
    
    /// check all pairs between the two lists, checking center-to-center distance
    static void setStericsT(Meca&, real stiff, BigLocusList const&, BigLocusList const&);
    
    /// check all pairs between the two lists, checking center-to-center distance
    static void setStericsU(Meca&, real stiff, BigLocusList const&, BigLocusList const&);
    
    /// check all pairs between the two lists, checking center-to-center distance
    static void setStericsX(Meca&, real stiff, BigLocusList const&, BigLocusList const&);

#if ( MAX_STERIC_PANES == 1 )
    
    /// cell corresponding to position `w`, and pane `p`
    BigLocusList& cell_list(Vector const& w) const
    {
        return pGrid.cell(w);
    }
    
    /// cell corresponding to index `w`, and pane `p`
    BigLocusList& cell_list(const size_t w) const
    {
        return pGrid.icell(w);
    }
    
    /// enter interactions into Meca with given stiffness
    void setSterics0(Meca&, real stiff) const;
    
    /// enter interactions into Meca with given stiffness
    void setStericsT(Meca&, real stiff) const;

#else
    
    /// cell corresponding to position `w`, and pane `p`
    BigLocusList& cell_list(Vector const& w, const size_t p) const
    {
        assert_true( 0 < p && p <= MAX_STERIC_PANES );
        return pGrid.cell(w).panes[p];
    }
    
    /// cell corresponding to index `c`, and pane `p`
    BigLocusList& cell_list(const size_t c, const size_t p) const
    {
        assert_true( 0 < p && p <= MAX_STERIC_PANES );
        return pGrid.icell(c).panes[p];
    }
    
    /// enter interactions into Meca in one panes with given parameters
    void setSterics0(Meca&, real stiff, size_t pan) const;
    
    /// enter interactions into Meca in one panes with given parameters
    void setStericsT(Meca&, real stiff, size_t pan) const;
    
    /// enter interactions into Meca between two panes with given parameters
    void setSterics0(Meca&, real stiff, size_t pan1, size_t pan2) const;
    
    /// enter interactions into Meca between two panes with given parameters
    void setStericsT(Meca&, real stiff, size_t pan1, size_t pan2) const;

#endif
    
public:
    
    /// creator
    LocusGrid();
    
    /// define grid covering specified Space, given a minimal cell size requirement
    size_t setGrid(Space const*, real min_width);
    
    /// allocate memory for grid
    void createCells();
    
    /// true if the grid was initialized by calling setGrid()
    size_t hasGrid() const { return pGrid.hasCells(); }
    
    void mark() const;

    /// sum of allocated size of lists for all cells
    size_t capacity() const;

    /// clear the grid
    void clear() { pGrid.clear(); }
    
#if ( MAX_STERIC_PANES == 1 )
    
    /// place Mecable vertex on the grid
    void add(Mecable const* m, size_t i, real rad)
    {
        Vector w = m->posPoint(i);
        cell_list(w).pane.emplace(m, i, rad, rad, w);
    }
    
    // link in the cell containing the middle of the segment:
    void add(Fiber const* f, size_t i, real rad, real rge)
    {
        Vector w = f->posPoint(i, 0.5);
        cell_list(w).pane.emplace(f, i, rad, rge, w);
    }
    
    /// enter interactions into Meca with given stiffness
    void setInteractions(Meca&, real stiff) const;
    
#else

    /// place Mecable vertex on the grid
    void add(size_t pan, Mecable const* mec, size_t inx, real rad)
    {
        if ( pan == 0 || pan > MAX_STERIC_PANES )
            throw InvalidParameter("point:steric is out-of-range");
        Vector w = mec->posPoint(inx);
        cell_list(w, pan).pane.emplace(mec, inx, rad, rad, w);
    }
    
    // link in the cell containing the middle of the segment:
    void add(size_t pan, Fiber const* fib, size_t inx, real rad, real sup)
    {
        if ( pan == 0 || pan > MAX_STERIC_PANES )
            throw InvalidParameter("line:steric is out-of-range");
        Vector w = fib->posPoint(inx, 0.5);
        cell_list(w, pan).pane.emplace(fib, inx, rad, sup, w);
    }
    
    /// enter interactions into Meca in one panes with given parameters
    void setInteractions(Meca&, real stiff, size_t pan) const;
    
    /// enter interactions into Meca between two panes with given parameters
    void setInteractions(Meca&, real stiff, size_t pan1, size_t pan2) const;
    
#endif

    /// underlying spatial grid
    Map<DIM> const& map() const { return pGrid; }
    
    /// OpenGL display function
    void drawGrid() const;
};


#endif
