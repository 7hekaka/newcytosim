// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef LOCUS_GRID_H
#define LOCUS_GRID_H

#include "grid.h"
#include "dim.h"
#include "vector.h"
#include "modulo.h"
#include "mecapoint.h"
#include "fiber.h"
#include "fiber_segment.h"
#include "interpolation.h"
#include "array.h"

class Space;
class Modulo;
class Simul;
class Mecable;
class FiberSegment;

/// number of panes in the steric engine
/** This should normally be set equal to 1, for optimal performance */
#define MAX_STERIC_PANES 1


/// Used for early exclusion of potential pairs, representing { position, interaction radius }
/** The test is done in single precision arithmetics, which is faster and hopefully sufficient */
class BigVector
{
public:

    float xx, yy, zz;
    float rr;
    
    BigVector() { xx = 0; yy = 0; zz = 0; rr = 0; }

    BigVector(Vector1 v, real r) { xx = float(v.XX); yy = 0; zz = 0; rr = float(r); }
    BigVector(Vector2 v, real r) { xx = float(v.XX); yy = float(v.YY); zz = 0; rr = float(r); }
    BigVector(Vector3 v, real r) { xx = float(v.XX); yy = float(v.YY); zz = float(v.ZZ); rr = float(r); }

    //float const* data() const { return &XX; }
    operator float*() { return &xx; }
    operator const float*() const { return &xx; }

    /// @return result of test `distance(this, arg) < sum_of_ranges`
    bool near(BigVector const& arg) const
    {
        float x = xx - arg.xx;
#if ( DIM == 1 )
        float r = rr + arg.rr;
        return ( std::fabs(x) <= r );
#elif ( DIM == 2 )
        float y = yy - arg.yy;
        float r = rr + arg.rr;
        return ( x*x + y*y <= r*r );
#else
        float y = yy - arg.yy;
        float z = zz - arg.zz;
        float r = rr + arg.rr;
        return ( x*x + y*y < r*r - z*z );
#endif
    }
};


/// represents the Segment of a Fiber, or one vertex of a Mecable
class BigLocus
{
    friend class LocusGrid;
    
public:
    
    /// position of center, and radius of interaction
    BigVector pos_;

    /// Fiber containing the segment, or Mecable containing a vertex
    Mecable const* obj_;
    
    /// equilibrium radius of the interaction (distance where force is zero)
    float   rad_;

    /// index of segment's first point
    index_t vix_;
    
public:
    
    BigLocus() : obj_(nullptr), vix_(0) {}
    
    BigLocus(Mecable const* f, index_t i, real r, real e, Vector const& w)
    : pos_(w, e), obj_(f), rad_(float(r)), vix_(i)
    {
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
    
    /// a cast ot a fiber
    Fiber const* fiber() const { return static_cast<Fiber const*>(obj_); }
    
    /// length of segment
    real seg() const { return fiber()->segmentation(); }

    /// should return 1.0 / seg()
    real segInv() const { return fiber()->segmentationInv(); }
    
    /// true if abscissa 'a', counted from point 0 is within the segment
    bool within(real a) const { return ( 0 <= a ) & ( a <= fiber()->segmentation() ); }

    /// true if the segment is the first of the Fiber
    bool isFirst() const { return ( vix_ == 0 ); }
    
    /// true if the segment is the last of the Fiber
    bool isLast() const { return ( vix_+2 == obj_->nbPoints() ); }

    /// Mecapoint to point 1
    Mecapoint vertex1() const { return Mecapoint(obj_, vix_); }
    
    /// Mecapoint to point 2
    Mecapoint vertex2() const { return Mecapoint(obj_, vix_+1); }
    
    /// FiberSegment
    FiberSegment segment() const { return FiberSegment(fiber(), vix_); }
    
    /// return interpolation at distance 'abs' from vertex 1
    Interpolation interpolation(real abs) const { return Interpolation(obj_, abs*segInv(), vix_, 1+vix_); }

};

/// we used this alias for clarity and backward compatibility
typedef BigLocus BigPoint;

//------------------------------------------------------------------------------

/// a list containing BigLocus and BigPoint in the same Array
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
    size_t locuses_;
    
public:
    
    /// constructor
    BigLocusList()
    {
        locuses_ = 0;
    }
    
    /// clear all panes
    void clear()
    {
        pane.clear();
        locuses_ = 0;
    }
    
    /// number of elements in list
    size_t size() const { return pane.size(); }
    
    /// number of BigLocus in list
    size_t num_locus() const { return locuses_; }
    
    /// number of BigPoints in list
    size_t num_points() const { return pane.size() - locuses_; }

    /// first element in list
    BigLocus const* begin() const { return pane.begin(); }

    /// first BigPoint in list
    BigLocus const* middle() const { return pane.addr(locuses_); }
    
    /// one past last element in list
    BigLocus const* end() const { return pane.end(); }
    
    /// reference to Object at index ii (val_[ii])
    BigLocus const& operator[](const index_t i) const { return pane.at(i); }

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
        for ( index_t p = 1; p <= MAX_STERIC_PANES; ++p )
            panes[p].clear();
    }
    
    BigLocusList& cell_list(index_t p)
    {
        assert_true( 0 < p && p <= MAX_STERIC_PANES );
        return panes[p];
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
 This implements a divide-and-conquer method to find particles that are within
 a certain cutoff distance from each other. In brief:
 - It covers the space with a Grid `pGrid`, initialized by `setGrid()`
 - A list of class `LocusGridCell` is associated with each cell of this grid.
 - `LocusGrid::add()` links a `BigPoint` or a `BigLocus` to the appropriate cell of the grid.
 - `LocusGrid::setSterics()` checks all pairs of Point/Locus that may overlap,
    calculating their distance, and calling Meca::addLink() if they are interacting
 .
 LocusGrid is meant to be faster than PointGrid, but only supports repulsive interactions.
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
    
    /// Meca
    Meca& meca_;
    
    /// stiffness
    real push_;

private:
    
    /// check two Spheres
    void checkPP(BigPoint const&, BigPoint const&) const;
    
    /// check Sphere against Line segment
    void checkPL(BigPoint const&, BigLocus const&) const;
    
    /// check Line segment against vertex1 of other segment
    void checkLL1(BigLocus const&, BigLocus const&) const;

    /// check Line segment against the vertex2 of other segment
    void checkLL2(BigLocus const&, BigLocus const&) const;

    /// check vertex1 of segment against vertex1 of Segment
    void checkLLP1(BigLocus const&, BigLocus const&, float, real, Vector const&) const;
    
    /// check vertex1 and vertex2 of segment against vertex2 of Segment
    void checkLLP2(BigLocus const&, BigLocus const&, float, real, Vector const&) const;
    
    /// check two Line segments
    void checkLL(BigLocus const&, BigLocus const&) const;
    
    
    /// check all pairs between the two lists
    void setSterics0(BigLocusList const&) const;
    
    /// check all pairs between the two lists
    void setSterics0(BigLocusList const&, BigLocusList const&) const;
    
    /// check all pairs between the two lists, checking center-to-center distance
    void setStericsT(BigLocusList const&) const;
    
    /// check all pairs between the two lists, checking center-to-center distance
    void setStericsT(BigLocusList const&, BigLocusList const&) const;
    
    /// check all pairs between the two lists, checking center-to-center distance
    void setStericsX(BigLocusList const&) const;

    /// check all pairs between the two lists, checking center-to-center distance
    void setStericsX(BigLocusList const&, BigLocusList const&) const;

    /// used for steric with periodic boundary conditions
    inline index_t direct_index(Vector const& w)
    {
#if ( DIM == 3 )
        return pGrid.direct_index3D(w.XX, w.YY, w.ZZ);
#elif ( DIM == 2 )
        return pGrid.direct_index2D(w.XX, w.YY);
#else
        return pGrid.index(w);
#endif
    }
    
#if ( MAX_STERIC_PANES == 1 )
    
    /// cell corresponding to position `w`
    BigLocusList& cell_pane(Vector const& w) const
    {
        return pGrid.cell(w);
    }
    
    /// cell corresponding to index `w`
    BigLocusList& cell_pane(const index_t c) const
    {
        return pGrid.icell(c);
    }
    
    /// cell corresponding to index `w`
    BigLocusList& cell_list(const index_t w) const
    {
        return pGrid.icell(w);
    }
    
    /// enter interactions into Meca
    void setSterics0() const;
    
    /// enter interactions into Meca
    void setStericsT() const;

#else
    
    /// cell corresponding to position `w`, and pane `p`
    BigLocusList& cell_pane(Vector const& w, const index_t p) const
    {
        assert_true( 0 < p && p <= MAX_STERIC_PANES );
        return pGrid.cell(w).panes[p];
    }
    
    /// cell corresponding to index `c`, and pane `p`
    BigLocusList& cell_pane(const index_t c, const index_t p) const
    {
        assert_true( 0 < p && p <= MAX_STERIC_PANES );
        return pGrid.icell(c).panes[p];
    }
    
    /// cell corresponding to index `c`, and pane `p`
    BigLocusList& cell_list(const index_t c, const index_t p) const
    {
        assert_true( 0 < p && p <= MAX_STERIC_PANES );
        return pGrid.icell(c).panes[p];
    }
    
    /// enter interactions into Meca in one pane
    void setSterics0(index_t pan) const;
    
    /// enter interactions into Meca in one pane
    void setStericsT(index_t pan) const;
    
    /// enter interactions into Meca between two panes
    void setSterics0(index_t pan1, index_t pan2) const;
    
    /// enter interactions into Meca between two panes
    void setStericsT(index_t pan1, index_t pan2) const;
    
#endif
    
public:
    
    /// creator
    LocusGrid(Meca& m) : meca_(m), push_(0) {}
    
    /// number of panes
    static int nbPanes() { return MAX_STERIC_PANES; }
    
    /// set stiffness
    void stiffness(real s) { push_ = s; }
               
    /// define grid covering specified volume, given a minimal cell size requirement
    index_t setGrid(Vector inf, Vector sup, real min_width);
    
    /// allocate memory for grid
    void createCells();
    
    /// true if the grid was initialized by calling setGrid()
    size_t hasGrid() const { return pGrid.hasCells(); }

    /// sum of allocated size of lists for all cells
    size_t capacity() const;

    /// clear the grid
    void clear() { pGrid.clearCells(); }
    
    /// link in the cell containing the middle of the segment
    void add(Fiber const* fib, unsigned inx, real rad, real rge)
    {
        Vector w = fib->midPoint(inx);
#if ( MAX_STERIC_PANES <= 1 )
        BigLocusList& bll = cell_pane(w);
#else
        BigLocusList& bll = cell_pane(w, fib->prop->steric_key);
#endif
        bll.pane.emplace(fib, inx, rad, rge, w);
        bll.locuses_ += 1;
    }
    
    /// place Mecable vertex on the grid
    template <typename MECABLE>
    void add(MECABLE const* mec, unsigned inx, real rad)
    {
        Vector w = mec->posPoint(inx);
#if ( MAX_STERIC_PANES <= 1 )
        BigLocusList& bll = cell_pane(w);
#else
        BigLocusList& bll = cell_pane(w, mec->prop->steric_key);
#endif
        bll.pane.emplace(mec, inx, rad, rad, w);
    }
    
#if ENABLE_PERIODIC_BOUNDARIES
    /// link in the cell containing the middle of the segment:
    void add_modulo(Fiber const* fib, unsigned inx, real rad, real rge)
    {
        Vector w = fib->midPoint(inx);
        modulo->fold(w);
        index_t c = direct_index(w);
        //assert_true( c == pGrid.index(w) );
#if ( MAX_STERIC_PANES <= 1 )
        BigLocusList& bll = cell_pane(c);
#else
        BigLocusList& bll = cell_pane(c, fib->prop->steric_key);
#endif
        bll.pane.emplace(fib, inx, rad, rge, w);
        bll.locuses_ += 1;
    }

    /// place Mecable vertex on the grid
    template <typename MECABLE>
    void add_modulo(MECABLE const* mec, unsigned inx, real rad)
    {
        Vector w = mec->posPoint(inx);
        modulo->fold(w);
        index_t c = direct_index(w);
        //assert_true( c == pGrid.index(w) );
#if ( MAX_STERIC_PANES <= 1 )
        BigLocusList& bll = cell_pane(c);
#else
        BigLocusList& bll = cell_pane(c, mec->prop->steric_key);
#endif
        bll.pane.emplace(mec, inx, rad, rad, w);
    }
#endif
    
    /// enter interactions into Meca
    void setSterics() const;
    
#if ( MAX_STERIC_PANES > 1 )

    /// enter interactions into Meca in one pane
    void setSterics(index_t pan) const;
    
    /// enter interactions into Meca between two panes
    void setSterics(index_t pan1, index_t pan2) const;
    
#endif

    /// underlying spatial grid
    Map<DIM> const& map() const { return pGrid; }
    
    /// OpenGL display function
    void drawGrid() const;
};


#endif
