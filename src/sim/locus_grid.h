// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef LOCUS_GRID_H
#define LOCUS_GRID_H

#include "grid.h"
#include "dim.h"
#include "vector.h"
#include "array.h"
#include "fiber.h"

class Space;
class Modulo;
class Simul;
class Mecable;
class Mecapoint;
class FiberSegment;



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
        return ( z*z + x*x <= r*r - y*y );
#endif
    }
};



/// represents the point of a Mecable for steric interactions
class BigPoint
{
    friend class LocusGrid;
    
public:
    
    /// position of center
    BigVector pos_;
    
    /// Mecable containing the point-of-interest
    Mecable const* mec_;

    /// Index of the point-of-interest in the Mecable
    unsigned pti_;
    
    /// key to exclude certain pairs from interacting
    unsigned key_;
    
public:
    
    BigPoint() {}
    
    BigPoint(Mecable const* m, size_t i, real r, Vector const& w)
    : pos_(w, r)
    {
        mec_ = m;
        pti_ = static_cast<unsigned>(i);
        key_ = 0;
        assert_true( i == pti_ );
    }
    
    /// construct Mecapoint
    inline Mecapoint point() const;
    
    /// position of center
    Vector cen() const { return mec_->posPoint(pti_); }
    //Vector cen() const { return Vector(&pos_.XX); }
    
    /// radius
    real rad() const { return pos_.RR; }
};


/// represents the Segment of a Fiber for steric interactions
class BigLocus
{
    friend class LocusGrid;
    
public:
    
    /// position of center, and radius of interaction
    BigVector pos_;

    /// Fiber to which the segment belongs to
    Fiber const* fib_;
    
    /// equilibrium radius of the interaction (distance where force is zero)
    real     rad_;

    /// index of segment's first point
    unsigned pti_;
    
    /// key to exclude certain pairs from interacting
    unsigned key_;
    
public:
    
    BigLocus() {}
    
    BigLocus(Fiber const*& f, size_t i, real r, real e, Vector const& w)
    : pos_(w, e)
    {
        fib_ = f;
        rad_ = r;
        pti_ = static_cast<unsigned>(i);
        key_ = 0;
        assert_true( i == pti_ );
    }
    
    /// construct FiberSegment
    FiberSegment segment() const;

    /// position of point 1
    Vector pos1() const { return fib_->posPoint(pti_); }
    
    /// position of point 2
    Vector pos2() const { return fib_->posPoint(pti_+1); }
    
    /// offset = point2 - point1
    Vector diff() const { return fib_->diffPoints(pti_); }
    
    /// offset = point1 - point0
    Vector prevDiff() const { return fib_->diffPoints(pti_-1); }
    
    /// position of point 2
    real len() const { return fib_->segmentation(); }

    /// true if the segment is the first of the Fiber
    bool isFirst() const { return ( pti_ == 0 ); }
    
    /// true if the segment is the last of the Fiber
    bool isLast() const { return ( pti_+2 == fib_->nbPoints() ); }

    /// Mecapoint to point 1
    Mecapoint point1() const;
    
    /// Mecapoint to point 2
    Mecapoint point2() const;
    
    /// to point 1
    BigPoint bigPoint1() const { return BigPoint(fib_, pti_, rad_, pos1()); }
    
    /// to point 2
    BigPoint bigPoint2() const { return BigPoint(fib_, pti_+1, rad_, pos2()); }
};


/// type for a list of FatPoint
typedef Array<BigPoint> BigPointList;

/// type for a list of FatLocus
typedef Array<BigLocus> BigLocusList;


/// number of panes in the steric engine
/** This should normally be set equal to 1, for optimal performance */
#define MAX_STERIC_PANES 1


/// a set of lists associated with the same location
class LocusGridCell
{
    friend class LocusGrid;
    
#if ( MAX_STERIC_PANES == 1 )
    
    /// unique steric pane
    BigPointList point_pane;
    
    /// unique steric pane
    BigLocusList locus_pane;
    
#else
    
    /// different steric panes
    BigPointList point_panes_0[MAX_STERIC_PANES];
    
    /// different steric panes
    BigLocusList locus_panes_0[MAX_STERIC_PANES];
    
    /// alias to the array of panes, with index 1 refering to point_panes_0[0]
    BigPointList * point_panes;
    
    /// alias to the array of panes, with index 1 refering to locus_panes_0[0]
    BigLocusList * locus_panes;
    
#endif
    
public:
    
#if ( MAX_STERIC_PANES == 1 )
    
    LocusGridCell()
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
    
    LocusGridCell() : point_panes(point_panes_0), locus_panes(locus_panes_0)
    {
        --point_panes;
        --locus_panes;
    }
    
    /// clear all panes
    void clear()
    {
        for ( size_t p = 1; p <= MAX_STERIC_PANES; ++p )
        {
            point_panes[p].clear();
            locus_panes[p].clear();
        }
    }
    
    BigPointList& point_list(size_t p)
    {
        assert_true( 0 < p && p <= MAX_STERIC_PANES );
        return point_panes[p];
    }
    
    
    BigLocusList& locus_list(size_t p)
    {
        assert_true( 0 < p && p <= MAX_STERIC_PANES );
        return locus_panes[p];
    }
    
    
    size_t capacity() const
    {
        size_t res = 0;
        for ( int i = 0; i < MAX_STERIC_PANES; ++i )
            res += point_panes[i].capacity() + locus_panes[i].capacity();
        return res;
    }

#endif
};


/// Divide-and-Conquer to implement steric interactions
/**
 A divide-and-conquer algorithm is used to find FatPoints that overlap:
 - It uses a grid 'pGrid' covering the space, initialized by setGrid()
 To each point on pGrid is associated a list of FatPoint* of class PointGridCell.
 - The functions 'add()' position the given FatPoints on the grid
 - Function setStericInteraction() uses pGrid to find pairs of FatPoints that may overlap.
 It then calculates their actual distance, and set a interaction from Meca if necessary
 .
 */
class LocusGrid
{
private:
    
    /// grid for divide-and-conquer strategies:
    Grid<LocusGridCell, DIM> pGrid;

private:
    
#if ( MAX_STERIC_PANES == 1 )
    
    /// cell corresponding to position `w`, and pane `p`
    BigPointList& point_list(Vector const& w) const
    {
        return pGrid.cell(w).point_pane;
    }
    
    /// cell corresponding to position `w`, and pane `p`
    BigLocusList& locus_list(Vector const& w) const
    {
        return pGrid.cell(w).locus_pane;
    }
    
    /// cell corresponding to index `w`, and pane `p`
    BigPointList& point_list(const size_t w) const
    {
        return pGrid.icell(w).point_pane;
    }
    
    /// cell corresponding to index `w`, and pane `p`
    BigLocusList& locus_list(const size_t w) const
    {
        return pGrid.icell(w).locus_pane;
    }
    
#else
    
    /// cell corresponding to position `w`, and pane `p`
    BigPointList& point_list(Vector const& w, const size_t p) const
    {
        assert_true( 0 < p && p <= MAX_STERIC_PANES );
        return pGrid.cell(w).point_panes[p];
    }
    
    /// cell corresponding to position `w`, and pane `p`
    BigLocusList& locus_list(Vector const& w, const size_t p) const
    {
        assert_true( 0 < p && p <= MAX_STERIC_PANES );
        return pGrid.cell(w).locus_panes[p];
    }
    
    /// cell corresponding to index `c`, and pane `p`
    BigPointList& point_list(const size_t c, const size_t p) const
    {
        assert_true( 0 < p && p <= MAX_STERIC_PANES );
        return pGrid.icell(c).point_panes[p];
    }
    
    /// cell corresponding to index `c`, and pane `p`
    BigLocusList& locus_list(const size_t c, const size_t p) const
    {
        assert_true( 0 < p && p <= MAX_STERIC_PANES );
        return pGrid.icell(c).locus_panes[p];
    }
    
#endif
    
public:
    
    /// creator
    LocusGrid();
    
    /// define grid covering specified Space, with cell of size min_step at least
    size_t setGrid(Space const*, real min_step);
    
    /// allocate memory for grid
    void createCells();
    
    /// true if the grid was initialized by calling setGrid()
    size_t hasGrid() const  { return pGrid.hasCells(); }

    /// true if Grid has some periodic direction
    bool isPeriodic() const { return pGrid.isPeriodic(); }
    
    /// sum of allocated size of lists for all cells
    size_t capacity() const;

    /// clear the grid
    void clear() { pGrid.clear(); }
    
#if ( MAX_STERIC_PANES == 1 )
    
    /// place Mecable vertex on the grid
    void add(Mecable const* m, size_t i, real rad)
    {
        Vector w = m->posPoint(i);
        point_list(w).emplace(m, i, rad, w);
    }
    
    /// place Fiber segment on the grid
    void add(Fiber const* f, size_t i, real rad, real rge)
    {
        // link in cell containing the middle of the segment
        Vector w = f->posPoint(i, 0.5);
        locus_list(w).emplace(f, i, rad, rge, w);
    }
    
    /// enter interactions into Meca with given stiffness
    void setInteractions(Meca&, real stiff) const;
    
#else
    
    /// place Mecable vertex on the grid
    void add(size_t pane, Mecable const*, size_t, real rad);
    
    /// place Fiber segment on the grid
    void add(size_t pane, Fiber const*, size_t, real rad, real rge);
    
    /// enter interactions into Meca in one panes with given parameters
    void setInteractions(Meca&, real stiff, size_t pan) const;
    
    /// enter interactions into Meca between two panes with given parameters
    void setInteractions(Meca&, real stiff, size_t pan1, size_t pan2) const;
    
#endif
    
    /// underlying spatial grid
    Map<DIM> const& map() const { return pGrid; }
    
    /// OpenGL display function
    void draw() const;
};


#endif
