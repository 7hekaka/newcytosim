// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef LOCUS_GRID_H
#define LOCUS_GRID_H

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


/// represents the point of a Mecable for steric interactions
class BigPoint
{
    friend class LocusGrid;
    
public:
    
    /// position of center
    Vector   pos_;
    
    /// equilibrium radius of the interaction (distance where force is zero)
    real     rad_;
    
    /// Mecable containing the point-of-interest
    Mecable const* mec_;

    /// Index of the point-of-interest in the Mecable
    unsigned pti_;
    
    /// key to exclude certain pairs from interacting
    unsigned key_;
    
public:
    
    BigPoint() {}
    
    BigPoint(Mecable const* m, size_t i, real r, Vector const& w)
    {
        mec_ = m;
        pos_ = w;
        rad_ = r;
        pti_ = i;
        key_ = 0;
    }
    
    /// construct Mecapoint
    inline Mecapoint point() const { return Mecapoint(mec_, pti_); }
};


/// represents the Segment of a Fiber for steric interactions
class BigLocus
{
    friend class LocusGrid;
    
public:
    
    /// position of center
    Vector   pos_;
    
    /// equilibrium radius of the interaction (distance where force is zero)
    real     rad_;

    /// Fiber to which the segment belongs to
    Fiber const* fib_;
    
    /// index of segment's first point
    unsigned pti_;
    
    /// key to exclude certain pairs from interacting
    unsigned key_;
    
public:
    
    BigLocus() {}
    
    BigLocus(Fiber const*& f, size_t i, real r, Vector const& w)
    {
        fib_ = f;
        pos_ = w;
        rad_ = r;
        pti_ = i;
        key_ = 0;
    }
    
    /// construct FiberSegment
    FiberSegment segment() const { return FiberSegment(fib_, pti_); }

    /// position of point 1
    Vector pos1() const { return fib_->posPoint(pti_); }
    
    /// position of point 2
    Vector pos2() const { return fib_->posPoint(pti_+1); }
    
    /// offset = point2 - point1
    Vector diff() const { return fib_->diffPoints(pti_); }
    
    /// offset = point2 - point1
    Vector prevDiff() const { return fib_->diffPoints(pti_-1); }
    
    /// position of point 2
    real len() const { return fib_->segmentation(); }

    /// true if the segment is the first of the Fiber
    bool isFirst() const { return ( pti_ == 0 ); }
    
    /// true if the segment is the last of the Fiber
    bool isLast() const { return ( pti_+2 == fib_->nbPoints() ); }

    /// Mecapoint to point 1
    Mecapoint point1() const { return Mecapoint(fib_, pti_); }
    
    /// Mecapoint to point 2
    Mecapoint point2() const { return Mecapoint(fib_, pti_+1); }
    
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
    
    /// max radius that can be included
    real max_diameter;
    
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
    static void setInteractions(Meca&, real stiff,
                                BigPointList &, BigLocusList &);
    
    /// check all pairs between the two lists
    static void setInteractions(Meca&, real stiff,
                                BigPointList &, BigLocusList &,
                                BigPointList &, BigLocusList &);
    
    /// check all pairs between the two lists, checking center-to-center distance
    static void setInteractions(Meca&, real stiff, real sup,
                                BigPointList &, BigLocusList &);
    
    /// check all pairs between the two lists, checking center-to-center distance
    static void setInteractions(Meca&, real stiff, real sup,
                                BigPointList &, BigLocusList &,
                                BigPointList &, BigLocusList &);

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
    
    /// clear the grid
    void clear()            { pGrid.clear(); }
    
#if ( MAX_STERIC_PANES == 1 )
    
    /// place Mecapoint on the grid
    void add(Mecable const* m, size_t i, real rad) const
    {
        Vector w = m->posPoint(i);
        point_list(w).emplace(m, i, rad, w);
    }
    
    /// place FiberSegment on the grid
    void add(Fiber const* f, size_t i, real rad) const
    {
        // link in cell containing the middle of the segment
        Vector w = f->posPoint(i, 0.5);
        locus_list(w).emplace(f, i, rad, w);
    }
    
    /// enter interactions into Meca with given stiffness
    void setInteractions(Meca&, real stiff) const;
    
#else
    
    /// place Mecapoint on the grid
    void add(size_t pane, Mecable const*, size_t, real rad) const;
    
    /// place FiberSegment on the grid
    void add(size_t pane, Fiber const*, size_t, real rad) const;
    
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
