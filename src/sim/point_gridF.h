// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef POINT_GRIDF_H
#define POINT_GRIDF_H

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


/// represents a Mecapoint for steric interactions
class FatPointF
{
    friend class PointGridF;
    
public:
    
    /// current position of center
    Vector         pos;
    
    /// equilibrium radius of the interaction (distance where force is zero)
    real           radius;
    
    /// indicates the central vertex
    Mecapoint      pnt;
    
public:
    
    FatPointF() {}
    
    
    FatPointF(Mecapoint const& p, real rd, Vector const& w)
    {
        pnt    = p;
        radius = rd;
        pos    = w;
    }
    
    /// set from Mecapoint p, with radius=rd and range=rd+erg
    void set(Mecapoint const& p, real rd, Vector const& w)
    {
        pnt    = p;
        radius = rd;
        pos    = w;
    }
};


/// represents the Segment of a Fiber for steric interactions
class FatLocusF
{
    friend class PointGridF;
    
public:
    
    /// equilibrium radius of the interaction (distance where force is zero)
    real           radius;
    
    /// indicates the fiber's segment
    FiberSegment   seg;
    
public:
    
    FatLocusF() {}
    
    FatLocusF(FiberSegment const& p, real rd)
    {
        seg    = p;
        radius = rd;
    }
    
    /// set from FiberSegment p, with radius=rd and range=rd+erg
    void set(FiberSegment const& p, real rd)
    {
        seg    = p;
        radius = rd;
    }
    
    /// true if the segment is the first of the Fiber
    bool isFirst() const
    {
        return seg.isFirst();
    }
    
    /// true if the segment is the last of the Fiber
    bool isLast() const
    {
        return seg.isLast();
    }
    
    FatPointF point1() const
    {
        return FatPointF(seg.exact1(), radius, seg.pos1());
    }
    
    FatPointF point2() const
    {
        return FatPointF(seg.exact2(), radius, seg.pos2());
    }
};


/// type for a list of FatPoint
typedef Array<FatPointF> FatPointListF;

/// type for a list of FatLocus
typedef Array<FatLocusF> FatLocusListF;


/// number of panes in the steric engine
/** This should normally be set equal to 1, for optimal performance */
#define MAX_STERIC_PANES 1


/// a set of lists associated with the same location
class PointGridCellF
{
    friend class PointGridF;
    
#if ( MAX_STERIC_PANES == 1 )
    
    /// unique steric pane
    FatPointListF point_pane;
    
    /// unique steric pane
    FatLocusListF locus_pane;
    
#else
    
    /// different steric panes
    FatPointListF point_panes_0[MAX_STERIC_PANES];
    
    /// different steric panes
    FatLocusListF locus_panes_0[MAX_STERIC_PANES];
    
    /// alias to the array of panes, with index 1 refering to point_panes_0[0]
    FatPointListF * point_panes;
    
    /// alias to the array of panes, with index 1 refering to locus_panes_0[0]
    FatLocusListF * locus_panes;
    
#endif
    
public:
    
#if ( MAX_STERIC_PANES == 1 )
    
    PointGridCellF()
    {
    }
    
    /// clear all panes
    void clear()
    {
        point_pane.clear();
        locus_pane.clear();
    }
    
#else
    
    PointGridCellF() : point_panes(point_panes_0), locus_panes(locus_panes_0)
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
    
    FatPointListF& point_list(size_t p)
    {
        assert_true( 0 < p && p <= MAX_STERIC_PANES );
        return point_panes[p];
    }
    
    
    FatLocusListF& locus_list(size_t p)
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
class PointGridF
{
private:
    
    /// grid for divide-and-conquer strategies:
    Grid<PointGridCellF, DIM> pGrid;
    
    /// max radius that can be included
    real max_diameter;
    
private:
    
    /// check two Spheres
    void checkPP(Meca&, real stiff, FatPointF const&, FatPointF const&) const;
    
    /// check Sphere against Line segment
    void checkPL(Meca&, real stiff, FatPointF const&, FatLocusF const&) const;
    
    /// check Line segment against Sphere
    void checkLL1(Meca&, real stiff, FatLocusF const&, FatLocusF const&) const;
    
    /// check Line segment against Sphere
    void checkLL2(Meca&, real stiff, FatLocusF const&, FatLocusF const&) const;
    
    /// check two Line segments
    void checkLL(Meca&, real stiff, FatLocusF const&, FatLocusF const&) const;
    
    /// check all interacting pairs between the two lists
    void setInteractions(Meca&, real stiff,
                         FatPointListF &, FatLocusListF &) const;
    
    /// check all interacting pairs between the two lists
    void setInteractions(Meca&, real stiff,
                         FatPointListF &, FatLocusListF &,
                         FatPointListF &, FatLocusListF &) const;
    
#if ( MAX_STERIC_PANES == 1 )
    
    /// cell corresponding to position `w`, and pane `p`
    FatPointListF& point_list(Vector const& w) const
    {
        return pGrid.cell(w).point_pane;
    }
    
    /// cell corresponding to position `w`, and pane `p`
    FatLocusListF& locus_list(Vector const& w) const
    {
        return pGrid.cell(w).locus_pane;
    }
    
    /// cell corresponding to index `w`, and pane `p`
    FatPointListF& point_list(const size_t w) const
    {
        return pGrid.icell(w).point_pane;
    }
    
    /// cell corresponding to index `w`, and pane `p`
    FatLocusListF& locus_list(const size_t w) const
    {
        return pGrid.icell(w).locus_pane;
    }
    
#else
    
    /// cell corresponding to position `w`, and pane `p`
    FatPointListF& point_list(Vector const& w, const size_t p) const
    {
        assert_true( 0 < p && p <= MAX_STERIC_PANES );
        return pGrid.cell(w).point_panes[p];
    }
    
    /// cell corresponding to position `w`, and pane `p`
    FatLocusListF& locus_list(Vector const& w, const size_t p) const
    {
        assert_true( 0 < p && p <= MAX_STERIC_PANES );
        return pGrid.cell(w).locus_panes[p];
    }
    
    /// cell corresponding to index `c`, and pane `p`
    FatPointListF& point_list(const size_t c, const size_t p) const
    {
        assert_true( 0 < p && p <= MAX_STERIC_PANES );
        return pGrid.icell(c).point_panes[p];
    }
    
    /// cell corresponding to index `c`, and pane `p`
    FatLocusListF& locus_list(const size_t c, const size_t p) const
    {
        assert_true( 0 < p && p <= MAX_STERIC_PANES );
        return pGrid.icell(c).locus_panes[p];
    }
    
#endif
    
public:
    
    /// creator
    PointGridF();
    
    /// define grid covering specified Space, with cell of size min_step at least
    size_t setGrid(Space const*, real min_step);
    
    /// allocate memory for grid
    void createCells();
    
    /// true if the grid was initialized by calling setGrid()
    size_t hasGrid() const  { return pGrid.hasCells(); }
    
    /// clear the grid
    void clear()            { pGrid.clear(); }
    
#if ( MAX_STERIC_PANES == 1 )
    
    /// place Mecapoint on the grid
    void add(Mecapoint const& p, real radius) const
    {
        Vector w = p.pos();
        point_list(w).emplace_back(p, radius, w);
    }
    
    /// place FiberSegment on the grid
    void add(FiberSegment const& p, real radius) const
    {
        //we use the middle of the segment (interpolation coefficient is ignored)
        Vector w = p.center();
        locus_list(w).emplace_back(p, radius);
    }
    
    /// enter interactions into Meca with given stiffness
    void setInteractions(Meca&, real stiff) const;
    
#else
    
    /// place Mecapoint on the grid
    void add(size_t pane, Mecapoint const&, real radius) const;
    
    /// place FiberSegment on the grid
    void add(size_t pane, FiberSegment const&, real radius) const;
    
    /// enter interactions into Meca in one panes with given parameters
    void setInteractions(Meca&, real stiff, size_t pan) const;
    
    /// enter interactions into Meca between two panes with given parameters
    void setInteractions(Meca&, real stiff, size_t pan1, size_t pan2) const;
    
#endif
    
    /// OpenGL display function
    void draw() const;
};


#endif
