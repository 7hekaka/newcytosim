// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef POINT_GRID_H
#define POINT_GRID_H

#include "dim.h"
#include "vector.h"
#include "grid.h"
#include "mecapoint.h"
#include "fiber_segment.h"
#include "array.h"

#ifdef DISPLAY
#include "grid_display.h"
#endif

class Space;
class Modulo;
class Simul;
class Fiber;


/// represents a Mecapoint for steric interactions
class FatPoint
{
    friend class PointGrid;
    
public:
    
    /// buffer for position
    Vector         pos;
    
    /// equilibrium radius of the interaction (distance where force is zero)
    real           radius;
    
    /// interaction range (maximum distance at which the force can operate)
    real           range;
    
    /// indicates the central vertex
    Mecapoint      pnt;
    
public:
    
    FatPoint() {}
    
    
    FatPoint(Mecapoint const& p, real rd, real rg, Vector const& w)
    {
        pnt    = p;
        radius = rd;
        range  = rg;
        pos    = w;
    }
    
    /// set from Mecapoint p, with radius=rd and range=rd+erg
    void set(Mecapoint const& p, real rd, real rg, Vector const& w)
    {
        pnt    = p;
        radius = rd;
        range  = rg;
        pos    = w;
    }
};


/// represents the Segment of a Fiber for steric interactions
class FatLocus
{
    friend class PointGrid;
    
public:
    
    /// equilibrium radius of the interaction (distance where force is zero)
    real           radius;
    
    /// interaction range (maximum distance at which the force can operate)
    real           range;
    
    /// this FatLocus represents the entire segment indicated by the FiberSegment
    FiberSegment   seg;
    
public:
    
    FatLocus() {}
    
    FatLocus(FiberSegment const& p, real rd, real rg)
    {
        seg    = p;
        radius = rd;
        range  = rg;
    }
    
    /// set from FiberSegment p, with radius=rd and range=rd+erg
    void set(FiberSegment const& p, real rd, real rg)
    {
        seg    = p;
        radius = rd;
        range  = rg;
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
    
    FatPoint point1() const
    {
        return FatPoint(seg.exact1(), radius, range, seg.pos1());
    }
    
    FatPoint point2() const
    {
        return FatPoint(seg.exact2(), radius, range, seg.pos2());
    }
};


/// type for a list of FatPoint
typedef Array<FatPoint> FatPointList;

/// type for a list of FatLocus
typedef Array<FatLocus> FatLocusList;


/// number of panes in the steric engine
/** This should normally be set equal to 1, for optimal performance */
#define NB_STERIC_PANES 3


/// a set of lists associated with the same location
class PointGridCell
{
    friend class PointGrid;
    
#if ( NB_STERIC_PANES == 1 )
    
    /// unique steric pane
    FatPointList point_pane;
    
    /// unique steric pane
    FatLocusList locus_pane;
    
#else
    
    /// different steric panes
    FatPointList point_panes_0[NB_STERIC_PANES];
    
    /// different steric panes
    FatLocusList locus_panes_0[NB_STERIC_PANES];
    
    /// alias to the array of panes, with index 1 refering to point_panes_0[0]
    FatPointList * point_panes;
    
    /// alias to the array of panes, with index 1 refering to locus_panes_0[0]
    FatLocusList * locus_panes;
    
#endif
    
public:
    
#if ( NB_STERIC_PANES == 1 )
    
    PointGridCell()
    {
    }
    
    /// clear all panes
    void clear()
    {
        point_pane.clear();
        locus_pane.clear();
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
        for ( unsigned p = 1; p <= NB_STERIC_PANES; ++p )
        {
            point_panes[p].clear();
            locus_panes[p].clear();
        }
    }
    
    FatPointList& point_list(unsigned p)
    {
        assert_true( 0 < p && p <= NB_STERIC_PANES );
        return point_panes[p];
    }
    
    
    FatLocusList& locus_list(unsigned p)
    {
        assert_true( 0 < p && p <= NB_STERIC_PANES );
        return locus_panes[p];
    }
    
#endif
};


/// Contains the stiffness parameters for the steric engine
class StericParam
{
public:
    real stiff_push;
    real stiff_pull;
    
    StericParam(real push, real pull)
    {
        stiff_push = push;
        stiff_pull = pull;
    }
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
class PointGrid
{
private:
    
    /// grid for divide-and-conquer strategies:
    Grid<PointGridCell, DIM> pGrid;
    
    /// max radius that can be included
    real max_diameter;
    
private:
    
    /// check two Spheres
    void checkPP(Meca&, StericParam const&, FatPoint const&, FatPoint const&) const;
    
    /// check Sphere against Line segment
    void checkPL(Meca&, StericParam const&, FatPoint const&, FatLocus const&) const;
    
    /// check Line segment against Sphere
    void checkLL1(Meca&, StericParam const&, FatLocus const&, FatLocus const&) const;
    
    /// check Line segment against Sphere
    void checkLL2(Meca&, StericParam const&, FatLocus const&, FatLocus const&) const;
    
    /// check two Line segments
    void checkLL(Meca&, StericParam const&, FatLocus const&, FatLocus const&) const;
    
    /// check all interacting pairs between the two lists
    void setInteractions(Meca&, StericParam const&,
                         FatPointList &, FatLocusList &) const;
    
    /// check all interacting pairs between the two lists
    void setInteractions(Meca&, StericParam const&,
                         FatPointList &, FatLocusList &,
                         FatPointList &, FatLocusList &) const;
    
#if ( NB_STERIC_PANES == 1 )
    
    /// cell corresponding to position `w`, and pane `p`
    FatPointList& point_list(Vector const& w) const
    {
        return pGrid.cell(w).point_pane;
    }
    
    /// cell corresponding to position `w`, and pane `p`
    FatLocusList& locus_list(Vector const& w) const
    {
        return pGrid.cell(w).locus_pane;
    }
    
    /// cell corresponding to index `w`, and pane `p`
    FatPointList& point_list(const unsigned w) const
    {
        return pGrid.icell(w).point_pane;
    }
    
    /// cell corresponding to index `w`, and pane `p`
    FatLocusList& locus_list(const unsigned w) const
    {
        return pGrid.icell(w).locus_pane;
    }
    
#else
    
    /// cell corresponding to position `w`, and pane `p`
    FatPointList& point_list(Vector const& w, const unsigned p) const
    {
        assert_true( 0 < p && p <= NB_STERIC_PANES );
        return pGrid.cell(w).point_panes[p];
    }
    
    /// cell corresponding to position `w`, and pane `p`
    FatLocusList& locus_list(Vector const& w, const unsigned p) const
    {
        assert_true( 0 < p && p <= NB_STERIC_PANES );
        return pGrid.cell(w).locus_panes[p];
    }
    
    /// cell corresponding to index `c`, and pane `p`
    FatPointList& point_list(const unsigned c, const unsigned p) const
    {
        assert_true( 0 < p && p <= NB_STERIC_PANES );
        return pGrid.icell(c).point_panes[p];
    }
    
    /// cell corresponding to index `c`, and pane `p`
    FatLocusList& locus_list(const unsigned c, const unsigned p) const
    {
        assert_true( 0 < p && p <= NB_STERIC_PANES );
        return pGrid.icell(c).locus_panes[p];
    }
    
#endif
    
public:
    
    /// creator
    PointGrid();
    
    /// define grid covering specified Space, with cell of size min_step at least
    size_t setGrid(Space const*, real min_step);
    
    /// allocate memory for grid
    void createCells();
    
    /// true if the grid was initialized by calling setGrid()
    size_t hasGrid() const  { return pGrid.hasCells(); }
    
    /// clear the grid
    void clear()            { pGrid.clear(); }
    
#if ( NB_STERIC_PANES == 1 )
    
    /// place Mecapoint on the grid
    void add(Mecapoint const& p, real radius, real extra_range) const
    {
        Vector w = p.pos();
        point_list(w).new_val().set(p, radius, extra_range, w);
    }
    
    /// place FiberSegment on the grid
    void add(FiberSegment const& p, real radius, real extra_range) const
    {
        //we use the middle of the segment (interpolation coefficient is ignored)
        Vector w = p.center();
        locus_list(w).new_val().set(p, radius, extra_range);
    }
    
    /// enter interactions into Meca with given stiffness
    void setInteractions(Meca&, StericParam const& pam) const;
    
#else
    
    /// place Mecapoint on the grid
    void add(unsigned pane, Mecapoint const&, real radius, real extra_range) const;
    
    /// place FiberSegment on the grid
    void add(unsigned pane, FiberSegment const&, real radius, real extra_range) const;
    
    /// enter interactions into Meca in one panes with given parameters
    void setInteractions(Meca&, StericParam const& pam, unsigned pan) const;
    
    /// enter interactions into Meca between two panes with given parameters
    void setInteractions(Meca&, StericParam const& pam, unsigned pan1, unsigned pan2) const;
    
#endif
    
#ifdef DISPLAY
    void draw() const
    {
        glPushAttrib(GL_LIGHTING_BIT);
        glDisable(GL_LIGHTING);
        drawEdges(pGrid);
        glPopAttrib();
    }
#endif
};


#endif
