// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
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


/// Stores the stiffness parameters for the steric engine
class Stiffness
{
public:
    real push;
    real pull;
    
    Stiffness(real h, real l)
    {
        push = h;
        pull = l;
    }
};


/// Used for early exclusing of potential pairs, representing { position, interaction radius }
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
    bool isFirst() const
    {
        return seg_.isFirst();
    }
    
    /// true if the segment is the last of the Fiber
    bool isLast() const
    {
        return seg_.isLast();
    }
    
    FatPoint fatPoint1() const
    {
        return FatPoint(seg_.exact1(), rad_, rge_, seg_.pos1());
    }
    
    FatPoint fatPoint2() const
    {
        return FatPoint(seg_.exact2(), rad_, rge_, seg_.pos2());
    }
    
    /// offset = point1 - point0
    Vector prevDiff() const { return seg_.fiber()->diffPoints(seg_.point()-1); }
};


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
    
private:
    
    /// check two Spheres
    static void checkPP(Meca&, Stiffness const&, FatPoint const&, FatPoint const&);
    
    /// check Sphere against Line segment
    static void checkPL(Meca&, Stiffness const&, FatPoint const&, FatLocus const&);
    
    /// check Line segment against Sphere
    static void checkLL1(Meca&, Stiffness const&, FatLocus const&, FatLocus const&);
    
    /// check Line segment against the terminal Sphere of a Fiber
    static void checkLL2(Meca&, Stiffness const&, FatLocus const&, FatLocus const&);
    
    /// check two Line segments
    static void checkLL(Meca&, Stiffness const&, FatLocus const&, FatLocus const&);
    
    /// check all pairs between the two lists
    static void setSterics0(Meca&, Stiffness const&,
                            FatPointList &, FatLocusList &);
    
    /// check all pairs between the two lists
    static void setSterics0(Meca&, Stiffness const&,
                            FatPointList &, FatLocusList &,
                            FatPointList &, FatLocusList &);
    
    /// check all pairs between two lists, checking center-to-center distance
    static void setStericsT(Meca&, Stiffness const&,
                            FatPointList &, FatLocusList &);
    
    /// check all pairs between two lists, checking center-to-center distance
    static void setStericsT(Meca&, Stiffness const&,
                            FatPointList &, FatLocusList &,
                            FatPointList &, FatLocusList &);

#if ( NUM_STERIC_PANES == 1 )
    
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
    FatPointList& point_list(const size_t w) const
    {
        return pGrid.icell(w).point_pane;
    }
    
    /// cell corresponding to index `w`, and pane `p`
    FatLocusList& locus_list(const size_t w) const
    {
        return pGrid.icell(w).locus_pane;
    }
    
    /// enter interactions into Meca with given stiffness
    void setSterics0(Meca&, Stiffness const&) const;
    
    /// enter interactions into Meca with given stiffness
    void setStericsT(Meca&, Stiffness const&) const;

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
    
    /// enter interactions into Meca in one panes with given parameters
    void setSterics0(Meca&, Stiffness const&, size_t pan) const;
    
    /// enter interactions into Meca in one panes with given parameters
    void setStericsT(Meca&, Stiffness const&, size_t pan) const;
    
    /// enter interactions into Meca between two panes with given parameters
    void setSterics0(Meca&, Stiffness const&, size_t pan1, size_t pan2) const;
    
    /// enter interactions into Meca between two panes with given parameters
    void setStericsT(Meca&, Stiffness const&, size_t pan1, size_t pan2) const;

#endif
    
public:
    
    /// creator
    PointGrid();
    
    /// define grid covering specified Space, with cell of size min_step at least
    size_t setGrid(Space const*, real min_step);
    
    /// allocate memory for grid
    void createCells();
    
    /// true if the grid was initialized by calling setGrid()
    size_t hasGrid() const { return pGrid.hasCells(); }
    
    /// sum of allocated size of lists for all cells
    size_t capacity() const;

    /// clear the grid
    void clear() { pGrid.clear(); }
    
#if ( NUM_STERIC_PANES == 1 )
    
    /// place Mecapoint on the grid
    void add(Mecable const* m, size_t i, real rad, real rge) const
    {
        Vector w = m->posPoint(i);
        point_list(w).emplace(Mecapoint(m, i), rad, rge, w);
    }
    
    /// place FiberSegment on the grid
    void add(Fiber const* f, size_t i, real rad, real rge, real sup) const
    {
        // link in cell containing the middle of the segment
        Vector w = f->posPoint(i, 0.5);
        locus_list(w).emplace(FiberSegment(f, i), rad, rge, sup, w);
    }
    
    /// enter interactions into Meca with given stiffness
    void setInteractions(Meca&, Stiffness const&) const;
    
#else
    
    /// place Mecapoint on the grid
    void add(size_t pane, Mecable const*, size_t, real rad, real rge) const;
    
    /// place FiberSegment on the grid
    void add(size_t pane, Fiber const*, size_t, real rad, real rge, real sup) const;
    
    /// enter interactions into Meca in one panes with given parameters
    void setInteractions(Meca&, Stiffness const&, size_t pan) const;
    
    /// enter interactions into Meca between two panes with given parameters
    void setInteractions(Meca&, Stiffness const&, size_t pan1, size_t pan2) const;
    
#endif
    
    /// underlying spatial grid
    Map<DIM> const& map() const { return pGrid; }
    
    /// OpenGL display function
    void draw() const;
};


#endif
