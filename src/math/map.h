// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Francois Nedelec; Created 07/03/2015. nedelec@embl.de

#ifndef MAP_H
#define MAP_H

#include "assert_macro.h"
#include "exceptions.h"
#include <cstdio>
#include <cmath>
#include "real.h"

///\def compile switch to disable support for periodic boundaries
/**
 Periodic boundaries are normally supported using a if statement called repeatedly.
 However, if GRID_HAS_PERIODIC==0, this test is ommitted, which might be faster,
 but Periodic boundaries cannot be supported henceforth.
 */
#define GRID_HAS_PERIODIC 1


/// Map divides a rectangle of dimensionality ORD into regular voxels
/** 
Map<ORD>, where ORD is an integer creates a regular grid over a rectangular region
of space of dimensionality ORD, initialized by setDimensions().

Functions are provided to convert from the space coordinates (of type real)
into an index usable to access the one-dimensional C-array of cells.
The cells are ordered successively, the first dimension (X) varying the fastest
i.e. cell[ii+1] will in most cases be located on the right of cell[ii], although
if cell[ii] is on the right edge, then cell[ii+1] is on the symmetric edge. 

\par Access:

Cells can be accessed in three ways:
 - Position:      a set of real       operator()( real[] ), or operator(real, real, real)
 - Index:         one integer         operator[](int index)
 - Coordinates:   a set of integer    function cell(int[]), or cell(int,int,int)
.
Valid indices are [0...nbCells()-1], where nbCells() is calculated by setDimensions().
If a position lies outside the rectangular region where the grid is defined,
index(real[]) returns the index of the closest voxel.

Functions to convert between the three types are provided:
 - index()
 - pack()
 - setCoordinatesFromIndex(),
 - setCoordinatesFromPosition()
 - setPositionFromCoordinates()
 - setPositionFromIndex()
.

\par Indices:

The grid is initialized by setDimensions(inf, sup, nbCells), which calculates:
  - cWidth[d] = ( sup[d] - inf[d] ) / nbCells[d], for d in [0, ORD[

The coordinates of a cell at position pos[] are:
  - c[d] = int(  ( pos[d] - inf[d] ) / cWidth[d] )

and its index is
  - with ORD==1: index = c[0]
  - with ORD==2: index = c[0] + nbcells[0] * c[1]
  - with ORD==3: index = c[0] + nbcells[0] * ( c[1] + nbcells[1] * c[2] )
  - etc.
.
    
For a 4x4 2D grid, the index are like this:

    12  13  14  15
    8    9  10  11
    4    5   6   7
    0    1   2   3

\par Neighborhood:

The class also provides information on which cells surround each cell:
 - createSquareRegions(range) calculates square regions of size range
   ( range==1 gives nearest neighbors ).
 - createRoundRegions(range) calculates round regions of size range
 - createSideRegions(range)
.
After calling one of the above function, getRegion(offsets, index) will set 'offsets'
to point to an array of 'index offsets' for the cell referred by 'index'.
A zero offset value (0) is always first in the list and refers to self.
In the example above:
    - for index = 0 it would return { 0 1 4 5 }
    - for index = 5 it would return { 0 -1 1 -5 -4 -3 3 4 5 }
.
You obtain the cell-indices of the neighboring cells by adding offsets[n] to 'index':
Example:

    CELL * cell = & map.icell(indx);
    int n_neighbors = map.getRegion(region, indx);
    for ( int n = 1; n < n_neighbors; ++n )
    {
        Cell & neighbor = cell[region[n]];
        ...
    }

*/

///\todo add Map<> copy constructor and copy assignment

template <int ORD>
class Map
{
public:

    /// Disabled copy constructor
    Map<ORD>(Map<ORD> const&);
    
    /// Disabled copy assignment
    Map<ORD>& operator=(Map<ORD> const&);

protected:
   
    /// Total number of cells in the map; size of cells[]
    size_t  mNbCells;
    
    /// The number of cells in each dimension
    size_t  mDim[ORD];
    
    /// Offset between two consecutive cells along each dimension
    size_t  mStride[ORD];
    
    /// The position of the inferior (min) edge in each dimension
    real    mInf[ORD];
    
    /// The position of the superior (max) edge in each dimension
    real    mSup[ORD];

    /// The size of a cell: cWidth[d] = ( mSup[d] - inf[d] ) / mDim[d]
    real    cWidth[ORD];
    
    /// cDelta[d] = 1.0 / cWidth[d]
    real    cDelta[ORD];
    
    /// mStart[d] = mInf[d] / cWidth[d]
    real    mStart[ORD];

    /// The volume of one cell
    real    cVolume;
    
    /// true if Map has periodic boundary conditions
    bool    mPeriodic[ORD];

protected:
    
    /// return closest integer to `c` in the segment [ 0, s-1 ]
    inline static size_t imagei_periodic(size_t s, long c)
    {
        ///@todo use remainder() function for branchless code?
        while ( c <  0 ) c += s;
        size_t u = (size_t) c;
        while ( u >= s ) u -= s;
        return u;
    }

    inline static size_t imagei_clamped(size_t s, long c)
    {
        return std::min((size_t)std::max(0L, c), s-1);
        //return c <= 0 ? 0 : ( c >= s ? s-1 : c );
    }
    
    /// return closest integer to `c` in the segment [ 0, mDim[d]-1 ]
    inline size_t image(const int d, long c) const
    {
#if GRID_HAS_PERIODIC
        if ( mPeriodic[d] )
            return imagei_periodic(mDim[d], c);
        else
#endif
            return imagei_clamped(mDim[d], c);
    }


    /// return f modulo s in [ 0, s-1 ]
    inline static size_t imagef_periodic(size_t s, real f)
    {
        while ( f <  0 )  f += (real)s;
        size_t u = (size_t) f;
        while ( u >= s )  u -= s;
        return u;
    }

    inline static size_t imagef_clamped(size_t s, real f)
    {
        if ( f > 0 )
        {
            size_t u = (size_t) f;
            //return ( u >= s ? s-1 : u );
            return std::min(u, s-1);
        }
        return 0;
    }
    
    
    /// returns  ( f - mInf[d] ) / cWidth[d]
    inline real map(const int d, real f) const
    {
        return f * cDelta[d] - mStart[d];
    }

    /// return closest integer to `c` in the segment [ 0, mDim[d]-1 ]
    inline size_t imagef(const int d, real f) const
    {
        real x = map(d, f);
#if GRID_HAS_PERIODIC
        if ( mPeriodic[d] )
            return imagef_periodic(mDim[d], x);
        else
#endif
            return imagef_clamped(mDim[d], x);
    }

    /// return closest integer to `c` in the segment [ 0, mDim[d]-1 ]
    inline size_t imagef(const int d, real f, real offset) const
    {
        real x = map(d, f) + offset;
#if GRID_HAS_PERIODIC
        if ( mPeriodic[d] )
            return imagef_periodic(mDim[d], x);
        else
#endif
            return imagef_clamped(mDim[d], x);
    }

//--------------------------------------------------------------------------
#pragma mark -
public:
    
    /// constructor
    Map() : mDim{0}, mInf{0}, mSup{0}, cWidth{0}, cDelta{0}, mStart{0}, mPeriodic{false}
    {
        mNbCells    = 0;
        regionsEdge = nullptr;
        regions     = nullptr;
        cVolume     = 0;
    }
    
    /// Free memory
    void destroy()
    {
        deleteRegions();
    }
    
    /// Destructor
    virtual ~Map()
    {
        destroy();
    }
    
    //--------------------------------------------------------------------------
    /// specifies the area covered by the Grid
    /**
     the edges of the area are specified in dimension `d` by 'infs[d]' and 'sups[d]',
     and the number of cells by 'nbcells[d]'.
     */
    void setDimensions(const real infs[ORD], real sups[ORD], const size_t cells[ORD])
    {
        cVolume = 1;
        size_t cnt = 1;

        for ( int d = 0; d < ORD; ++d )
        {
            if ( cells[d] <= 0 )
                throw InvalidParameter("Cannot build grid as nbcells[] is <= 0");
            
            if ( infs[d] >= sups[d] )
                throw InvalidParameter("Cannot build grid as sup[] <= inf[]");
            
            mStride[d] = cnt;
            cnt      *= cells[d];
            mDim[d]   = cells[d];
            mInf[d]   = infs[d];
            mSup[d]   = sups[d];
            cWidth[d] = ( mSup[d] - mInf[d] ) / real(mDim[d]);
            // inverse of cell width:
            cDelta[d] = real(mDim[d]) / ( mSup[d] - mInf[d] );
            mStart[d] = ( mDim[d] * mInf[d] ) / ( mSup[d] - mInf[d] );
            cVolume  *= cWidth[d];
        }
        mNbCells = cnt;
    }
    
    ///true if setDimensions() was called
    bool hasDimensions() const
    {
        return mNbCells > 0;
    }
    
    /// true if dimension `d` has periodic boundary conditions
    bool isPeriodic(int d) const
    {
#if GRID_HAS_PERIODIC
        if ( d < ORD )
            return mPeriodic[d];
#endif
        return false;
    }
    
    /// change boundary conditions
    void setPeriodic(size_t d, bool p)
    {
#if GRID_HAS_PERIODIC
        if ( d < ORD )
            mPeriodic[d] = p;
#else
        if ( p )
            throw InvalidParameter("grid.h was compiled with GRID_HAS_PERIODIC = 0");
#endif
    }
    
    /// true if boundary conditions are periodic
    bool isPeriodic() const
    {
#if GRID_HAS_PERIODIC
        for ( int d = 0; d < ORD; ++d )
            if ( mPeriodic[d] )
                return true;
#endif
        return false;
    }

    //--------------------------------------------------------------------------
#pragma mark -

    /// total number of cells in the map
    size_t       nbCells()           const { return mNbCells; }

    /// number of cells in dimensionality `d`
    size_t       breadth(size_t d)   const { return mDim[d]; }
    
    /// offset to the next cell in the direction `d`
    size_t       stride(size_t d)    const { return mStride[d]; }
    
    /// position of the inferior (left/bottom/etc) edge
    const real*  inf()               const { return mInf;    }
    real         inf(size_t d)       const { return mInf[d]; }
    
    /// position of the superior (right/top/etc) edge
    const real*  sup()               const { return mSup;    }
    real         sup(size_t d)       const { return mSup[d]; }
    
    /// the widths of a cell
    const real*  delta()             const { return cDelta;    }
    real         delta(size_t d)     const { return cDelta[d]; }
    
    const real*  cellWidth()         const { return cWidth;    }
    real         cellWidth(size_t d) const { return cWidth[d]; }
    
    /// the volume of a cell
    real         cellVolume()        const { return cVolume; }

    /// position in dimension `d`, of the cell of index `c`
    real         position(size_t d, real c) const { return mInf[d] + c * cWidth[d]; }
    
    /// index in dimension `d` corresponding to position `w`
    int          index(size_t d, real w) const { return (int)map(d, w); }

    /// half the diagonal length of the unit cell
    real cellRadius() const
    {
        real res = cWidth[0] * cWidth[0];
        for ( int d = 1; d < ORD; ++d )
            res += cWidth[d] * cWidth[d];
        return 0.5 * std::sqrt(res);
    }
    
    /// the smallest cell width, along dimensions that have more than `min_size` cells
    real minimumWidth(size_t min_size) const
    {
        real res = INFINITY;
        for ( int d = 0; d < ORD; ++d )
        {
            if ( mDim[d] > min_size )
                res = std::min(res, cWidth[d]);
        }
        return res;
    }
    
    /// radius of the minimal sphere placed in (0,0,0) that entirely covers all cells
    real radius() const
    {
        real res = 0;
        for ( int d = 0; d < ORD; ++d )
        {
            real m = std::max(mSup[d], -mInf[d]);
            res += m * m;
        }
        return std::sqrt(res);
    }

    //--------------------------------------------------------------------------
#pragma mark - Conversion

    /// checks if coordinates are inside the box
    bool inside(const int coord[ORD]) const
    {
        for ( int d = 0; d < ORD; ++d )
        {
            if ( coord[d] < 0 || (size_t)coord[d] >= mDim[d] )
                return false;
        }
        return true;
    }
    
    /// checks if point is inside the box
    bool inside(const real w[ORD]) const
    {
        for ( int d = 0; d < ORD; ++d )
        {
            if ( w[d] < mInf[d] || w[d] >= mSup[d] )
                return false;
        }
        return true;
    }
    
    /// periodic image
    void bringInside(int coord[ORD]) const
    {
        for ( int d = 0; d < ORD; ++d )
            coord[d] = image(d, coord[d]);
    }
    
    /// conversion from index to coordinates
    void setCoordinatesFromIndex(int coord[ORD], size_t indx) const
    {
        for ( int d = 0; d < ORD; ++d )
        {
            coord[d] = indx % mDim[d];
            indx    /= mDim[d];
        }
    }
    
    /// conversion from Position to coordinates (offset should be in [0,1])
    void setCoordinatesFromPosition(int coord[ORD], const real w[ORD], const real offset=0) const
    {
        for ( int d = 0; d < ORD; ++d )
            coord[d] = imagef(d, w[d], offset);
    }

    /// conversion from Index to Position (offset should be in [0,1])
    void setPositionFromIndex(real w[ORD], size_t indx, real offset) const
    {
        for ( int d = 0; d < ORD; ++d )
        {
            w[d] = mInf[d] + cWidth[d] * ( offset + indx % mDim[d] );
            indx /= mDim[d];
        }
    }
    
    /// conversion from Coordinates to Position (offset should be in [0,1])
    void setPositionFromCoordinates(real w[ORD], const int coord[ORD], real offset=0) const
    {
        for ( int d = 0; d < ORD; ++d )
            w[d] = mInf[d] + cWidth[d] * ( offset + coord[d] );
    }

    /// conversion from coordinates to index
    size_t pack(const int coord[ORD]) const
    {
        size_t inx = image(ORD-1, coord[ORD-1]);
        
        for ( int d = ORD-2; d >= 0; --d )
            inx = mDim[d] * inx + image(d, coord[d]);
        
        return inx;
    }
    
    
    /// returns the index of the cell whose center is closest to the point w[]
    size_t index(const real w[ORD]) const
    {
        size_t inx = imagef(ORD-1, w[ORD-1]);
        
        for ( int d = ORD-2; d >= 0; --d )
            inx = mDim[d] * inx + imagef(d, w[d]);
        
        return inx;
    }

    
    /// returns the index of the cell whose center is closest to the point w[]
    size_t index(const real w[ORD], const real offset) const
    {
        size_t inx = imagef(ORD-1, w[ORD-1], offset);
        
        for ( int d = ORD-2; d >= 0; --d )
            inx = mDim[d] * inx + imagef(d, w[d], offset);
        
        return inx;
    }

    
    /// return cell that is next to `c` in the direction `dim`
    size_t next(size_t c, int dim) const
    {
        size_t s[ORD];
        for ( int d = 0; d < ORD; ++d )
        {
            s[d] = c % mDim[d];
            c   /= mDim[d];
        }

        s[dim] = image(dim, s[dim]+1);

        c = s[ORD-1];
        for ( int d = ORD-2; d >= 0; --d )
            c = mDim[d] * c + s[d];
        return c;
    }
    
    /// convert coordinate to array index, if ORD==1
    size_t pack1D(const int x) const
    {
        return image(0, x);
    }
    
    /// convert coordinate to array index, if ORD==2
    size_t pack2D(const int x, const int y) const
    {
        return image(0, x) + mDim[0]*image(1, y);
    }
    
    /// convert coordinate to array index, if ORD==3
    size_t pack3D(const int x, const int y, const int z) const
    {
        return image(0, x) + mDim[0]*( image(1, y) + mDim[1]*image(2, z) );
    }

    
    /// return index of cell corresponding to position (x), if ORD==1
    size_t index1D(const real x) const
    {
        return image(0, map(0, x));
    }
    
    /// return index of cell corresponding to position (x, y), if ORD==2
    size_t index2D(const real x, const real y) const
    {
        return image(0, map(0, x))
               + mDim[0] * image(1, map(1, y));
    }
    
    /// return index of cell corresponding to position (x, y, z), if ORD==3
    size_t index3D(const real x, const real y, const real z) const
    {
        return image(0, map(0, x))
               + mDim[0]*( image(1, map(1, y))
               + mDim[1]*  image(2, map(2, z)) );
    }

    //--------------------------------------------------------------------------

#pragma mark - Regions

    /** For any cell, we can find the adjacent cells by adding 'index offsets'
    However, the valid offsets depends on wether the cell is on a border or not.
    For each 'edge', a list of offsets and its mDim are stored.*/
    
private:
    
    /// array of index offset to neighbors, for each edge type
    int  * regionsEdge;
    
    /// pointers to regionsEdge[], as a function of cell index
    int ** regions;
    
private:
    
    /// calculate the edge-characteristic from the size `s`, coordinate `c` and range `r`
    static size_t edge_signature(const size_t s, const int r, const int c)
    {
        if ( c < r )
            return (size_t)( r - c );
        else if ( (size_t)( c + r + 1 ) > s )
            return (size_t)( 2 * r + c + 1 ) - s;
        else
            return 0;
    }
    
    /// caculate the edge characteristic from the coordinates of a cell and the range vector
    size_t edgeFromCoordinates(const int coord[ORD], const size_t range[ORD]) const
    {
        size_t e = 0;
        for ( int d = ORD; d > 0; --d )
        {
            e *= 2 * range[d-1] + 1;
            e += edge_signature(mDim[d-1], range[d-1], coord[d-1]);
        }
        return e;
    }
    
    
    int * newRectangularGrid(size_t& cmx, const size_t range[ORD])
    {
        cmx = 1;
        for ( int d = 0; d < ORD; ++d )
            cmx *= ( 2 * range[d] + 1 );
        int * ccc = new int[ORD*cmx];
        
        size_t nb = 1;
        for ( int d = 0; d < ORD; ++d )
        {
            size_t h = 0;
            for ( ; h < nb && h < cmx; ++h )
                ccc[ORD*h+d] = 0;
            for ( size_t s = 1; s <= range[d]; ++s )
            {
                for ( size_t n = 0; n < nb; ++n )
                {
                    for ( int e = 0; e < d; ++e )
                        ccc[ORD*h+e] = ccc[ORD*n+e];
                    ccc[ORD*h+d] = (int)s;
                    ++h;
                    for ( int e = 0; e < d; ++e )
                        ccc[ORD*h+e] = ccc[ORD*n+e];
                    ccc[ORD*h+d] = -(int)s;
                    ++h;
                }
            }
            nb = h;
        }
        assert_true(nb==cmx);
        return ccc;
    }
    
    
    /// calculate cell index offsets between 'ori' and 'ori+shift'
    int calculateOffsets(int offsets[], int shift[], size_t cnt, int ori[], bool positive)
    {
        int nb = 0;
        int cc[ORD];
        int ori_indx = (int)pack(ori);
        for ( size_t ii = 0; ii < cnt; ++ii )
        {
            for ( int d = 0; d < ORD; ++d )
                cc[d] = ori[d] + shift[ORD*ii+d];
            int off = (int)pack(cc) - ori_indx;
            
            bool add = ( positive ? off >= 0 : true );
            if ( isPeriodic() )
            {
                //check that cell is not already included:
                for ( int n = 0; n < nb; ++n )
                    if ( offsets[n] == off )
                    {
                        add = false;
                        break;
                    }
            }
            else 
                add &= inside(cc);
            
            if ( add )
                offsets[nb++] = off;
        }
        return nb;
    }
    
   
    /// create regions in the offsets buffer
    /**
     Note: the range is taken specified in units of cells: 1 = 1 cell
     @todo: specify range in calculateRegion() as real distance!
     */
    void createRegions(int * ccc, const size_t regMax, const size_t range[ORD], bool positive)
    {
        size_t edgeMax = 0;
        for ( int d = ORD-1; d >= 0; --d )
        {
            edgeMax *= 2 * range[d] + 1;
            edgeMax += 2 * range[d];
        }
        ++edgeMax;
        
        //allocate and reset arrays:
        deleteRegions();
        
        regions     = new int*[mNbCells];
        regionsEdge = new int[edgeMax*(regMax+1)];
        for ( size_t e = 0; e < edgeMax*(regMax+1); ++e )
            regionsEdge[e] = 0;
        
        int ori[ORD];
        for ( size_t indx = 0; indx < mNbCells; ++indx )
        {
            setCoordinatesFromIndex(ori, indx);
            size_t e = edgeFromCoordinates(ori, range);
            assert_true( e < edgeMax );
            int * reg = regionsEdge + e * ( regMax + 1 );
            if ( reg[0] == 0 )
            {
                // calculate the region for this new edge-characteristic
                reg[0] = calculateOffsets(reg+1, ccc, regMax, ori, positive);
                //printf("edge %i has region of %i cells\n", e, reg[0]);
            }
            else if ( 0 )
            {
                // compare result for a different cell of the same edge-characteristic
                int * rig = new int[regMax+1];
                rig[0] = calculateOffsets(rig+1, ccc, regMax, ori, positive);
                if ( rig[0] != reg[0] )
                    ABORT_NOW("inconsistent region size");
                for ( int s = 1; s < rig[0]+1; ++s )
                    if ( rig[s] != reg[s] )
                        ABORT_NOW("inconsistent region offsets");
                delete[] rig;
            }
            regions[indx] = reg;
        }
    }
    
    /// accept within a certain diameter
    bool reject_disc(const int c[ORD], real radius)
    {
        real dsq = 0;
        for ( int d = 0; d < ORD; ++d ) 
            dsq += cWidth[d] * c[d] * cWidth[d] * c[d];
        return ( dsq > radius * radius );
    }
    
    /// accept within a certain diameter
    bool reject_square(const int c[ORD], real radius)
    {
        for ( int d = 0; d < ORD; ++d ) 
            if ( abs_real( cWidth[d] * c[d] ) > radius )
                return true;
        return false;
    }
    
    /// used to remove ORD coordinates in array `ccc`
    void remove_entry(int * ccc, const size_t s, size_t& cmx)
    {
        --cmx;
        for ( size_t x = ORD*s; x < ORD*cmx; ++x )
            ccc[x] = ccc[x+ORD];
    }
    
public:
    
    /// create regions which contains cells at a distance 'range' or less
    /**
     Note: the range is specified in real units.
     the region will cover an area of space that is approximately square.
     */
    void createSquareRegions(const real radius)
    {
        size_t range[ORD];
        for ( int d = 0; d < ORD; ++d )
            range[d] = std::ceil( radius / cWidth[d] );
        size_t cmx = 0;
        int * ccc = newRectangularGrid(cmx, range);
        
        for ( size_t s = cmx; s > 0 ; --s )
            if ( reject_square(ccc+ORD*(s-1), radius) )
                remove_entry(ccc, s-1, cmx);
        
        createRegions(ccc, cmx, range, false);
        delete[] ccc;
    }
    
    /// create regions which contains cells at a distance 'range' or less
    /**
     Note: the range is specified in real units.
     The region covers an area of space that is approximately circular.
     */
    void createRoundRegions(const real radius)
    {
        size_t range[ORD];
        for ( int d = 0; d < ORD; ++d )
        {
            assert_true( cWidth[d] > 0 );
            range[d] = std::ceil( radius / cWidth[d] );
        }
        size_t cmx = 0;
        int * ccc = newRectangularGrid(cmx, range);
       
        for ( size_t s = cmx; s > 0 ; --s )
            if ( reject_disc(ccc+ORD*(s-1), radius) )
                remove_entry(ccc, s-1, cmx);
        
        createRegions(ccc, cmx, range, false);
        delete[] ccc;
    }

    /// regions that only contain cells of greater index.
    /**
     This is suitable for pair-wise interaction of particles, since it can
     be used to go through the cells one by one such that at the end, all
     pairs of cells have been considered only once.

     Note: the radius is taken specified in units of cells: 1 = 1 cell
     */
    void createSideRegions(const int radius)
    {
        size_t range[ORD];
        for ( int d = 0; d < ORD; ++d )
            range[d] = radius;
        size_t cmx = 0;
        int * ccc = newRectangularGrid(cmx, range);
        createRegions(ccc, cmx, range, true);
        delete[] ccc;
    }
    
    /// true if createRegions() or createRoundRegions() was called
    bool hasRegions() const
    {
        return ( regions && regionsEdge );
    }
    
    /// set region array 'offsets' for given cell index
    /**
     A zero offset is always first in the list.
     //\returns the size of the list.
     
     \par Example:

         CELL * cell = & map.icell(indx);
         n_offset = map.getRegion(offset, indx);
         for ( int n = 1; n < n_offset; ++n )
         {
             Cell & neighbor = cell[offset[n]];
             ...
         }
     
     Note: createRegions() must be called first
    */
    int getRegion(int*& offsets, const size_t indx) const
    {
        assert_true( hasRegions() );
        offsets = regions[indx]+1;
        assert_true( offsets[0] == 0 );
        return regions[indx][0];
    }
    
    /// free memory claimed by the regions
    void deleteRegions()
    {
        delete[] regions;
        regions = nullptr;
        
        delete[] regionsEdge;
        regionsEdge = nullptr;
    }

#pragma mark -

    /// write total number of cells and number of subdivision in each dimension
    void printSummary(std::ostream& os, const char arg[])
    {
        char str[512], *ptr = str;
        char*const end = str+sizeof(str);
        ptr += snprintf(ptr, end-ptr, "%s of dim %i has %lu cells:", arg, ORD, mNbCells);
        for ( int d = 0; d < ORD; ++d )
            ptr += snprintf(ptr, end-ptr, " [ %9.3f %+9.3f ] / %lu = %4.3f", mInf[d], mSup[d], mDim[d], cWidth[d]);
        os << str << std::endl;
    }
};


#endif
