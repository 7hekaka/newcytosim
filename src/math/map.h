// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Francois Nedelec; Created 07/03/2015. nedelec@embl.de

#ifndef MAP_H
#define MAP_H

#include "assert_macro.h"
#include "exceptions.h"
#include <cstdio>
#include <cmath>
#include "real.h"
#include "modulo.h"

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
     8   9  10  11
     4   5   6   7
     0   1   2   3

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
    Map<ORD>& operator = (Map<ORD> const&);

    /// type for the edge signature of cells
    typedef unsigned char edge_type;
    
protected:
   
    /// Total number of cells in the map; size of cells[]
    size_t mNbCells;
    
    /// The number of cells in each dimension
    size_t mDim[ORD<4?4:ORD];
    
    /// Offset between two consecutive cells along each dimension
    size_t mStride[ORD];
    
    /// The position of the inferior (min) edge in each dimension
    real mInf[ORD];
    
    /// The position of the superior (max) edge in each dimension
    real mSup[ORD];

    /// The size of a cell: cWidth[d] = ( mSup[d] - inf[d] ) / mDim[d]
    real cWidth[ORD];
    
    /// cDelta[d] = 1.0 / cWidth[d]
    real cDelta[ORD];
    
    /// mStart[d] = mInf[d] / cWidth[d]
    real mStart[ORD];

    /// The volume of one cell
    real cVolume;
    
    /// true if Map has periodic boundary conditions
    bool mPeriodic[ORD];

protected:
    
    /// return closest integer to `c` in the segment [ 0, s-1 ]
    inline size_t image_(int d, long x) const
    {
        size_t S(mDim[d]);
        ///@todo use remainder() function for branchless code?
        while ( x <  0 ) x += S;
        size_t u(static_cast<size_t>(x));
        while ( u >= S ) u -= S;
        return u;
    }

    inline size_t clamp_(int d, long x) const
    {
        size_t S(mDim[d]-1);
        size_t u(static_cast<size_t>(std::max(0L, x)));
        return std::min(u, S);
        //return c <= 0 ? 0 : ( c >= s ? s-1 : c );
    }
    
    /// return closest integer to `c` in the segment [ 0, mDim[d]-1 ]
    inline size_t ind_(const int d, long x) const
    {
#if ENABLE_PERIODIC_BOUNDARIES
        if ( mPeriodic[d] )
            return image_(d, x);
        else
#endif
            return clamp_(d, x);
    }


    /// return f modulo s in [ 0, s-1 ]
    static inline size_t imagef_periodic(size_t s, real f)
    {
        while ( f <  0 ) f += (real)s;
        size_t u((size_t)f);
        while ( u >= s ) u -= s;
        return u;
    }

    static inline size_t imagef_clamped(size_t s, real f)
    {
        if ( f > 0 )
        {
            return std::min(size_t(f), s);
        }
        return 0;
    }
    
    
    /// returns  ( f - mInf[d] ) / cWidth[d]
    inline long map_(const int d, real f) const
    {
        return static_cast<long>( f * cDelta[d] - mStart[d] );
    }
    
    /// returns  0.5 + ( f - mInf[d] ) / cWidth[d]
    inline real mapC(const int d, real f) const
    {
        return f * cDelta[d] - ( mStart[d] - 0.5 );
    }

    /// return closest integer to `c` in the segment [ 0, mDim[d]-1 ]
    inline size_t imagef(const int d, real f) const
    {
        real x = f * cDelta[d] - mStart[d];
#if ENABLE_PERIODIC_BOUNDARIES
        if ( mPeriodic[d] )
            return imagef_periodic(mDim[d], x);
        else
#endif
            return imagef_clamped(mDim[d]-1, x);
    }

    /// return closest integer to `c` in the segment [ 0, mDim[d]-1 ]
    inline size_t imagef(const int d, real f, real offset) const
    {
        real x = f * cDelta[d] - ( mStart[d] - offset );
#if ENABLE_PERIODIC_BOUNDARIES
        if ( mPeriodic[d] )
            return imagef_periodic(mDim[d], x);
        else
#endif
            return imagef_clamped(mDim[d]-1, x);
    }

//--------------------------------------------------------------------------
#pragma mark -
public:
    
    /// constructor
    Map() : mDim{0}, mInf{0}, mSup{0}, cWidth{0}, cDelta{0}, mStart{0}, mPeriodic{false}
    {
        mNbCells = 0;
        eRegion = nullptr;
        region_ = nullptr;
        cVolume = 0;
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
        bool reshaped = false;
        
        for ( int d = 0; d < ORD; ++d )
        {
            if ( cells[d] <= 0 )
                throw InvalidParameter("Cannot build grid as nbcells[] is <= 0");
            
            if ( infs[d] >= sups[d] )
                throw InvalidParameter("Cannot build grid as sup[] <= inf[]");
            reshaped |= ( mDim[d] != cells[d] );
            
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
        if ( reshaped )
            deleteRegions();
    }
    
    ///true if setDimensions() was called
    bool hasDimensions() const
    {
        return mNbCells > 0;
    }
    
    /// true if dimension `d` has periodic boundary conditions
    bool isPeriodic(int d) const
    {
#if ENABLE_PERIODIC_BOUNDARIES
        if ( d < ORD )
            return mPeriodic[d];
#endif
        return false;
    }
    
    /// change boundary conditions
    void setPeriodic(int d, bool p)
    {
#if ENABLE_PERIODIC_BOUNDARIES
        if ( d < ORD )
            mPeriodic[d] = p;
#else
        if ( p )
            throw InvalidParameter("grid.h was compiled with ENABLE_PERIODIC_BOUNDARIES=0");
#endif
    }
    
    /// true if boundary conditions are periodic
    bool isPeriodic() const
    {
#if ENABLE_PERIODIC_BOUNDARIES
        for ( int d = 0; d < ORD; ++d )
            if ( mPeriodic[d] )
                return true;
#endif
        return false;
    }

    //--------------------------------------------------------------------------
#pragma mark -

    /// total number of cells in the map
    size_t nbCells()      const { return mNbCells; }

    /// number of cells in dimensionality `d`
    size_t breadth(int d) const { return mDim[d]; }
    
    /// offset to the next cell in the direction `d`
    size_t stride(int d)  const { return mStride[d]; }

    /// position of the inferior (left/bottom/front) edge
    real inf(int d)   const { return mInf[d]; }
    
    /// position of the superior (right/top/back) edge
    real sup(int d)   const { return mSup[d]; }
    
    /// the widths of a cell
    real cellWidth(int d) const { return cWidth[d]; }
    
    /// inverse of the widths of a cell
    real delta(int d)     const { return cDelta[d]; }
    
    /// access to data vectors
    const real * inf()       const { return mInf; }
    const real * sup()       const { return mSup; }
    const real * cellWidth() const { return cWidth; }
    const real * delta()     const { return cDelta; }

    /// the volume of a cell
    real cellVolume() const { return cVolume; }

    /// position in dimension `d`, of the cell of index `c`
    real position(int d, real c) const { return mInf[d] + c * cWidth[d]; }
    
    /// index in dimension `d` corresponding to position `w`
    long index(int d, real w) const { return map_(d, w); }

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
            coord[d] = ind_(d, coord[d]);
    }
    
    /// conversion from index to coordinates
    void setCoordinatesFromIndex(int coord[ORD], size_t indx) const
    {
        for ( int d = 0; d < ORD; ++d )
        {
            coord[d] = indx % mDim[d];
            indx /= mDim[d];
        }
    }
    
    /// conversion from Position to coordinates (offset should be in [0,1])
    void setCoordinatesFromPosition(int coord[ORD], const real w[ORD], const real offset=0) const
    {
        for ( int d = 0; d < ORD; ++d )
            coord[d] = imagef(d, w[d], offset);
    }
    
    void setPositionFromIndex(real res[ORD], size_t indx, real offset) const
    {
        for ( int d = 0; d < ORD; ++d )
        {
            res[d] = mInf[d] + cWidth[d] * ( offset + indx % mDim[d] );
            indx /= mDim[d];
        }
    }

    /// conversion from Index to Position (offset should be in [0,1])
    template < typename REAL >
    void setPositionFromIndex(REAL res[ORD], size_t indx, real offset) const
    {
        for ( int d = 0; d < ORD; ++d )
        {
            res[d] = mInf[d] + cWidth[d] * ( offset + indx % mDim[d] );
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
        size_t inx = ind_(ORD-1, coord[ORD-1]);
        
        for ( int d = ORD-2; d >= 0; --d )
            inx = mDim[d] * inx + ind_(d, coord[d]);
        
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

        s[dim] = ind_(dim, s[dim]+1);

        c = s[ORD-1];
        for ( int d = ORD-2; d >= 0; --d )
            c = mDim[d] * c + s[d];
        return c;
    }
    
    /// convert coordinate to array index, if ORD==1
    size_t pack1D(const int x) const
    {
        return ind_(0, x);
    }
    
    /// convert coordinate to array index, if ORD==2
    size_t pack2D(const int x, const int y) const
    {
        return ind_(1, y) * mDim[0] + ind_(0, x);
    }
    
    /// convert coordinate to array index, if ORD==3
    size_t pack3D(const int x, const int y, const int z) const
    {
        return ( ind_(2, z) * mDim[1] + ind_(1, y) ) * mDim[0] + ind_(0, x);
    }
    
    /// convert coordinate to array index, if ORD==1
    size_t pack1D_clamped(const int x) const
    {
        return clamp_(0, x);
    }
    
    /// convert coordinate to array index, if ORD==2
    size_t pack2D_clamped(const int x, const int y) const
    {
        return clamp_(1, y) * mDim[0] + clamp_(0, x);
    }
    
    /// convert coordinate to array index, if ORD==3
    size_t pack3D_clamped(const int x, const int y, const int z) const
    {
        return ( clamp_(2, z) * mDim[1] + clamp_(1, y) ) * mDim[0] + clamp_(0, x);
    }

    
    /// return index of cell corresponding to position (x), if ORD==1
    size_t index1D(const real x) const
    {
        return ind_(0, map_(0, x));
    }
    
    /// return index of cell corresponding to position (x, y), if ORD==2
    size_t index2D(const real x, const real y) const
    {
        size_t X = ind_(0, map_(0, x));
        size_t Y = ind_(1, map_(1, y));
        return  Y * mDim[0] + X;
    }
    
    /// return index of cell corresponding to position (x, y, z), if ORD==3
    size_t index3D(const real x, const real y, const real z) const
    {
        size_t X = ind_(0, map_(0, x));
        size_t Y = ind_(1, map_(1, y));
        size_t Z = ind_(2, map_(2, z));
        return ( Z * mDim[1] + Y ) * mDim[0] + X;
    }
    
    /// return index of cell corresponding to position (x, y), if ORD==2
    size_t direct_index2D(const real x, const real y) const
    {
        size_t X = static_cast<size_t>(map_(0, x)); assert_true( X < mDim[0] );
        // with semi-periodic conditions, Y may not be inside, and clamping is necessary:
        size_t Y = static_cast<size_t>(clamp_(1, map_(1, y))); assert_true( Y < mDim[1] );
        return Y * mDim[0] + X;
    }

    /// return index of cell corresponding to position (x, y, z), if ORD==3
    size_t direct_index3D(const real x, const real y, const real z) const
    {
        size_t X = static_cast<size_t>(map_(0, x)); assert_true( X < mDim[0] );
        // with semi-periodic conditions, Y may not be inside, and clamping is necessary:
        size_t Y = static_cast<size_t>(clamp_(1, map_(1, y))); assert_true( Y < mDim[1] );
        size_t Z = static_cast<size_t>(clamp_(2, map_(2, z))); assert_true( Z < mDim[2] );
        return ( Z * mDim[1] + Y ) * mDim[0] + X;
    }

    //--------------------------------------------------------------------------

#pragma mark - Regions

    /** For any cell, we can find the adjacent cells by adding 'index offsets'
    However, the valid offsets depends on wether the cell is on a border or not.
    For each 'edge', a list of offsets and its mDim are stored.*/
    
private:
    
    /// array of index offset to neighbors, for each edge type
    int * eRegion;
    
    /// index into eRegion[], as a function of cell index
    edge_type * region_;
    
private:
    
    /// calculate the edge-characteristic from the size `s`, range `r` and coordinate `c`
    static int edge_signature(const int s, const int r, const int c)
    {
        /* the characteristic is positive, and given that 0 <= c <= s-1,
         if will be in [ 0, 2*r+1 ], with 0 far away from any edge.
         For example with r = 2, we get characteristics as follows:
         [ 4, 3, 0, 0 ... 0, 1, 2 ]
         */
        if ( c < r )
            return ( 2 * r - c );
        else
            return std::max( 0, r + c + 1 - s );
    }
    
    /// caculate the edge characteristic from the coordinates of a cell and the range vector
    edge_type edgeFromCoordinates(const int coord[ORD], const size_t range[ORD], bool symmetric) const
    {
        int e = 0;
        if ( symmetric )
        {
            for ( int d = ORD-1; d >= 0; --d )
            {
                e *= 2 * range[d] + 1;
                e += edge_signature(mDim[d], range[d], coord[d]);
            }
        }
        else
        {
            for ( int d = ORD-1; d >= 0; --d )
            {
                e *= range[d] + 1;
                e += std::max( 0, coord[d] + (int)range[d] + 1 - (int)mDim[d] );
            }
        }
        assert_true(static_cast<edge_type>(e) == e);
        return static_cast<edge_type>(e);
    }
    
    /// return array of dimensionality ORD, containing indices with reference to the center cell
    static size_t initRectangularGrid(int * ccc, size_t ccc_size, const size_t range[ORD], bool symmetric)
    {
        size_t res = 1;
        for ( unsigned d = 0; d < ORD; ++d )
            ccc[d] = 0;
        for ( unsigned d = 0; d < ORD; ++d )
        {
            size_t h = res;
            for ( int s = 1; s <= range[d]; ++s )
            {
                for ( size_t n = 0; n < res; ++n )
                {
                    for ( unsigned e = 0; e < d; ++e )
                        ccc[ORD*h+e] = ccc[ORD*n+e];
                    ccc[ORD*h+d] = (int)s;
                    ++h;
                    if ( symmetric )
                    {
                        for ( unsigned e = 0; e < d; ++e )
                            ccc[ORD*h+e] = ccc[ORD*n+e];
                        ccc[ORD*h+d] = -(int)s;
                        ++h;
                    }
                }
            }
            res = h;
        }
        //printf("initRectangularGrid %i has %i neighbors\n", ccc_size, res);
        return res;
    }
    
    
    /// calculate cell index offsets between 'ori' and 'ori+shift'
    int calculateOffsets(int offsets[], int shift[], size_t cnt, int ori[], bool symmetric)
    {
        int res = 0;
        int cc[ORD];
        int ori_indx = (int)pack(ori);
        for ( size_t ii = 0; ii < cnt; ++ii )
        {
            for ( unsigned d = 0; d < ORD; ++d )
                cc[d] = ori[d] + shift[ORD*ii+d];
            int off = (int)pack(cc) - ori_indx;
            
            if ( symmetric || off >= 0 )
            {
                bool add = true;
                if ( isPeriodic() )
                {
                    //check that cell is not already included:
                    for ( int n = 0; n < res; ++n )
                        if ( offsets[n] == off )
                        {
                            add = false;
                            break;
                        }
                }
                else
                    add &= inside(cc);
                
                if ( add )
                    offsets[res++] = off;
            }
        }
        return res;
    }
    
   
    /// create regions in the offsets buffer
    /**
     Note: the range is taken specified in units of cells: 1 = 1 cell
     @todo: specify range in calculateRegion() as real distance!
     */
    void createRegions(int * ccc, const edge_type regMax, const size_t range[ORD], bool symmetric)
    {
        //allocate and reset arrays:
        deleteRegions();
        
        region_ = new edge_type[mNbCells]{0};
        eRegion = new int[regMax*(regMax+1)]{0};
        
        int ori[ORD]{0};
        for ( size_t i = 0; i < mNbCells; ++i )
        {
            setCoordinatesFromIndex(ori, i);
            edge_type e = edgeFromCoordinates(ori, range, symmetric);
            assert_true( e < regMax );
            edge_type off = e * regMax + e;
            int * reg = eRegion + off;
            if ( reg[0] == 0 )
            {
                // calculate the region for this new edge-characteristic
                reg[0] = calculateOffsets(reg+1, ccc, regMax, ori, symmetric);
                //printf("edge type %i has %i neighbors\n", e, reg[0]);
            }
            region_[i] = off;
#if ( 0 )
            // compare result for a different cell of the same edge-characteristic
            int * rig = new int[regMax+1];
            rig[0] = calculateOffsets(rig+1, ccc, regMax, ori, symmetric);
            if ( rig[0] != reg[0] )
                ABORT_NOW("inconsistent region size");
            for ( int s = 1; s < rig[0]+1; ++s )
                if ( rig[s] != reg[s] )
                    ABORT_NOW("inconsistent region offsets");
            delete[] rig;
#endif
        }
    }
    
    /// accept within a certain diameter
    bool reject_disc(const int c[ORD], real radius)
    {
        real dsq = 0;
        for ( int d = 0; d < ORD; ++d ) 
            dsq += square(cWidth[d] * c[d]);
        return ( dsq > square(radius) );
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
        size_t cmx = 1;
        size_t range[ORD];
        for ( int d = 0; d < ORD; ++d )
        {
            assert_true( cWidth[d] > REAL_EPSILON );
            size_t R = std::ceil( radius / cWidth[d] );
            cmx *= ( 2 * R + 1 );
            range[d] = R;
        }
        if ( cmx != (edge_type)cmx )
            throw InvalidParameter("Region size exceeds edge_type capacity");

        int * ccc = new int[ORD*cmx]{0};
        initRectangularGrid(ccc, cmx, range, true);
        
        for ( size_t s = cmx; s-- > 0 ; )
            if ( reject_square(ccc+ORD*s, radius) )
                remove_entry(ccc, s, cmx);
        
        createRegions(ccc, cmx, range, true);
        delete[] ccc;
    }
    
    /// create regions which contains cells at a distance 'radius' or less
    /**
     Note: the range is specified in real units.
     The region covers an area of space that is approximately circular.
     */
    void createRoundRegions(const real radius)
    {
        size_t cmx = 1;
        size_t range[ORD];
        for ( int d = 0; d < ORD; ++d )
        {
            assert_true( cWidth[d] > REAL_EPSILON );
            size_t R = std::ceil( radius / cWidth[d] );
            cmx *= ( 2 * R + 1 );
            range[d] = R;
        }
        if ( cmx != (edge_type)cmx )
            throw InvalidParameter("Region size exceeds edge_type capacity");

        int * ccc = new int[ORD*cmx]{0};
        initRectangularGrid(ccc, cmx, range, true);

        for ( size_t s = cmx; s-- > 0 ; )
            if ( reject_disc(ccc+ORD*s, radius) )
                remove_entry(ccc, s, cmx);
        
        createRegions(ccc, cmx, range, true);
        delete[] ccc;
    }

    /// regions that only contain cells of greater index.
    /**
     This is suitable for pair-wise interaction of particles, since it can
     be used to go through the cells one by one such that at the end, all
     pairs of cells have been considered only once.

     Note: the radius is taken specified in units of cells: 1 = 1 cell
     */
    void createSideRegions(const unsigned num_cells_radius)
    {
        size_t cmx = 1;
        size_t range[ORD];
        for ( int d = 0; d < ORD; ++d )
        {
            cmx *= ( num_cells_radius + 1 );
            range[d] = num_cells_radius;
        }
        if ( cmx != (edge_type)cmx )
            throw InvalidParameter("Region size exceeds edge_type capacity");
        int * ccc = new int[ORD*cmx]{0};
        initRectangularGrid(ccc, cmx, range, false);
        createRegions(ccc, cmx, range, false);
        delete[] ccc;
    }
    
    
    /// true if createRegions() or createRoundRegions() was called
    bool hasRegions() const
    {
        return ( region_ && eRegion );
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
    int getRegion(int const*& offsets, const size_t indx) const
    {
        assert_true( hasRegions() );
        int * R = eRegion + region_[indx];
        offsets = R + 1;
        assert_true( offsets[0] == 0 );
        return R[0];
    }
    
    /// free memory claimed by the regions
    void deleteRegions()
    {
        delete[] region_;
        region_ = nullptr;
        
        delete[] eRegion;
        eRegion = nullptr;
    }

#pragma mark -

    /// write total number of cells and number of subdivision in each dimension
    void printSummary(std::ostream& os, const char arg[])
    {
        os << arg << " of dim " << ORD << " has " << mNbCells << " cells: ";
        for ( int d = 0; d < ORD; ++d )
        {
            char o = '[', c = ']';
            if ( mPeriodic[d] ) { o = ']'; c = '['; }
            os << " " << o << mInf[d] << ", " << mSup[d];
            os << c << "/" << mDim[d] << " = " << cWidth[d];
        }
        os << std::endl;
    }
};


#endif
