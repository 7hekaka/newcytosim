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
    
For a 4x4 2D grid, the indices are distributed like this:

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
In the example above:
    - for index = 0 it would return { 1 4 5 }
    - for index = 5 it would return { -1 1 -5 -4 -3 3 4 5 }
.
Note that, as of 16.04.2025, ** the zero offset for self is not included anymore. **

The indices of the neighboring cells are calculated by adding offsets[n] to 'index':
Example:

        CELL * cell = & map.icell(indx);
        int cnt = map.getRegion(region, indx);
        for ( int n = 0; n < cnt; ++n )
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
    index_t mNbCells;
    
    /// The number of cells in each dimension
    index_t mDim[ORD<4?4:ORD];
    
    /// Offset between two consecutive cells along each dimension
    index_t mStride[ORD];
    
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
    inline index_t image_(int d, int x) const
    {
        index_t S(mDim[d]);
        ///@todo use remainder() function for branchless code?
        while ( x <  0 ) x += S;
        index_t u(static_cast<index_t>(x));
        while ( u >= S ) u -= S;
        return u;
    }

    inline index_t clamp_(int d, int x) const
    {
        index_t S(mDim[d]-1);
        index_t u(static_cast<index_t>(std::max(0, x)));
        return std::min(u, S);
        //return c <= 0 ? 0 : ( c >= s ? s-1 : c );
    }
    
    /// return closest integer to `c` in the segment [ 0, mDim[d]-1 ]
    inline index_t ind_(const int d, int x) const
    {
#if ENABLE_PERIODIC_BOUNDARIES
        if ( mPeriodic[d] )
            return image_(d, x);
        else
#endif
            return clamp_(d, x);
    }


    /// return f modulo s in [ 0, s-1 ]
    static inline index_t imagef_periodic(index_t s, real f)
    {
        while ( f <  0 ) f += (real)s;
        index_t u((index_t)f);
        while ( u >= s ) u -= s;
        return u;
    }

    static inline index_t imagef_clamped(index_t s, real f)
    {
        if ( f > 0 )
        {
            return std::min(index_t(f), s);
        }
        return 0;
    }
    
    
    /// returns  ( f - mInf[d] ) / cWidth[d]
    inline int map_(const int d, real f) const
    {
        return static_cast<int>( f * cDelta[d] - mStart[d] );
    }
    
    /// returns  0.5 + ( f - mInf[d] ) / cWidth[d]
    inline real mapC(const int d, real f) const
    {
        return f * cDelta[d] - ( mStart[d] - 0.5 );
    }

    /// return closest integer to `c` in the segment [ 0, mDim[d]-1 ]
    inline index_t imagef(const int d, real f) const
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
    inline index_t imagef(const int d, real f, real offset) const
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
        chunk_ = 0;
        region_ = nullptr;
        border_ = nullptr;
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
    void setDimensions(const real infs[ORD], real sups[ORD], const index_t cells[ORD])
    {
        cVolume = 1;
        index_t cnt = 1;
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
    index_t nbCells()      const { return mNbCells; }

    /// number of cells in dimensionality `d`
    index_t breadth(int d) const { return mDim[d]; }
    
    /// offset to the next cell in the direction `d`
    index_t stride(int d)  const { return mStride[d]; }

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
    int index(int d, real w) const { return map_(d, w); }

    /// half the diagonal length of the unit cell
    real cellRadius() const
    {
        real res = cWidth[0] * cWidth[0];
        for ( int d = 1; d < ORD; ++d )
            res += cWidth[d] * cWidth[d];
        return 0.5 * std::sqrt(res);
    }
    
    /// the smallest cell width, along dimensions that have more than `min_size` cells
    real minimumWidth(index_t min_size) const
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
            if ( coord[d] < 0 || (index_t)coord[d] >= mDim[d] )
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
    void setCoordinatesFromIndex(int coord[ORD], index_t indx) const
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
    
    void setPositionFromIndex(real res[ORD], index_t indx, real offset) const
    {
        for ( int d = 0; d < ORD; ++d )
        {
            res[d] = mInf[d] + cWidth[d] * ( offset + indx % mDim[d] );
            indx /= mDim[d];
        }
    }

    /// conversion from Index to Position (offset should be in [0,1])
    template < typename REAL >
    void setPositionFromIndex(REAL res[ORD], index_t indx, real offset) const
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
    index_t pack(const int coord[ORD]) const
    {
        index_t inx = ind_(ORD-1, coord[ORD-1]);
        
        for ( int d = ORD-2; d >= 0; --d )
            inx = mDim[d] * inx + ind_(d, coord[d]);
        
        return inx;
    }
    
    
    /// returns the index of the cell whose center is closest to the point w[]
    index_t index(const real w[ORD]) const
    {
        index_t inx = imagef(ORD-1, w[ORD-1]);
        
        for ( int d = ORD-2; d >= 0; --d )
            inx = mDim[d] * inx + imagef(d, w[d]);
        
        return inx;
    }

    
    /// returns the index of the cell whose center is closest to the point w[]
    index_t index(const real w[ORD], const real offset) const
    {
        index_t inx = imagef(ORD-1, w[ORD-1], offset);
        
        for ( int d = ORD-2; d >= 0; --d )
            inx = mDim[d] * inx + imagef(d, w[d], offset);
        
        return inx;
    }

    
    /// return cell that is next to `c` in the direction `dim`
    index_t next(index_t c, int dim) const
    {
        index_t s[ORD];
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
    index_t pack1D(const int x) const
    {
        return ind_(0, x);
    }
    
    /// convert coordinate to array index, if ORD==2
    index_t pack2D(const int x, const int y) const
    {
        return ind_(1, y) * mDim[0] + ind_(0, x);
    }
    
    /// convert coordinate to array index, if ORD==3
    index_t pack3D(const int x, const int y, const int z) const
    {
        return ( ind_(2, z) * mDim[1] + ind_(1, y) ) * mDim[0] + ind_(0, x);
    }
    
    /// convert coordinate to array index, if ORD==1
    index_t pack1D_clamped(const int x) const
    {
        return clamp_(0, x);
    }
    
    /// convert coordinate to array index, if ORD==2
    index_t pack2D_clamped(const int x, const int y) const
    {
        return clamp_(1, y) * mDim[0] + clamp_(0, x);
    }
    
    /// convert coordinate to array index, if ORD==3
    index_t pack3D_clamped(const int x, const int y, const int z) const
    {
        return ( clamp_(2, z) * mDim[1] + clamp_(1, y) ) * mDim[0] + clamp_(0, x);
    }

    
    /// return index of cell corresponding to position (x), if ORD==1
    index_t index1D(const real x) const
    {
        return ind_(0, map_(0, x));
    }
    
    /// return index of cell corresponding to position (x, y), if ORD==2
    index_t index2D(const real x, const real y) const
    {
        index_t X = ind_(0, map_(0, x));
        index_t Y = ind_(1, map_(1, y));
        return  Y * mDim[0] + X;
    }
    
    /// return index of cell corresponding to position (x, y, z), if ORD==3
    index_t index3D(const real x, const real y, const real z) const
    {
        index_t X = ind_(0, map_(0, x));
        index_t Y = ind_(1, map_(1, y));
        index_t Z = ind_(2, map_(2, z));
        return ( Z * mDim[1] + Y ) * mDim[0] + X;
    }
    
    /// return index of cell corresponding to position (x, y), if ORD==2
    index_t direct_index2D(const real x, const real y) const
    {
        index_t X = static_cast<index_t>(map_(0, x)); assert_true( X < mDim[0] );
        // with semi-periodic conditions, Y may not be inside, and clamping is necessary:
        index_t Y = static_cast<index_t>(clamp_(1, map_(1, y))); assert_true( Y < mDim[1] );
        return Y * mDim[0] + X;
    }

    /// return index of cell corresponding to position (x, y, z), if ORD==3
    index_t direct_index3D(const real x, const real y, const real z) const
    {
        index_t X = static_cast<index_t>(map_(0, x)); assert_true( X < mDim[0] );
        // with semi-periodic conditions, Y may not be inside, and clamping is necessary:
        index_t Y = static_cast<index_t>(clamp_(1, map_(1, y))); assert_true( Y < mDim[1] );
        index_t Z = static_cast<index_t>(clamp_(2, map_(2, z))); assert_true( Z < mDim[2] );
        return ( Z * mDim[1] + Y ) * mDim[0] + X;
    }

    //--------------------------------------------------------------------------

#pragma mark - Regions

    /** For any cell, we can find the adjacent cells by adding 'index offsets'
    However, the valid offsets depends on wether the cell is on a border or not.
    For each 'edge', a list of offsets and its mDim are stored.*/
    
private:
    
    /// array of index offset to neighbors, for each edge type
    int * region_;
    
    /// the edge type, as a function of cell index
    edge_type * border_;
    
    /// size of a chunk of data in region_
    size_t chunk_;
    
private:
    
    /// calculate the edge-characteristic from the size `s`, range `r` and coordinate `c`
    static int edge_signature(const int s, const int r, const int c)
    {
        /* the characteristic is positive, and given that 0 <= c <= s-1,
         if will be in [ 0, 2*r+1 ], with 0 far away from any edge.
         For example with r = 2, we get characteristics as follows:
         [ 3, 1, 0, 0 ... 0, 2, 4 ]
         */
        if ( c < r ) return 2 * ( r - c ) - 1;
        return std::max(0, 2 * ( r + c - s + 1 ));
    }
    
    static edge_type numberEdgeTypes(const index_t range[ORD])
    {
        int r = 0;
        for ( int d = ORD-1; d >= 0; --d )
        {
            r *= 2 * range[d] + 1;
            r += 2 * range[d];
        }
        ++r;
        assert_true(static_cast<edge_type>(r) == r);
        return static_cast<edge_type>(r);
    }

    /// caculate the edge characteristic from the coordinates of a cell and the range vector
    edge_type edgeFromCoordinates(const int coord[ORD], const index_t range[ORD]) const
    {
        int r = 0;
        for ( int d = ORD-1; d >= 0; --d )
        {
            int e = edge_signature(mDim[d], range[d], coord[d]);
            assert_true( e < 2 * range[d] + 1 );
            r = ( 2 * range[d] + 1 ) * r + e;
        }
        assert_true(static_cast<edge_type>(r) == r);
        return static_cast<edge_type>(r);
    }
    
    /// return array of dimensionality ORD, containing indices with reference to the center cell
    static index_t initRectangularGrid(int * ccc, index_t sup, const index_t range[ORD])
    {
        index_t res = 1;
        for ( unsigned d = 0; d < ORD; ++d )
            ccc[d] = 0;
        for ( unsigned d = 0; d < ORD; ++d )
        {
            index_t row = res;
            for ( int s = 1; s <= range[d]; ++s )
            {
                for ( index_t n = 0; n < row; ++n )
                {
                    for ( unsigned e = 0; e < d; ++e )
                        ccc[ORD*res+e] = ccc[ORD*n+e];
                    ccc[ORD*res+d] = -s;
                    ++res;
                    for ( unsigned e = 0; e < d; ++e )
                        ccc[ORD*res+e] = ccc[ORD*n+e];
                    ccc[ORD*res+d] = s;
                    ++res;
                }
            }
        }
        assert_true(res <= sup);
        //printf("initRectangularGrid %i found %i neighbors\n", sup, res);
        return res;
    }
    
    
    /// return array of dimensionality ORD, containing indices with reference to the center cell
    static index_t keepPositiveOffsets(int * ccc, index_t sup)
    {
        index_t cnt = 1;
        for ( int i = 1; i < sup; ++i )
        {
            int * C = ccc + ORD * i;
            bool keep = false;
            for ( int d = ORD-1; d >= 0; --d )
            {
                if ( C[d] > 0 ) keep = true;
                else if ( C[d] < 0 ) break;
            }
            if ( keep )
            {
                for( int d = 0; d < ORD; ++d )
                    ccc[ORD*cnt+d] = C[d];
                ++cnt;
            }
        }
        //printf("keepPositiveOffsets kept %i offset\n", cnt);
        return cnt;
    }

    /// print grid
    void printGrid(std::ostream& os, int * ccc, index_t sup)
    {
        for ( int i = 0; i < sup; ++i )
        {
            os << "\n" << std::noshowpos << i << " :  ";
            for ( int d = 0; d < ORD; ++d )
                os << "   " << std::showpos << ccc[ORD*i+d];
        }
        os << std::noshowpos << "\n";
    }
    
    /// calculate cell index offsets between 'ori' and 'ori+shift'
    void calculateOffsets(int offsets[], const int shift[], index_t cnt, const int ori[], int ori_indx)
    {
        int res = 0;
        int cc[ORD] {0};
        for ( index_t i = 0; i < cnt; ++i )
        {
            for ( unsigned d = 0; d < ORD; ++d )
                cc[d] = ori[d] + shift[ORD*i+d];
            int off = (int)pack(cc) - ori_indx;
            
            if ( isPeriodic() )
            {
                //check that cell was not already included:
                bool add = true;
                for ( int n = 0; n < res; ++n )
                {
                    if ( offsets[n] == off )
                    {
                        add = false;
                        break;
                    }
                }
                if ( add )
                    offsets[res++] = off;
            }
            else if ( inside(cc) )
                offsets[res++] = off;
        }
        assert_true(res <= cnt);
        assert_true(offsets[0] == 0);
        offsets[0] = res - 1;
    }
    
   
    /// create regions in the offsets buffer
    /**
     Note: the range is taken specified in units of cells: 1 = 1 cell
     @todo: specify range in calculateRegion() as real distance!
     */
    void createRegions(int * ccc, const edge_type regSize, const index_t range[ORD])
    {
        deleteRegions();
        
        edge_type edgeMax = numberEdgeTypes(range);
        //printf("edgeMax %i regSize %i\n", edgeMax, regSize);

        chunk_ = regSize;
        border_ = new edge_type[mNbCells]{0};
        region_ = new int[regSize*edgeMax]{0};
        
        int ori[ORD]{0};
        for ( index_t i = 0; i < mNbCells; ++i )
        {
            setCoordinatesFromIndex(ori, i);
            edge_type e = edgeFromCoordinates(ori, range);
            assert_true( e < edgeMax );
            border_[i] = e;
            int * reg = region_ + chunk_ * (size_t)e;
            if ( reg[0] <= 0 )
            {
                // calculate the region for this new edge-characteristic
                calculateOffsets(reg, ccc, regSize, ori, i);
                //printf("cell %i with edge type %i has %i neighbors\n", i, e, reg[0]);
            }
#if ( 0 )
            else
            {
                printf("checking cell %i with edge type %i\n", i, e);
                // compare result for a different cell with the same edge-characteristic
                int * rig = new int[regSize];
                calculateOffsets(rig, ccc, regSize, ori, i);
                bool different = false;
                for ( int s = 0; s < regSize; ++s )
                    if ( rig[s] != reg[s] ) different = true;
                if ( different )
                {
                    for ( int s = 0; s < regSize; ++s )
                        printf(" %2i", reg[s]);
                    printf("\n");
                    for ( int s = 0; s < regSize; ++s )
                        printf(" %2i", rig[s]);
                    printf("\n");
                    ABORT_NOW("inconsistent regions");
                }
                delete[] rig;
            }
#endif
        }
    }
    
    /// accept within a certain diameter
    bool in_disc(const int c[ORD], real radius)
    {
        real dsq = 0;
        for ( int d = 0; d < ORD; ++d ) 
            dsq += square(cWidth[d] * c[d]);
        return ( dsq <= square(radius) );
    }
    
    /// accept within a certain diameter
    bool in_square(const int c[ORD], real radius)
    {
        for ( int d = 0; d < ORD; ++d ) 
            if ( abs_real( cWidth[d] * c[d] ) > radius )
                return false;
        return true;
    }
    
    static void checkEdgeTypeRange(size_t arg)
    {
        if ( arg != (edge_type)arg )
        {
            // if this error occurs, recompile using a larger type for edge_type
            fprintf(stderr, "FATAL ERROR: edge_type cannot hold %lu values\n", arg);
            throw InvalidParameter("Region size exceeds edge_type capacity");
        }
    }
    
public:
    
    /// create regions which contains cells at a distance 'range' or less
    /**
     Note: the range is specified in real units.
     the region will cover an area of space that is approximately square.
     */
    void createSquareRegions(const real radius)
    {
        index_t cmx = 1;
        index_t range[ORD];
        for ( int d = 0; d < ORD; ++d )
        {
            assert_true( cWidth[d] > REAL_EPSILON );
            index_t R = std::ceil( radius / cWidth[d] );
            cmx *= ( 2 * R + 1 );
            range[d] = R;
        }
        checkEdgeTypeRange(cmx);
        int * ccc = new int[ORD*cmx]{0};
        initRectangularGrid(ccc, cmx, range);
        
        int cnt = 0;
        for ( index_t s = 0; s < cmx; ++s )
        {
            if ( in_square(ccc+ORD*s, radius) )
            {
                for ( int d = 0; d < ORD; ++d )
                    ccc[ORD*cnt+d] = ccc[ORD*s+d];
                ++cnt;
            }
        }
        createRegions(ccc, cnt, range);
        delete[] ccc;
    }
    
    /// create regions which contains cells at a distance 'radius' or less
    /**
     Note: the range is specified in real units.
     The region covers an area of space that is approximately circular.
     */
    void createRoundRegions(const real radius)
    {
        index_t cmx = 1;
        index_t range[ORD];
        for ( int d = 0; d < ORD; ++d )
        {
            assert_true( cWidth[d] > REAL_EPSILON );
            index_t R = std::ceil( radius / cWidth[d] );
            cmx *= ( 2 * R + 1 );
            range[d] = R;
        }
        checkEdgeTypeRange(cmx);
        int * ccc = new int[ORD*cmx]{0};
        initRectangularGrid(ccc, cmx, range);

        int cnt = 0;
        for ( index_t s = 0; s < cmx; ++s )
        {
            if ( in_disc(ccc+ORD*s, radius) )
            {
                for ( int d = 0; d < ORD; ++d )
                    ccc[ORD*cnt+d] = ccc[ORD*s+d];
                ++cnt;
            }
        }
        createRegions(ccc, cnt, range);
        delete[] ccc;
    }

    /// regions that only contain cells of greater index.
    /**
     This is suitable for pair-wise interaction of particles, since it can
     be used to go through the cells one by one such that at the end, all
     pairs of cells have been considered only once.

     Note: the radius is taken specified in units of cells: 1 = 1 cell
     */
    void createSideRegions(const index_t num_cells_radius)
    {
        index_t cmx = 1;
        index_t range[ORD];
        for ( int d = 0; d < ORD; ++d )
        {
            cmx *= ( 2 * num_cells_radius + 1 );
            range[d] = num_cells_radius;
        }
        checkEdgeTypeRange(cmx);
        int * ccc = new int[ORD*cmx]{0};
        initRectangularGrid(ccc, cmx, range);
        cmx = keepPositiveOffsets(ccc, cmx);
        //printGrid(std::clog, ccc, cmx);
        createRegions(ccc, cmx, range);
        //printRegions(std::clog, "SideRegions");
        delete[] ccc;
    }
    
    
    /// true if createRegions() or createRoundRegions() was called
    bool hasRegions() const
    {
        return ( region_ && border_ );
    }
    
    /// set region array 'offsets' for given cell index
    /**
     A zero offset is always first in the list.
     //\returns the size of the list.
     
     \par Example:

         CELL * cell = & map.icell(indx);
         int n_offset = map.getRegion(offset, indx);
         for ( int n = 0; n < n_offset; ++n )
         {
             Cell & neighbor = cell[offset[n]];
             ...
         }
     
     Note: createRegions() must be called first
    */
    int getRegion(int const*& offsets, const index_t indx) const
    {
        assert_true( hasRegions() );
        int const * R = region_ + chunk_ * (size_t)border_[indx];
        offsets = R + 1;
        return R[0];
    }
    
    /// free memory claimed by the regions
    void deleteRegions()
    {
        delete[] region_;
        region_ = nullptr;
        
        delete[] border_;
        border_ = nullptr;
    }

#pragma mark -

    /// write total number of cells and number of subdivision in each dimension
    void printSummary(std::ostream& os, const char arg[]) const
    {
        os << arg << " of dim " << ORD << " has " << mNbCells << " voxels: ";
        for ( int d = 0; d < ORD; ++d )
        {
            char o = '[', c = ']';
            if ( mPeriodic[d] ) { o = ']'; c = '['; }
            os << " " << o << mInf[d] << ", " << mSup[d];
            os << c << "/" << mDim[d] << " = " << cWidth[d];
        }
        os << std::endl;
    }
    
  
    /// write the list of neigboring cells for each cell
    void printRegions(std::ostream& os, const char arg[]) const
    {
        os << arg << " of dim " << ORD << " has " << mNbCells << " voxels";
        if ( hasRegions() )
        {
            int const* region;
            for ( index_t i = 0; i < mNbCells; ++i )
            {
                int nR = getRegion(region, i);
                os << "\n  cell " << i << " : ";
                for ( int r = 0; r < nR; ++r )
                    os << i+region[r] << " ";
            }
        }
        else
            os << "\n undefined regions\n";
        os << std::endl;
    }
};


#endif
