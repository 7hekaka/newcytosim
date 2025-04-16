// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Francois Nedelec; Last updated Jan 2008. nedelec@embl.de

#ifndef GRID_H
#define GRID_H

#include "map.h"


/// A Map with an instantiation of class CELL in each voxel
/** 
 Most of the functionality is provided by the parent class Map.
 This in addition provides functions to interpolate numerical values,
 and to clear the values of the cells.
*/
///\todo add Grid<> copy constructor and copy assignment

template <typename CELL, int ORD>
class Grid : public Map<ORD>
{
    
    /// Disabled copy constructor
    Grid<CELL, ORD>(Grid<CELL, ORD> const&);
    
    /// Disabled copy assignment
    Grid<CELL, ORD>& operator = (Grid<CELL, ORD> const&);

protected:
    
    /// The array of pointers to cells
    CELL * gCell;
    
    /// allocated size of array gCells[]
    size_t gAllocated;

public:
    
    /// Type of the parent class
    typedef Map<ORD> MAP;

    /// The type of cells (=CELL)
    typedef CELL value_type;
    
    /// constructor
    Grid()
    {
        gAllocated = 0;
        gCell = nullptr;
    }
    
    /// Free memory
    void destroy()
    {
        deleteCells();
        MAP::destroy();
    }
    
    /// Destructor
    virtual ~Grid()
    {
        destroy();
    }
    
    /// allocate the array of cells
    size_t createCells()
    {
        if ( MAP::mNbCells == 0 )
            printf("Map has no cells: call Map::setDimensions() first\n");
        
        if ( MAP::mNbCells > gAllocated )
        {
            delete[] gCell;
            // allocate one more than necessary, and use a multiple of 4:
            gAllocated = ( MAP::mNbCells + 4 ) & ~3UL;
            //printf("Grid %p allocated %lu cells\n", this, gAllocated);
            gCell = new CELL[gAllocated];
            return gAllocated;
        }
        return 0;
    }
    
    /// returns true if cells have been allocated
    size_t hasCells() const
    {
        if ( gCell )
            return gAllocated;
        return 0;
    }
    
    /// deallocate array of cells
    void deleteCells()
    {
        delete[] gCell;
        gCell = nullptr;
        gAllocated = 0;
    }
    
    /// call function clear() for all cells
    void clearCells()
    {
        if ( !gCell )
            return;
        
        CELL * c = gCell;
        const CELL * end = gCell + MAP::mNbCells;
#if ( 1 )
        //we unroll the loop for speed
        const CELL * stop = end - 7;
        for ( ; c < stop; c += 8 )
        {
            c[0].clear();
            c[1].clear();
            c[2].clear();
            c[3].clear();
            c[4].clear();
            c[5].clear();
            c[6].clear();
            c[7].clear();
        }
#endif
        for ( ; c < end; ++c )
            c->clear();
    }

    //--------------------------------------------------------------------------

    /// address of the cell array ( equivalent to &cell(0) )
    CELL * data() const
    {
        return gCell;
    }
    
    /// return cell at index 'inx'
    CELL& icell(const index_t inx) const
    {
        assert_true( inx < MAP::mNbCells );
        return gCell[ inx ];
    }
    
    /// reference to CELL whose center is closest to w[]
    CELL& cell(const real w[ORD]) const
    {
        index_t inx = MAP::index(w);
        assert_true( inx < MAP::mNbCells );
        return gCell[ inx ];
    }
    
    /// reference to CELL of coordinates c[]
    CELL& cell(const int c[ORD]) const
    {
        assert_true( MAP::pack(c) < MAP::mNbCells );
        return gCell[ MAP::pack(c) ];
    }
   
    /// operator access to a cell by index
    CELL& operator[](const index_t inx) const
    {
        assert_true( inx <= MAP::mNbCells );
        return gCell[ inx ];
    }

    /// short-hand access to a cell by coordinates
    CELL& operator()(const int c[ORD]) const
    {
        assert_true( MAP::pack(c) < MAP::mNbCells );
        return gCell[ MAP::pack(c) ];
    }
    
    /// operator access to a cell by position
    CELL& operator()(const real w[ORD]) const
    {
        assert_true( MAP::index(w) < MAP::mNbCells );
        return gCell[ MAP::index(w) ];
    }
    
    
    //--------------------------------------------------------------
#pragma mark -
    
    /// create a 1D-map
    void create1D(real i, real s, index_t d)
    {
        assert_true( ORD == 1 );
        MAP::setDimensions(&i, &s, &d);
        createCells();
    }

    /// access to cell for ORD==1
    CELL& icell1D(const int x) const
    {
        return gCell[MAP::pack1D(x)];
    }
    
    /// access to cell for ORD==2
    CELL& icell2D(const int x, const int y) const
    {
        return gCell[MAP::pack2D(x,y)];
    }
    
    /// access to cell for ORD==3
    CELL& icell3D(const int x, const int y, const int z) const
    {
        return gCell[MAP::pack3D(x,y,z)];
    }
    
    /// access to cell for ORD==1
    CELL& icell1D_clamped(const int x) const
    {
        index_t inx = MAP::pack1D_clamped(x);
        assert_true( inx < MAP::mNbCells );
        return gCell[inx];
    }
    
    /// access to cell for ORD==2
    CELL& icell2D_clamped(const int x, const int y) const
    {
        index_t inx = MAP::pack2D_clamped(x,y);
        assert_true( inx < MAP::mNbCells );
        return gCell[inx];
    }
    
    /// access to cell for ORD==3
    CELL& icell3D_clamped(const int x, const int y, const int z) const
    {
        index_t inx = MAP::pack3D_clamped(x,y,z);
        assert_true( inx < MAP::mNbCells );
        return gCell[inx];
    }

    //-----------------------------------------------------------------------
#pragma mark - Interpolate
    
    /// return linear interpolation of values stored at the center of each cell
    CELL interpolate( const real w[ORD] ) const
    {
        //we have 2^ORD corner cells
        const int sz = 1 << ORD;
        index_t inx[sz];   //incides of the corner cells
        real    alp[sz];   //coefficients of interpolation
        
        int nb = 0;
        for ( int d = ORD-1; d >= 0; --d )
        {
            real a = MAP::mapC(d, w[d]);
            int ia = (int)std::floor(a);
            a -= ia;
            int  l = MAP::ind_(d, ia-1);
            int  u = MAP::ind_(d, ia  );
            
            if ( nb == 0 )
            {
                //initialize the edges ( l and u ) and appropriate coefficients
                inx[1] = u;  alp[1] = a;
                inx[0] = l;  alp[0] = 1-a;
                nb = 2;
            }
            else
            {
                //double the amount of edges at each round,
                //with the indices and coefficients for lower and upper bounds
                for ( int c = 0; c < nb; ++c )
                {
                    inx[c+nb] = MAP::breadth(d) * inx[c] + u;
                    alp[c+nb] = alp[c] * a;
                    inx[c   ] = MAP::breadth(d) * inx[c] + l;
                    alp[c   ] = alp[c] * (1-a);
                }
                nb *= 2;
            }
        }
        assert_true( nb == sz );
        
        //sum weighted cells to interpolate
        CELL res = 0;
        for ( int c = 0; c < sz; ++c ) 
            res += alp[c] * gCell[inx[c]];
        return res;
    }
    
    
    /// return linear interpolation of values stored at the center of each cell, if ORD==1
    CELL interpolate1D( const real xx ) const
    {
        assert_true( ORD == 1 );
        
        real ax = MAP::mapC(0, xx);
        
#if ENABLE_PERIODIC_BOUNDARIES
        int ix = (int)std::floor(ax);
#else
        int ix = (int)ax;
#endif
        
        ax -= ix;
        
        index_t lx = MAP::ind_(0, ix-1);
        index_t ux = MAP::ind_(0, ix  );
        
        return gCell[lx] + ax * ( gCell[ux] - gCell[lx] );
    }

    
    /// return linear interpolation of values stored at the center of each cell, if ORD==2
    CELL interpolate2D( const real w[ORD] ) const
    {
        assert_true( ORD == 2 );
        
        real ax = MAP::mapC(0, w[0]);
        real ay = MAP::mapC(1, w[1]);
        
#if ENABLE_PERIODIC_BOUNDARIES
        int ix = (int)std::floor(ax);
        int iy = (int)std::floor(ay);
#else
        int ix = (int)ax;
        int iy = (int)ay;
#endif
        
        ax -= ix;
        ay -= iy;
        
        index_t ly = MAP::ind_(1, iy-1) * MAP::breadth(0);
        index_t uy = MAP::ind_(1, iy  ) * MAP::breadth(0);

        index_t lx = MAP::ind_(0, ix-1);
        index_t ux = MAP::ind_(0, ix  );
        
        //sum weighted cells to get interpolation
        CELL  rl = gCell[lx+ly] + ay * ( gCell[lx+uy] - gCell[lx+ly] );
        CELL  ru = gCell[ux+ly] + ay * ( gCell[ux+uy] - gCell[ux+ly] );

        return rl + ax * ( ru - rl );
    }

    
    /// return linear interpolation of values stored at the center of each cell, if ORD==3
    CELL interpolate3D( const real w[ORD] ) const
    {
        assert_true( ORD == 3 );
        
        real ax = MAP::mapC(0, w[0]);
        real ay = MAP::mapC(1, w[1]);
        real az = MAP::mapC(2, w[2]);
        
#if ENABLE_PERIODIC_BOUNDARIES
        int ix = (int)std::floor(ax);
        int iy = (int)std::floor(ay);
        int iz = (int)std::floor(az);
#else
        int ix = (int)ax;
        int iy = (int)ay;
        int iz = (int)az;
#endif
        
        ax -= ix;
        ay -= iy;
        az -= iz;
        
        index_t lz = MAP::ind_(2, iz-1) * MAP::breadth(1) * MAP::breadth(0);
        index_t uz = MAP::ind_(2, iz  ) * MAP::breadth(1) * MAP::breadth(0);
        
        index_t ly = MAP::ind_(1, iy-1) * MAP::breadth(0);
        index_t uy = MAP::ind_(1, iy  ) * MAP::breadth(0);

        index_t lx = MAP::ind_(0, ix-1);
        index_t ux = MAP::ind_(0, ix  );

        CELL * cul = gCell + (ux+ly), rul = cul[lz] + ( cul[uz] - cul[lz] ) * az;
        CELL * cuu = gCell + (ux+uy), ruu = cuu[lz] + ( cuu[uz] - cuu[lz] ) * az;
        CELL * cll = gCell + (lx+ly), rll = cll[lz] + ( cll[uz] - cll[lz] ) * az;
        CELL * clu = gCell + (lx+uy), rlu = clu[lz] + ( clu[uz] - clu[lz] ) * az;

        CELL ru = rul + ( ruu - rul ) * ay;
        CELL rl = rll + ( rlu - rll ) * ay;

        return rl + ( ru - rl ) * ax;
    }


    //--------------------------------------------------------------------------
#pragma mark - For Numerical Cells

    /// set all cells to `val`
    void setValues(const CELL val)
    {
        assert_true( MAP::mNbCells <= gAllocated );
        for ( index_t i = 0; i < MAP::mNbCells; ++i )
            gCell[i] = val;
    }
    
    /// multiply all cells by `val`
    void scaleValues(const CELL val)
    {
        assert_true ( MAP::mNbCells <= gAllocated );
        for ( index_t i = 0; i < MAP::mNbCells; ++i )
            gCell[i] *= val;
    }
    
    /// get sum, minimum and maximum value over all cells
    void infoValues(CELL& sum, CELL& avg, CELL& mn, CELL& mx) const
    {
        assert_true( MAP::mNbCells <= gAllocated );
        sum = 0;
        mn = gCell[0];
        mx = gCell[0];
        for ( index_t i = 0; i < MAP::mNbCells; ++i )
        {
            //mn = std::min(mn, gCell[i]);
            //mx = std::max(mx, gCell[i]);
            sum += gCell[i];
        }
        avg = sum * ( 1.0 / MAP::mNbCells );
    }

    /// sum of all values, if CELL supports the addition
    CELL sumValues() const
    {
        assert_true( MAP::mNbCells <= gAllocated );
        CELL result = 0;
        for ( index_t i = 0; i < MAP::mNbCells; ++i )
            result += gCell[i];
        return result;
    }
    
    /// maximum value over all cells
    CELL maxValue() const
    {
        assert_true( MAP::mNbCells <= gAllocated );
        CELL res = gCell[0];
        for ( index_t i = 1; i < MAP::mNbCells; ++i )
            res.e_min(gCell[i]);
        return res;
    }
    
    /// minimum value over all cells
    CELL minValue() const
    {
        assert_true( MAP::mNbCells <= gAllocated );
        CELL res = gCell[0];
        for ( index_t i = 1; i < MAP::mNbCells; ++i )
            res.e_max(gCell[i]);
        return res;
    }
    
    /// true if any( cells[] < 0 )
    bool hasNegativeValue() const
    {
        assert_true( MAP::mNbCells <= gAllocated );
        for ( index_t i = 0; i < MAP::mNbCells; ++i )
            if ( gCell[i].negative() )
                return true;
        return false;
    }
    
#pragma mark -
    
    /// the sum of the values in the region around cell referred by 'inx'
    CELL sumValuesInRegion(const index_t inx) const
    {
        int const* offsets = nullptr;
        const CELL * ce = gCell + inx;
        int nR = MAP::getRegion(offsets, inx);
        CELL result = ce[0];
        for ( int c = 0; c < nR; ++c )
            result += ce[ offsets[c] ];
        return result;
    }
    
    /// the sum of the values in the region around cell referred by 'inx'
    CELL avgValueInRegion(const index_t inx) const
    {
        int const* offsets = nullptr;
        const CELL * ce = gCell + inx;
        int nR = MAP::getRegion(offsets, inx);
        CELL result = ce[0];
        for ( int c = 0; c < nR; ++c )
            result += ce[ offsets[c] ];
        return result / (real)(nR+1);
    }
    
    /// the maximum of the values in the region around cell referred by 'inx'
    CELL maxValueInRegion(const index_t inx) const
    {
        assert_true( MAP::mNbCells <= gAllocated );
        int const* offsets = nullptr;
        const CELL * ce = gCell + inx;
        int nR = MAP::getRegion(offsets, inx);
        CELL result = gCell[inx];
        for ( int c = 0; c < nR; ++c )
            result = std::max(result, ce[ offsets[c] ]);
        return result;
    }
    
#pragma mark -
    
    /// write values to a file, with the position for each cell (file can be stdout)
    void printValues(FILE* file, const real offset) const
    {
        assert_true( MAP::mNbCells <= gAllocated );
        real w[ORD];
        for ( index_t i = 0; i < MAP::mNbCells; ++i )
        {
            MAP::setPositionFromIndex(w, i, offset);
            for ( int d=0; d < ORD; ++d )
                fprintf(file, "%7.2f ", w[d]);
            fprintf(file,"  %f\n", gCell[i]);
        }
    }
    
    /// write values to a file, with the range for each cell (file can be stdout)
    void printValuesWithRange(FILE* file) const
    {
        assert_true( MAP::mNbCells <= gAllocated );
        real l[ORD], r[ORD];
        for ( index_t i = 0; i < MAP::mNbCells; ++i )
        {
            MAP::setPositionFromIndex(l, i, 0.0);
            MAP::setPositionFromIndex(r, i, 1.0);
            for ( int d=0; d < ORD; ++d )
                fprintf(file, "%7.2f %7.2f  ", l[d], r[d]);
            fprintf(file,"  %f\n", gCell[i]);
        }
    }
};

#endif
