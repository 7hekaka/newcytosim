// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#ifndef SPARMATSYMBLK_H
#define SPARMATSYMBLK_H

#include "dim.h"
#include "real.h"
#include <cstdio>
#include <iostream>
#include "assert_macro.h"


/**
 The block size 'S_BLOCK_SIZE' can be defined on the command line during compilation,
 and is otherwise set here, to match the dimensionality of the simulation
 */

#define S_BLOCK_SIZE ( DIM < 3 ? DIM : 3 )

#if ( S_BLOCK_SIZE == 1 )
#   include "matrix11.h"
#elif ( S_BLOCK_SIZE == 2 )
#   include "matrix22.h"
#elif ( S_BLOCK_SIZE == 3 )
#   include "matrix33.h"
#elif ( S_BLOCK_SIZE == 4 )
#   include "matrix44.h"
#endif

///real symmetric sparse Matrix
/**
 The lower triangle of the matrix is stored.
 Elements are stored in no particular order in each column.

 SparMatSymBlk uses a sparse storage, with arrays of elements for each column.
 Each element is a full square block of size DIM x DIM.
 
 F. Nedelec, 17--27 March 2017, revised entirely June 2018
 */
class SparMatSymBlk final
{
public:

#if ( S_BLOCK_SIZE == 1 )
    typedef Matrix11 Block;
#elif ( S_BLOCK_SIZE == 2 )
    typedef Matrix22 Block;
#elif ( S_BLOCK_SIZE == 3 )
    typedef Matrix33 Block;
#elif ( S_BLOCK_SIZE == 4 )
    typedef Matrix44 Block;
#endif

    /// accessory class used to sort columns
    class Element;
    
    /// A column of the sparse matrix
    class Column
    {
        friend class SparMatSymBlk;
        friend class Meca;

        size_t  allo_;   ///< allocated size of array
        size_t   nbb_;   ///< number of blocks in column
        size_t * inx_;   ///< line index for each element
        Block  * blk_;   ///< all blocks
        
    public:
        
        /// constructor
        Column();
        
        /// the assignment operator will transfer memory
        void operator =(Column&);
        
        /// destructor
        ~Column() { deallocate(); }
        
        /// allocate to hold 'nb' elements
        void allocate(size_t nb);
        
        /// deallocate memory
        void deallocate();

        /// set as zero
        void reset();
        
        /// sort element by increasing indices, using given temporary array
        void sortElements(Element[], size_t);
        
        /// print
        void printBlocks(std::ostream&) const;
        
        /// true if column is empty
        bool isNotZero() const { return ( nbb_ > 0 ); }

        /// return n-th block (not necessarily, located at line inx_[n]
        Block& operator[](size_t n) const { return blk_[n]; }

        /// return block located at line 'i' and column 'j'
        Block& block(size_t i, size_t j);
        
        /// multiplication of a vector: Y <- Y + M * X, block_size = 1
        void vecMulAdd1D(const real* X, real* Y, size_t j) const;
        /// multiplication of a vector: Y <- Y + M * X, block_size = 2
        void vecMulAdd2D(const real* X, real* Y, size_t j) const;
        /// multiplication of a vector: Y <- Y + M * X, block_size = 3
        void vecMulAdd3D(const real* X, real* Y, size_t j) const;
        /// multiplication of a vector: Y <- Y + M * X, block_size = 4
        void vecMulAdd4D(const real* X, real* Y, size_t j) const;

        
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 2
        void vecMulAdd2D_SSE(const double* X, double* Y, size_t j) const;
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 2
        void vecMulAdd2D_AVX(const double* X, double* Y, size_t j) const;
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 2
        void vecMulAdd2D_AVXU(const double* X, double* Y, size_t j) const;
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 2
        void vecMulAdd2D_AVXUU(const double* X, double* Y, size_t j) const;
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 3
        void vecMulAdd3D_SSE(const float* X, float* Y, size_t j) const;
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 3
        void vecMulAdd3D_SSEU(const float* X, float* Y, size_t j) const;
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 3
        void vecMulAdd3D_AVX(const double* X, double* Y, size_t j) const;
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 3
        void vecMulAdd3D_AVXU(const double* X, double* Y, size_t j) const;
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 4
        void vecMulAdd4D_AVX(const double* X, double* Y, size_t j) const;
    };

private:

    /// create Elements
    static size_t newElements(Element*& ptr, size_t);
    
    /// sort matrix block in increasing index order
    void sortElements();

private:
    
    /// size of matrix
    size_t   size_;
    
    /// amount of memory which has been allocated
    size_t   alloc_;

    /// array col_[c][] holds Elements of column 'c'
    Column * column_;
    
    /// colidx_[i] is the index of the first non-empty column of index >= i
    size_t * colidx_;

public:
    
    /// return the size of the matrix
    size_t size() const { return size_; }
    
    /// change the size of the matrix
    void resize(size_t s) { allocate(s); size_=s; }

    /// base for destructor
    void deallocate();
    
    /// default constructor
    SparMatSymBlk();
    
    /// default destructor
    ~SparMatSymBlk()  { deallocate(); }
    
    /// set to zero
    void reset();
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(size_t alc);
    
    /// number of elements in j-th column
    size_t column_size(size_t j) const { assert_true(j<size_); return column_[j].nbb_; }
    
    /// line index of n-th element in j-th column
    size_t column_index(size_t j, size_t n) const { return column_[j].inx_[n]; }

    /// returns element at (i, i)
    Block& diag_block(size_t i);

    /// returns element stored at line ii and column jj, if ( ii > jj )
    Block& block(const size_t ii, const size_t jj)
    {
        assert_true( ii < size_ );
        assert_true( jj < size_ );
        assert_true( ii % S_BLOCK_SIZE == 0 );
        assert_true( jj % S_BLOCK_SIZE == 0 );
#if ( 0 )
        // safe swap, with branchless code:
        size_t i = std::max(ii, jj);
        size_t j = std::min(ii, jj);
        return column_[j].block(i, j);
#else
        assert_true( ii >= jj );
        return column_[jj].block(ii, jj);
#endif
    }
    
    /// returns the address of element at (x, y), no allocation is done
    real* addr(size_t x, size_t y) const;
    
    /// returns the address of element at (x, y), allocating if necessary
    real& operator()(size_t x, size_t y);
    
    /// scale the matrix by a scalar factor
    void scale(real);
    
    /// add terms with `i` and `j` in [start, start+cnt[ to `mat`
    void addDiagonalBlock(real* mat, size_t ldd, size_t start, size_t cnt) const;
    
    /// add scaled terms with `i` in [start, start+cnt[ if ( j > i ) and ( j <= i + rank ) to `mat`
    void addLowerBand(real alpha, real* mat, size_t ldd, size_t start, size_t cnt, size_t rank) const;

    /// add `alpha*trace()` for blocks within [start, start+cnt[ if ( j <= i + rank ) to `mat`
    void addDiagonalTrace(real alpha, real* mat, size_t ldd, size_t start, size_t cnt, size_t rank, bool sym) const;
    
    
    /// prepare matrix for multiplications by a vector (must be called)
    bool prepareForMultiply(int);

    /// multiplication of a vector, for columns within [start, stop[
    void vecMulAdd(const real*, real* Y, size_t start, size_t stop) const;
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd_ALT(const real* X, real* Y, size_t start, size_t stop) const;
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd(const real* X, real* Y) const { vecMulAdd(X, Y, 0, size_); }
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd_TIME(const real* X, real* Y) const;
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd_ALT(const real* X, real* Y) const;
    
    /// 2D isotropic multiplication (not implemented)
    void vecMulAddIso2D(const real* X, real* Y) const {};
    /// 3D isotropic multiplication (not implemented)
    void vecMulAddIso3D(const real*, real*) const {};
    
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M)
    void vecMul(const real* X, real* Y) const;

    /// true if matrix is non-zero
    bool isNotZero() const;
    
    /// number of blocks in columns [start, stop[. Set allocated size
    size_t nbElements(size_t start, size_t stop, size_t& alc) const;
    
    /// total number of blocks currently in use
    size_t nbElements() const { size_t alc=0; return nbElements(0, size_, alc); }

    /// returns a string which a description of the type of matrix
    std::string what() const;
    
    /// print matrix columns in sparse mode: ( i, j : value ) if |value| >= inf
    void printSparse(std::ostream&, real inf, size_t start, size_t stop) const;
    
    /// print matrix in sparse mode: ( i, j : value ) if |value| >= inf
    void printSparse(std::ostream& os, real inf) const { printSparse(os, inf, 0, size_); }

    /// print size of columns
    void printSummary(std::ostream&, size_t start, size_t stop);
    
    /// print
    void printBlocks(std::ostream&) const;

    /// debug function
    int bad() const;
};


#endif

