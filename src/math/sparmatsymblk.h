// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#ifndef SPARMATSYMBLK_H
#define SPARMATSYMBLK_H

#include "dim.h"
#include "real.h"
#include <cstdio>
#include <iostream>
#include "assert_macro.h"


/**
 The block size 'BLOCK_SIZE' can be defined on the command line during compilation,
 and is otherwise set here, to match the dimensionality of the simulation
 */

#define BLOCK_SIZE DIM

#if ( BLOCK_SIZE == 1 )
#   include "matrix11.h"
#elif ( BLOCK_SIZE == 2 )
#   include "matrix22.h"
#elif ( BLOCK_SIZE == 3 )
#   include "matrix33.h"
#elif ( BLOCK_SIZE == 4 )
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

#if ( BLOCK_SIZE == 1 )
    typedef Matrix11 Block;
#elif ( BLOCK_SIZE == 2 )
    typedef Matrix22 Block;
#elif ( BLOCK_SIZE == 3 )
    typedef Matrix33 Block;
#elif ( BLOCK_SIZE == 4 )
    typedef Matrix44 Block;
#endif

    /// accessory class used to sort columns
    class Element;
    
    /// number of real in a block
    static constexpr size_t SB = sizeof(Block) / sizeof(real);

private:
    
    /// A column of the sparse matrix
    class Column
    {
        friend class SparMatSymBlk;
        friend class Meca;

        size_t  allo_;   ///< allocated size of array
        size_t  size_;   ///< number of blocks in column
        size_t * inx_;   ///< line index for each element
        Block  * blk_;   ///< all blocks
        
    public:
        
        /// constructor
        Column() { size_=0; allo_=0; inx_=nullptr; blk_=nullptr; }
        
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
        void print(std::ostream&) const;
        
        /// true if column is empty
        bool isNotZero() const { return ( size_ > 0 ); }

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
        void vecMulAdd2D_SSE(const real* X, real* Y, size_t j) const;
        
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 2
        void vecMulAdd2D_AVX(const real* X, real* Y, size_t j) const;
        
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 2
        void vecMulAdd2D_AVXU(const real* X, real* Y, size_t j) const;
        
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 2
        void vecMulAdd2D_AVXUU(const real* X, real* Y, size_t j) const;
        
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 3
        void vecMulAdd3D_SSE(const real* X, real* Y, size_t j) const;
        
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 3
        void vecMulAdd3D_SSEU(const real* X, real* Y, size_t j) const;

        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 3
        void vecMulAdd3D_AVX(const real* X, real* Y, size_t j) const;

        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 3
        void vecMulAdd3D_AVXU(const real* X, real* Y, size_t j) const;

        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 4
        void vecMulAdd4D_AVX(const real* X, real* Y, size_t j) const;
    };
    
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
    
    /// next_[ii] is the index of the first non-empty column of index >= ii
    size_t * next_;

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
    
    /// return column at index j
    Column const& column(size_t j) const { return column_[j]; }
    
    /// returns element at (i, i)
    Block& diag_block(size_t i);

    /// returns element stored at line ii and column jj, if ( ii > jj )
    Block& block(const size_t ii, const size_t jj)
    {
        assert_true( ii < size_ );
        assert_true( jj < size_ );
        assert_true( ii % BLOCK_SIZE == 0 );
        assert_true( jj % BLOCK_SIZE == 0 );
#if ( 1 )
        // safe swap, with branchless code:
        size_t i = std::max(ii, jj);
        size_t j = std::min(ii, jj);
        return column_[j].block(i, j);
#else
        assert_true( ii > jj );
        return column_[jj].block(ii, jj);
#endif
    }
    
    /// returns the address of element at (x, y), no allocation is done
    real* addr(size_t x, size_t y) const;
    
    /// returns the address of element at (x, y), allocating if necessary
    real& operator()(size_t x, size_t y);
    
    /// scale the matrix by a scalar factor
    void scale(real);
    
    /// add the diagonal block ( start, start+nb ) to `mat`
    void addDiagonalBlock(real* mat, size_t ldd, size_t start, size_t nb) const;
    
    /// add `alpha*trace()` for sub blocks within ( start, start+nb ) to `mat`
    void addDiagonalTrace(real alpha, real* mat, size_t ldd, size_t start, size_t nb) const;
    
    /// add `alpha*trace()` for sub blocks within ( start, start+nb ) to `mat`
    void addDiagonalTraceBanded(real alpha, real* mat, size_t ldd, size_t start, size_t nb, size_t rank) const;

    /// add upper triangular half of 'this' block ( start, start+nb ) to `mat`
    void addTriangularBlock(real* mat, size_t ldd, size_t start, size_t nb) const;
    
    
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
    
    /// number of blocks in columns within [start, stop[
    size_t nbElements(size_t start, size_t stop) const;
    
    /// number of blocks which are not null
    size_t nbElements() const { return nbElements(0, size_); }

    /// returns a string which a description of the type of matrix
    std::string what() const;
    
    /// printf debug function in sparse mode: i, j : value
    void printSparse(std::ostream&, real) const;

    /// print size of columns
    void printColumns(std::ostream&);
    
    /// debug function
    int bad() const;
};


#endif

