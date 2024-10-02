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

        Block* blk_;    ///< all blocks
        unsigned* inx_; ///< line index for each element
        unsigned alo_;    ///< allocated size of array
        unsigned nbb_;    ///< number of blocks in column

    public:
        
        /// constructor
        Column();
        
        /// the assignment operator will transfer memory
        void operator = (Column&);
        
        /// destructor
        ~Column() { deallocate(); }
        
        /// allocate to hold 'nb' elements
        void allocate(index_t nb);
        
        /// deallocate memory
        void deallocate();

        /// set as zero
        void reset();
        
        /// sort element by increasing indices, using given temporary array
        void sortElements(Element[], index_t);
        
        /// print
        void printBlocks(std::ostream&) const;
        
        /// true if column is empty
        bool notEmpty() const { return ( nbb_ > 0 ); }
        
        /// return n-th block (not necessarily, located at line inx_[n]
        Block& operator[](index_t n) const { return blk_[n]; }

        /// return block corresponding to index
        Block* find_block(index_t j) const;

        /// return block located at line 'i' and column 'j'
        Block& block(index_t i, index_t j);
        
        /// multiplication of a vector: Y <- Y + M * X, block_size = 1
        void vecMulAdd1D(const real* X, real* Y, index_t j) const;
        /// multiplication of a vector: Y <- Y + M * X, block_size = 2
        void vecMulAdd2D(const real* X, real* Y, index_t j) const;
        /// multiplication of a vector: Y <- Y + M * X, block_size = 3
        void vecMulAdd3D(const real* X, real* Y, index_t j) const;
        /// multiplication of a vector: Y <- Y + M * X, block_size = 4
        void vecMulAdd4D(const real* X, real* Y, index_t j) const;

        
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 2
        void vecMulAdd2D_SSE(const double* X, double* Y, index_t j) const;
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 2
        void vecMulAdd2D_AVX(const double* X, double* Y, index_t j) const;
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 2
        void vecMulAdd2D_AVXU(const double* X, double* Y, index_t j) const;
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 2
        void vecMulAdd2D_AVXUU(const double* X, double* Y, index_t j) const;
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 3
        void vecMulAdd3D_SSE(const float* X, float* Y, index_t j) const;
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 3
        void vecMulAdd3D_SSEU(const float* X, float* Y, index_t j) const;
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 3
        void vecMulAdd3D_AVX(const double* X, double* Y, index_t j) const;
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 3
        void vecMulAdd3D_AVXU(const double* X, double* Y, index_t j) const;
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 4
        void vecMulAdd4D_AVX(const double* X, double* Y, index_t j) const;
    };

private:

    /// create Elements
    static index_t newElements(Element*& ptr, index_t);
    
    /// sort matrix block in increasing index order
    void sortElements();

private:
    
    /// size of matrix
    index_t rsize_;
    
    /// amount of memory which has been allocated
    index_t alloc_;

    /// array col_[c][] holds Elements of column 'c'
    Column * column_;
    
    /// colidx_[i] is the index of the first non-empty column of index >= i
    unsigned * colidx_;

public:
    
    /// return the size of the matrix
    index_t size() const { return rsize_ * S_BLOCK_SIZE; }
    
    /// change the size of the matrix
    void resize(index_t s) { rsize_ = s / S_BLOCK_SIZE; allocate(rsize_); }

    /// base for destructor
    void deallocate();
    
    /// default constructor
    SparMatSymBlk();
    
    /// default destructor
    ~SparMatSymBlk()  { deallocate(); }
    
    /// set to zero
    void reset();
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(index_t alc);
    
    /// number of columns
    index_t num_columns() const { return rsize_; }

    /// number of elements in j-th column
    index_t column_size(index_t j) const { assert_true(j<rsize_); return column_[j].nbb_; }
    
    /// line index of n-th element in j-th column (not multiplied by BLOCK_SIZE)
    index_t column_index(index_t j, index_t n) const { assert_true(j<rsize_); return column_[j].inx_[n]; }

    /// returns element at (i, i)
    Block& diag_block(index_t i);

    /// returns element stored at line ii and column jj, if ( ii > jj )
    Block& block(const index_t ii, const index_t jj)
    {
        assert_true( ii < rsize_ );
#if ( 1 )
        assert_true( ii >= jj );
        return column_[jj].block(ii, jj);
#else
        assert_true( jj < rsize_ );
        index_t i = std::max(ii, jj);
        index_t j = std::min(ii, jj);
        return column_[j].block(i, j);
#endif
    }
    
    /// returns the address of element at line i, column j, no allocation is done
    real* address(index_t i, index_t j) const;

    /// returns the address of element at line i, column j, allocating if necessary
    real& element(index_t i, index_t j);

    /// returns the address of element at line i, column j, allocating if necessary
    real& operator()(index_t i, index_t j) { return element(i,j); }

    /// scale the matrix by a scalar factor
    void scale(real);
    
    /// add terms with `i` and `j` in [start, start+cnt[ to `mat`
    void addDiagonalBlock(real* mat, index_t ldd, index_t start, index_t cnt, index_t mul) const;
    
    /// add scaled terms with `i` in [start, start+cnt[ if ( j > i ) and ( j <= i + rank ) to `mat`
    void addLowerBand(real alpha, real* mat, index_t ldd, index_t start, index_t cnt, index_t mul, index_t rank) const;

    /// add `alpha*trace()` for blocks within [start, start+cnt[ if ( j <= i + rank ) to `mat`
    void addDiagonalTrace(real alpha, real* mat, index_t ldd, index_t start, index_t cnt, index_t mul, index_t rank, bool sym) const;
    
    
    /// prepare matrix for multiplications by a vector (must be called)
    bool prepareForMultiply(int);

    /// multiplication of a vector, for columns within [start, stop[
    void vecMulAdd(const real*, real* Y, index_t start, index_t stop) const;
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd(const real* X, real* Y) const { vecMulAdd(X, Y, 0, rsize_); }

    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd_ALT(const real* X, real* Y, index_t start, index_t stop) const;
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd_ALT(const real* X, real* Y) const { vecMulAdd_ALT(X, Y, 0, rsize_); }
    
    /// 2D isotropic multiplication (not implemented)
    void vecMulAddIso2D(const real* X, real* Y) const {};
    /// 3D isotropic multiplication (not implemented)
    void vecMulAddIso3D(const real*, real*) const {};
    
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M)
    void vecMul(const real* X, real* Y) const;

    /// true if matrix is non-zero
    bool notZero() const;
    
    /// number of blocks in columns [start, stop[. Set allocated size
    size_t nbElements(index_t start, index_t stop, size_t& alc) const;
    
    /// total number of blocks currently in use
    size_t nbElements() const { size_t alc=0; return nbElements(0, rsize_, alc); }

    /// returns a string which a description of the type of matrix
    std::string what() const;
    
    /// print matrix columns in sparse mode: ( i, j : value ) if |value| >= inf
    void printSparse(std::ostream&, real inf, index_t start, index_t stop) const;
    
    /// print matrix in sparse mode: ( i, j : value ) if |value| >= inf
    void printSparse(std::ostream& os, real inf) const { printSparse(os, inf, 0, rsize_); }

    /// print size of columns
    void printSummary(std::ostream&, index_t start, index_t stop);
    
    /// print
    void printBlocks(std::ostream&) const;

    /// debug function
    int bad() const;
};


#endif

