// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#ifndef SPARMATBLK_H
#define SPARMATBLK_H

#include <cstdio>
#include <iostream>

#include "dim.h"
#include "real.h"
#include "vector.h"
#include "assert_macro.h"

/**
 The block size 'BLOCK_SIZE' can be defined on the command line during compilation,
 and is otherwise set here, to match the dimensionality of the simulation
 */

#define BLOCK_SIZE ( DIM < 3 ? DIM : 3 )

#if ( BLOCK_SIZE == 1 )
#  include "matrix11.h"
#elif ( BLOCK_SIZE == 2 )
#  include "matrix22.h"
#elif ( BLOCK_SIZE == 3 )
#  include "matrix34.h"
#elif ( BLOCK_SIZE == 4 )
#  include "matrix44.h"
#endif


// Flag to enable AVX implementation
#ifdef __AVX__
#  define SPARMATBLK_USES_AVX REAL_IS_DOUBLE
#else
#  define SPARMATBLK_USES_AVX 0
#endif

/// Sparse Matrix with block elements
/**
 The lower triangle of the matrix is stored.
 Elements are stored in no particular order in each column.

 SparMatBlk uses a sparse storage, with arrays of elements for each column.
 Each element is a full square block of size DIM x DIM.
 
 FJN @ Cambridge, August-September 2019
 */
class SparMatBlk final
{
public:
    
#if ( BLOCK_SIZE == 1 )
    typedef Matrix11 Block;
#elif ( BLOCK_SIZE == 2 )
    typedef Matrix22 Block;
#elif ( BLOCK_SIZE == 3 )
    typedef Matrix34 Block;
#elif ( BLOCK_SIZE == 4 )
    typedef Matrix44 Block;
#endif

    /// accessory class
    class Element;
    
    /// A line of the sparse matrix
    class Line
    {
        friend class SparMatBlk;
        friend class Meca;

        Block *blk_;   ///< block elements
        Block *sbk_;   ///< pointer for consolidated elements
        index_t *inx_; ///< column index for each element
        index_t rlen_; ///< number of elements in row
        index_t allo_; ///< allocated size

    public:
        
        /// constructor
        Line() : blk_(nullptr), sbk_(nullptr), inx_(nullptr), rlen_(0), allo_(0) { }
        
        /// the assignment operator will transfer memory
        void operator = (Line&);
        
        /// destructor
        ~Line() { deallocate(); }
        
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
        bool notEmpty() const { return ( rlen_ > 0U ); }
        
        /// return n-th block (not necessarily, located at line inx_[n]
        Block& operator[](index_t n) const { return blk_[n]; }

        /// return block corresponding to index
        Block* find_block(index_t j) const;

        /// return block located at column 'j', allocating if necessary
        Block& block(index_t j);
        
        /// multiplication of a vector: L * X
        void vecMulLine(const real* X, real* Y) const;
        
        /// multiplication of a vector: L * X
        real vecMul1D(const real* X) const;

#if SPARMATBLK_USES_AVX
        /// multiplication of a vector: L * X
        vec2 vecMul2D(const double* X) const;
        
        /// multiplication of a vector: L * X
        vec2 vecMul2DU(const double* X) const;
        
        /// multiplication of a vector: L * X
        vec4 vecMul3DU(const double* X) const;
        
        /// multiplication of a vector: L * X
        vec4 vecMul3DUU(const double* X) const;

        /// multiplication of a vector: L * X
        vec4 vecMul3D(const double* X) const;
        
        /// multiplication of a vector: L * X
        vec4 vecMul4D(const double* X) const;
#endif
    };

private:

    /// create Elements
    static index_t newElements(Element*& ptr, index_t);
    
    /// sort matrix block in increasing index order
    void sortElements();
    
    /// reallocate to use contiguous memory
    void consolidate();
    
    /// copy lower triangle into upper side
    void symmetrize();
    
    /// this is 'true' if symmetrize() was called
    bool is_symmetric;
    
private:
    
    /// number of lines in the matrix, divided by block size
    index_t rsize_;
    
    /// amount of memory which has been allocated
    index_t alloc_;
    
    /// array row_[c][] holds Elements of line 'c'
    Line * row_;
    
    /// colidx_[i] is the index of the first non-empty row of index >= i
    unsigned* colidx_;

    /// memory for consolidated version
    Block * blocks_;
    
public:
    
    /// return the size of the matrix
    index_t size() const { return rsize_ * BLOCK_SIZE; }
    
    /// change the size of the matrix
    void resize(index_t s) { rsize_ = s / BLOCK_SIZE; allocate(rsize_); }

    /// base for destructor
    void deallocate();
    
    /// default constructor
    SparMatBlk();
    
    /// default destructor
    ~SparMatBlk() { deallocate(); }
    
    /// set to zero
    void reset();
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(index_t alc);
    
    /// number of columns
    index_t num_columns() const { return rsize_; }

    /// number of elements in j-th column
    index_t column_size(index_t j) const { assert_true(j<rsize_); return row_[j].rlen_; }
    
    /// reduced line index of n-th element in j-th column (not multiplied by BLOCK_SIZE)
    index_t column_index(index_t j, index_t n) const { assert_true(j<rsize_); return row_[j].inx_[n]; }

    /// returns element stored at line ii and column jj, if ( ii > jj )
    Block& block(const index_t ii, const index_t jj)
    {
        assert_true( ii >= jj );
        assert_true( ii < rsize_ );
        assert_true( jj < rsize_ );
        assert_false( is_symmetric );
        return row_[ii].block(jj);
    }
    
    /// returns element stored at line ii and column jj, if ( ii > jj )
    Block& diag_block(const index_t ii)
    {
        assert_true( ii < rsize_ );
        assert_false( is_symmetric );
        return row_[ii].block(ii);
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
    void vecMulAdd2D(const real* X, real* Y, index_t start, index_t stop) const;

    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd3D(const real* X, real* Y, index_t start, index_t stop) const;
    
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd_ALT(const real* X, real* Y, index_t start, index_t stop) const;

    
    /// multiplication of a vector, for columns within [start, stop[
    void vecMul(const real*, real* Y, index_t start, index_t stop) const;

    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMul2D(const real* X, real* Y, index_t start, index_t stop) const;

    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMul3D(const real* X, real* Y, index_t start, index_t stop) const;

    
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd(const real* X, real* Y) const { vecMulAdd(X, Y, 0, rsize_); }
    
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd_ALT(const real* X, real* Y) const { vecMulAdd_ALT(X, Y, 0, rsize_); }

    /// 2D isotropic multiplication (not implemented)
    void vecMulAddIso2D(const real* X, real* Y) const {};
    
    /// 3D isotropic multiplication (not implemented)
    void vecMulAddIso3D(const real*, real*) const {};
    
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M)
    void vecMul(const real* X, real* Y) const { vecMul(X, Y, 0, rsize_); }

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

    /// print content of one column
    void printSummary(std::ostream&);
    
    /// print
    void printBlocks(std::ostream&) const;

    /// debug function
    int bad() const;
};


#endif

