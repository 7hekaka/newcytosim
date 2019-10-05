// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MATSPARSEBLK_H
#define MATSPARSEBLK_H

#include "real.h"
#include <cstdio>
#include <iostream>
#include "assert_macro.h"

/**
 The block size 'BLOCK_SIZE' can be defined on the command line during compilation,
 and is otherwise set here, to match the dimensionality of the simulation
 */

#ifndef BLOCK_SIZE
#  include "dim.h"
#  define BLOCK_SIZE DIM
#endif

#include "vector.h"

#if ( BLOCK_SIZE == 1 )
#  include "matrix11.h"
typedef Matrix11 SubBlock;
#elif ( BLOCK_SIZE == 2 )
#  include "matrix22.h"
typedef Matrix22 SubBlock;
#elif ( BLOCK_SIZE == 3 )
#  include "matrix34.h"
typedef Matrix34 SubBlock;
#elif ( BLOCK_SIZE == 4 )
#  include "matrix44.h"
typedef Matrix44 SubBlock;
#endif

/// type for indices
typedef unsigned index_t;

// Flag to enable AVX implementation
#ifdef __AVX__
#  define MATRIXSB_USES_AVX 0 //REAL_IS_DOUBLE
#else
#  define MATRIXSB_USES_AVX 0
#endif

/// Sparse Matrix with block elements
/**
 MatrixSparseBlock uses a sparse storage, with arrays of elements for each column.
 Each element is a full square block of size DIM x DIM.
 
 Elements are stored in random order in the column.
 The lower triangle of the matrix is stored.
 
 FJN @ Cambridge, August-September 2019
 */
class MatrixSparseBlock
{
public:
    
    /// accessory class
    class Element;

private:
    
    /// A line of the sparse matrix
    class Line
    {
        friend class MatrixSparseBlock;

        size_t     size_;   ///< number of elements
        size_t     allo_;   ///< allocated size
        SubBlock * blk_;    ///< block elements
        SubBlock * sbk_;    ///< pointer for consolidate elements
        index_t  * inx_;    ///< column indices for each element
        
    public:
        
        /// constructor
        Line() { size_=0; allo_=0; inx_=nullptr; blk_=nullptr; sbk_=nullptr; }
        
        /// the assignment operator will transfer memory
        void operator =(Line&);
        
        /// destructor
        ~Line() { deallocate(); }
        
        /// allocate to hold 'nb' elements
        void allocate(size_t nb);
        
        /// deallocate memory
        void deallocate();

        /// set as zero
        void reset();
        
        /// sort element by increasing indices, using provided temporary array
        void sort(Element*&, size_t);
        
        /// print
        void print(std::ostream&) const;

        /// return n-th block (not necessarily, located at line inx_[n]
        SubBlock& operator[](int n) { return blk_[n]; }

        /// return block located at column 'j'
        SubBlock& block(index_t j);
        
        /// multiplication of a vector: L * X
        Vector vecMul(const real* X) const;
        
        /// multiplication of a vector: L * X
        real vecMul1D(const real* X) const;

#if MATRIXSB_USES_AVX
        /// multiplication of a vector: L * X
        vec2 vecMul2D(const real* X) const;
        
        /// multiplication of a vector: L * X
        vec2 vecMul2DU(const real* X) const;
        
        /// multiplication of a vector: L * X
        vec4 vecMul3DU(const real* X) const;
        
        /// multiplication of a vector: L * X
        vec4 vecMul3DU4(const real* X) const;

        /// multiplication of a vector: L * X
        vec4 vecMul3D(const real* X) const;
        
        /// multiplication of a vector: L * X
        vec4 vecMul4D(const real* X) const;
#endif
    };
    
    /// create Elements
    static size_t newElements(Element*& ptr, size_t);
    
    /// sort matrix block in increasing index order
    void sortElements();
    
    /// reallocate to use contiguous memory
    void consolidate();
    
    /// copy lower triangle into upper side
    void symmetrize();
    
    /// this is 'true' if symmetrize() was called
    bool already_symmetric;
    
private:
    
    /// number of lines in the matrix
    index_t   size_;
    
    /// amount of memory which has been allocated
    size_t    allocated_;
    
    /// array col_[c][] holds Elements of column 'c'
    Line *    row_;
    
    /// next_[ii] is the index of the first non-empty column of index >= ii
    index_t * next_;

    /// memory for consolidated version
    SubBlock* blocks_;
    
public:
    
    /// return the size of the matrix
    index_t size() const { return size_; }
    
    /// change the size of the matrix
    void resize(index_t s) { allocate(s); size_=s; }

    /// base for destructor
    void deallocate();
    
    /// default constructor
    MatrixSparseBlock();
    
    /// default destructor
    virtual ~MatrixSparseBlock()  { deallocate(); }
    
    /// set all the element to zero
    void reset();
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(size_t alc);
    
    /// returns element stored at line ii and column jj, if ( ii > jj )
    SubBlock& block(const index_t ii, const index_t jj)
    {
        assert_true( ii >= jj );
        assert_true( ii < size_ );
        assert_true( jj < size_ );
        assert_true( ii % BLOCK_SIZE == 0 );
        assert_true( jj % BLOCK_SIZE == 0 );
        return row_[ii].block(jj);
    }
    
    /// returns element stored at line ii and column jj, if ( ii > jj )
    SubBlock& diag_block(const index_t ii)
    {
        assert_true( ii < size_ );
        assert_true( ii % BLOCK_SIZE == 0 );
        return row_[ii].block(ii);
    }

    /// returns the address of element at (x, y), no allocation is done
    real* addr(index_t x, index_t y) const;
    
    /// returns the address of element at (x, y), allocating if necessary
    real& operator()(index_t x, index_t y);
    
    /// scale the matrix by a scalar factor
    void scale(real);
    
    /// add the diagonal block ( start, start+nb ) from this matrix to M
    void addDiagonalBlock(real* mat, unsigned ldd, index_t start, unsigned nb) const;
    
    /// add upper triangular half of 'this' block ( start, start+nb ) to `mat`
    void addTriangularBlock(real* mat, index_t ldd, index_t start, unsigned nb, unsigned dim) const;
    
    
    ///optional optimization that may accelerate multiplications by a vector
    void prepareForMultiply(int dim);

    /// multiplication of a vector, for columns within [start, stop[
    void vecMulAdd(const real*, real* Y, index_t start, index_t stop) const;
    
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd_SCAL(const real* X, real* Y, index_t start, index_t stop) const;

    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd2D(const real* X, real* Y, index_t start, index_t stop) const;

    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd3D(const real* X, real* Y, index_t start, index_t stop) const;
    
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd_ALT(const real* X, real* Y, index_t start, index_t stop) const;
    
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd_TIME(const real* X, real* Y, index_t start, index_t stop) const;

    
    /// multiplication of a vector, for columns within [start, stop[
    void vecMul(const real*, real* Y, index_t start, index_t stop) const;

    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMul2D(const real* X, real* Y, index_t start, index_t stop) const;

    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMul3D(const real* X, real* Y, index_t start, index_t stop) const;

    
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd(const real* X, real* Y) const { vecMulAdd(X, Y, 0, size_); }
    
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd_ALT(const real* X, real* Y) const { vecMulAdd_ALT(X, Y, 0, size_); }

    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd_SCAL(const real* X, real* Y) const { vecMulAdd_SCAL(X, Y, 0, size_); }

    /// 2D isotropic multiplication (not implemented)
    void vecMulAddIso2D(const real* X, real* Y) const {};
    
    /// 3D isotropic multiplication (not implemented)
    void vecMulAddIso3D(const real*, real*) const {};
    
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M)
    void vecMul(const real* X, real* Y) const { vecMul(X, Y, 0, size_); }

    /// true if matrix is non-zero
    bool nonZero() const;
    
    /// number of blocks which are not null
    size_t nbElements() const;

    /// returns a string which a description of the type of matrix
    std::string what() const;
    
    /// printf debug function in sparse mode: i, j : value
    void printSparse(std::ostream&) const;

    /// print content of one column
    void printLines(std::ostream&);
    
    /// debug function
    int bad() const;
};


#endif

