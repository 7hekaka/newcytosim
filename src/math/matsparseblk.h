// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MATSPARSEBLK_H
#define MATSPARSEBLK_H

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

#define BLOCK_SIZE DIM

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


// Flag to enable AVX implementation
#ifdef __AVX__
#  define MATRIXSB_USES_AVX REAL_IS_DOUBLE
#else
#  define MATRIXSB_USES_AVX 0
#endif

/// Sparse Matrix with block elements
/**
 The lower triangle of the matrix is stored.
 Elements are stored in no particular order in each column.

 MatrixSparseBlock uses a sparse storage, with arrays of elements for each column.
 Each element is a full square block of size DIM x DIM.
 
 FJN @ Cambridge, August-September 2019
 */
class MatrixSparseBlock final
{
public:
    
    /// accessory class
    class Element;

    /// number of real in a block
    static constexpr size_t SB = sizeof(SubBlock) / sizeof(real);

private:
    
    /// A line of the sparse matrix
    class Line
    {
        friend class MatrixSparseBlock;

        size_t    size_;    ///< number of elements
        size_t    allo_;    ///< allocated size
        SubBlock * blk_;    ///< block elements
        SubBlock * sbk_;    ///< pointer for consolidate elements
        size_t   * inx_;    ///< column indices for each element
        
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
        
        /// sort element by increasing indices, using given temporary array
        void sortElements(Element[], size_t);
        
        /// print
        void print(std::ostream&) const;

        /// return n-th block (not necessarily, located at line inx_[n]
        SubBlock& operator[](size_t n) { return blk_[n]; }

        /// return block located at column 'j'
        SubBlock& block(size_t j);
        
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
    size_t    size_;
    
    /// amount of memory which has been allocated
    size_t    alloc_;
    
    /// array col_[c][] holds Elements of column 'c'
    Line *    row_;
    
    /// next_[ii] is the index of the first non-empty column of index >= ii
    size_t *  next_;

    /// memory for consolidated version
    SubBlock* blocks_;
    
public:
    
    /// return the size of the matrix
    size_t size() const { return size_; }
    
    /// change the size of the matrix
    void resize(size_t s) { allocate(s); size_=s; }

    /// base for destructor
    void deallocate();
    
    /// default constructor
    MatrixSparseBlock();
    
    /// default destructor
    ~MatrixSparseBlock()  { deallocate(); }
    
    /// set to zero
    void reset();
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(size_t alc);
    
    /// returns element stored at line ii and column jj, if ( ii > jj )
    SubBlock& block(const size_t ii, const size_t jj)
    {
        assert_true( ii >= jj );
        assert_true( ii < size_ );
        assert_true( jj < size_ );
        assert_true( ii % BLOCK_SIZE == 0 );
        assert_true( jj % BLOCK_SIZE == 0 );
        return row_[ii].block(jj);
    }
    
    /// returns element stored at line ii and column jj, if ( ii > jj )
    SubBlock& diag_block(const size_t ii)
    {
        assert_true( ii < size_ );
        assert_true( ii % BLOCK_SIZE == 0 );
        return row_[ii].block(ii);
    }

    /// returns the address of element at (x, y), no allocation is done
    real* addr(size_t x, size_t y) const;
    
    /// returns the address of element at (x, y), allocating if necessary
    real& operator()(size_t x, size_t y);
    
    /// scale the matrix by a scalar factor
    void scale(real);
    
    /// add terms within ( start, start+nb ) to `mat`
    void addDiagonalBlock(real* mat, size_t ldd, size_t start, size_t nb) const;
    
    /// add `alpha*trace()` for sub blocks within ( start, start+nb ) to `mat`
    void addDiagonalTrace(real alpha, real* mat, size_t ldd, size_t start, size_t nb) const;

    /// add lower triangle within ( start, start+nb ) to `mat`
    void addTriangularBlock(real* mat, size_t ldd, size_t start, size_t nb) const;
    
    
    /// prepare matrix for multiplications by a vector (must be called)
    bool prepareForMultiply(int);

    /// multiplication of a vector, for columns within [start, stop[
    void vecMulAdd(const real*, real* Y, size_t start, size_t stop) const;
    
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd_SCAL(const real* X, real* Y, size_t start, size_t stop) const;

    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd2D(const real* X, real* Y, size_t start, size_t stop) const;

    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd3D(const real* X, real* Y, size_t start, size_t stop) const;
    
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd_ALT(const real* X, real* Y, size_t start, size_t stop) const;
    
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMulAdd_TIME(const real* X, real* Y, size_t start, size_t stop) const;

    
    /// multiplication of a vector, for columns within [start, stop[
    void vecMul(const real*, real* Y, size_t start, size_t stop) const;

    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMul2D(const real* X, real* Y, size_t start, size_t stop) const;

    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(Y) = dim(M)
    void vecMul3D(const real* X, real* Y, size_t start, size_t stop) const;

    
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
    bool isNotZero() const;
    
    /// number of blocks which are not null
    size_t nbElements() const;

    /// returns a string which a description of the type of matrix
    std::string what() const;
    
    /// printf debug function in sparse mode: i, j : value
    void printSparse(std::ostream&, real) const;

    /// print content of one column
    void printLines(std::ostream&);
    
    /// debug function
    int bad() const;
};


#endif

