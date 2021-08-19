// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#ifndef SPARMATSYMBLKDIAG_H
#define SPARMATSYMBLKDIAG_H

#include "dim.h"
#include "real.h"
#include <cstdio>
#include <iostream>
#include "assert_macro.h"


/**
 The block size 'BLOCK_SIZE' can be defined on the command line during compilation,
 and is otherwise set here, to match the dimensionality of the simulation
 */

#define BLOCK_SIZE ( DIM < 3 ? DIM : 3 )

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

 MatrixSparseSymmetricBlock uses a sparse storage, with arrays of elements for each column.
 Each element is a full square block of size DIM x DIM.
 
 F. Nedelec, 17--27 March 2017, revised entirely June 2018, Nov 2020
 */
class SparMatSymBlkDiag final
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
    
    /// A column of the sparse matrix
    class alignas(32) Column
    {
        friend class SparMatSymBlkDiag;
        friend class Meca;

        Block    dia_;   ///< diagonal block
        size_t * inx_;   ///< line index for each element
        Block  * blk_;   ///< off-diagonal blocks
        size_t  allo_;   ///< allocated size of array
        size_t  size_;   ///< number of blocks in column

    public:
        
        /// constructor
        Column();
        
        /// the assignment operator will transfer memory ownership
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
        bool isNotZero() const { return ( dia_ != 0.0 ) || ( size_ > 0 ); }
        
        /// return block located on the diagonal
        Block& diag_block() { return dia_; }

        /// return block located at line 'i' and column 'j'
        Block& block(size_t i, size_t j);
        
        /// return n-th off-diagonal block, located at line inx_[n]
        Block& operator[](size_t n) const { return blk_[n]; }

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
        void vecMulAddTriangle3D_SSE(const float* X, float* Y, size_t j) const;

        
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 3
        void vecMulAdd3D_AVX(const double* X, double* Y, size_t j) const;
        
        /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M), block_size = 3
        void vecMulAddTriangle3D_AVX(const double* X, double* Y, size_t j) const;

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
    Column * pilar_;
    
    /// colix_[i] is the index of the first non-empty column of index >= i
    size_t * colix_;

public:
    
    /// return the size of the matrix
    size_t size() const { return size_; }
    
    /// change the size of the matrix
    void resize(size_t s) { allocate(s); size_=s; }

    /// base for destructor
    void deallocate();
    
    /// default constructor
    SparMatSymBlkDiag();
    
    /// default destructor
    ~SparMatSymBlkDiag()  { deallocate(); }
    
    /// set to zero
    void reset();
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(size_t alc);
    
    /// return column at index j
    Column const& column(size_t j) const { return pilar_[j/BLOCK_SIZE]; }
    
    /// return column at index j
    Column&       column(size_t j)       { return pilar_[j/BLOCK_SIZE]; }

    /// returns element stored at line ii and column jj, if ( ii > jj )
    Block& block(const size_t ii, const size_t jj)
    {
        assert_true( ii < size_ );
        assert_true( jj < size_ );
        assert_true( ii % BLOCK_SIZE == 0 );
        assert_true( jj % BLOCK_SIZE == 0 );
#if ( 0 )
        // safe swap, with branchless code:
        size_t i = std::max(ii, jj);
        size_t j = std::min(ii, jj);
        return column(j).block(i, j);
#else
        assert_true( ii > jj );
        return column(jj).block(ii, jj);
#endif
    }
    
    /// returns element at (i, i)
    Block& diag_block(size_t i) { return column(i).diag_block(); }
    
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
    void vecMulDiagonal3D(const real* X, real* Y) const;

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
    void printColumns(std::ostream&, size_t start, size_t stop);
    
    /// debug function
    int bad() const;
};


#endif

