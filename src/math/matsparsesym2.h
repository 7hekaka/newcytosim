// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MATSPARSESYM2_H
#define MATSPARSESYM2_H

#include "real.h"
#include <cstdio>
#include <string>

#define MATRIX2_OPTIMIZED_MULTIPLY 1

///real symmetric sparse Matrix, with optimized multiplication
/**
 MatrixSparseSymmetric2 uses a sparse storage, with arrays of elements for each column.
 For multiplication, it uses a another format, from Numerical Recipes.
 The conversion is done when prepareForMultiply() is called
 
 Elements are stored in random order in the column.
*/
class MatrixSparseSymmetric2
{
public:
    
    /// An element of the sparse matrix
    struct Element
    {
        real     val;   ///< The value of the element
        size_t  inx;   ///< The index of the line
        
        void reset(size_t i)
        {
            inx = i;
            val = 0.0;
        }
    };
    
private:
    
    /// size of matrix
    size_t  size_;

    /// amount of memory which has been allocated
    size_t  allocated_;
    
    /// array col_[c][] holds Elements of column 'c'
    Element ** col_;
    
    /// col_size_[c] is the number of Elements in column 'c'
    size_t   * col_size_;
    
    /// col_max_[c] is the number of Elements allocated in column 'c'
    size_t   * col_max_;
    
    /// allocate column to hold specified number of values
    void allocateColumn(size_t col, size_t nb);
   
#if MATRIX2_OPTIMIZED_MULTIPLY
    
    /// next_[ii] is the index of the first non-empty column of index >= ii
    size_t * next_;
    
    /// update next_[], a pointer to the next non-empty column
    void setNextColumn();

    ///array of index for the optmized multiplication
    ///@todo migrate to DSS Symmetric Matrix Storage format
    size_t     nmax_;
    size_t   * ija_;
    real     * sa_;

#endif

public:
    
    /// return the size of the matrix
    size_t size() const { return size_; }
    
    /// change the size of the matrix
    void resize(size_t s) { allocate(s); size_=s; }

    /// base for destructor
    void deallocate();
    
    /// default constructor
    MatrixSparseSymmetric2();
    
    /// default destructor
    virtual ~MatrixSparseSymmetric2()  { deallocate(); }
    
    /// set to zero
    void reset();
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(size_t sz);
    
    /// returns the address of element at (x, y), no allocation is done
    real* addr(size_t x, size_t y) const;
    
    /// returns the address of element at (x, y), allocating if necessary
    real& operator()(size_t x, size_t y);
    
    /// scale the matrix by a scalar factor
    void scale(real);
    
    /// add the diagonal block ( x, x, x+sx, x+sx ) from this matrix to M
    void addDiagonalBlock(real* mat, size_t ldd, size_t si, size_t nb) const;
    
    /// add upper triangular half of 'this' block ( idx, idx, idx+siz, idx+siz ) to `mat`
    void addTriangularBlock(real* mat, size_t ldd, size_t si, size_t nb, size_t dim) const;
    
    ///optional optimization that may accelerate multiplications by a vector
    void prepareForMultiply(int dim);
    
    /// multiplication of a vector: Y = Y + M * X with dim(X) = dim(M)
    void vecMulAdd(const real* X, real* Y, size_t start, size_t stop) const;

    /// 2D isotropic multiplication of a vector: Y = Y + M * X with dim(X) = 2 * dim(M)
    void vecMulAddIso2D(const real* X, real* Y, size_t start, size_t stop) const;
    
    /// 3D isotropic multiplication of a vector: Y = Y + M * X with dim(X) = 3 * dim(M)
    void vecMulAddIso3D(const real* X, real* Y, size_t start, size_t stop) const;
    
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M)
    void vecMulAdd(const real* X, real* Y)      const { vecMulAdd(X, Y, 0, size_); }
    
    /// 2D isotropic multiplication of a vector: Y <- Y + M * X with dim(X) = 2 * dim(M)
    void vecMulAddIso2D(const real* X, real* Y) const { vecMulAddIso2D(X, Y, 0, size_); }
    
    /// 3D isotropic multiplication of a vector: Y <- Y + M * X with dim(X) = 3 * dim(M)
    void vecMulAddIso3D(const real* X, real* Y) const { vecMulAddIso3D(X, Y, 0, size_); }

    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M)
    void vecMulAdd_ALT(const real* X, real* Y)  const { vecMulAdd(X, Y, 0, size_); }

    /// true if matrix is non-zero
    bool isNotZero() const;
    
    /// number of element which are not null
    size_t nbElements(size_t start, size_t stop) const;
    
    /// number of blocks which are not null
    size_t nbElements() const { return nbElements(0, size_); }

    /// returns a string which a description of the type of matrix
    std::string what() const;
    
    /// printf debug function in sparse mode: i, j : value
    void printSparse(std::ostream&) const;
    
    /// print content of one column
    void printColumn(std::ostream&, size_t);
    
    /// print content of one column
    void printColumns(std::ostream&);

    /// debug function
    int bad() const;
};


#endif

