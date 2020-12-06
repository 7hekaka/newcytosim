// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#ifndef SPARMATSYM_H
#define SPARMATSYM_H

#include "real.h"
#include <cstdio>
#include <string>

///real symmetric sparse Matrix
/**
 SparMatSym uses a sparse storage, with arrays of elements for each column.
 */
class SparMatSym final
{
public:
    
    /// An element of the sparse matrix
    struct Element
    {
        real     val;  ///< The value of the element
        unsigned inx;  ///< The index of the line
        
        void reset(size_t i)
        {
            inx = static_cast<unsigned>(i);
            val = 0;
        }
    };
    
private:
    
    /// size of matrix
    size_t size_;
    
    /// amount of memory allocated
    size_t alloc_;
    
    /// array column_[c][] holds Elements of column 'c'
    Element ** column_;
    
    /// colsiz_[c] is the number of Elements in column 'c'
    unsigned * colsiz_;
    
    /// colmax_[c] is the number of Elements allocated in column 'c'
    unsigned * colmax_;
    
    /// allocate column to hold specified number of values
    void allocateColumn(size_t col, size_t nb);
    
public:
    
    /// return the size of the matrix
    size_t size() const { return size_; }
    
    /// change the size of the matrix
    void resize(size_t s) { allocate(s); size_=s; }

    /// base for destructor
    void deallocate();
    
    /// default constructor
    SparMatSym();
    
    /// default destructor
    ~SparMatSym()  { deallocate(); }
    
    /// set to zero
    void reset();
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(size_t sz);
    
    /// returns the address of element at (x, y), no allocation is done
    real* addr(size_t x, size_t y) const;
    
    /// returns a modifiable reference to the diagonal term at given index
    real& diagonal(size_t ix);
    
    /// returns the address of element at (x, y), allocating if necessary
    real& operator()(size_t x, size_t y);
    
    /// scale the matrix by a scalar factor
    void scale(real);
    
    /// add lower triangular half of 'this' block ( idx, idx, idx+siz, idx+siz ) to `mat`
    void addDiagonalBlock(real* mat, size_t ldd, size_t start, size_t cnt, size_t amp=1) const;
    
    /// prepare matrix for multiplications by a vector (must be called)
    void prepareForMultiply(int dim);
    
    /// multiplication of a vector: Y = Y + M * X with dim(X) = dim(M)
    void vecMulAdd(const real* X, real* Y) const;
    
    /// multiplication of a vector: Y = Y + M * X with dim(X) = dim(M)
    void vecMulAdd_ALT(const real* X, real* Y) const { vecMulAdd(X, Y); }

    /// 2D isotropic multiplication of a vector: Y = Y + M * X with dim(X) = 2 * dim(M)
    void vecMulAddIso2D(const real* X, real* Y) const;
    
    /// 3D isotropic multiplication of a vector: Y = Y + M * X with dim(X) = 3 * dim(M)
    void vecMulAddIso3D(const real* X, real* Y) const;
    
    /// true if matrix is non-zero
    bool isNotZero() const;
    
    /// number of elements in columns [start, stop[
    size_t nbElements(size_t start, size_t stop) const;
    
    /// number of elements in matrix
    size_t nbElements() const { return nbElements(0, size_); }

    /// returns a string which a description of the type of matrix
    std::string what() const;
    
    /// print matrix columns in sparse mode: ( i, j : value ) if |value| >= inf
    void printSparse(std::ostream&, real inf, size_t start, size_t stop) const;
    
    /// print matrix in sparse mode: ( i, j : value ) if |value| >= inf
    void printSparse(std::ostream& os, real inf) const { printSparse(os, inf, 0, size_); }

    /// print content of one column
    void printColumn(std::ostream&, size_t);
    
    /// print content of one column
    void printColumns(std::ostream&, size_t start, size_t stop);

    /// debug function
    int bad() const;
};


#endif

