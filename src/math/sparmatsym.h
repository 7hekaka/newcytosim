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
        real    val;  ///< The value of the element
        index_t inx;  ///< The index of the line
        
        void reset(index_t i)
        {
            inx = i;
            val = 0;
        }
    };
    
private:
    
    /// size of matrix
    index_t size_;
    
    /// amount of memory allocated
    index_t alloc_;
    
    /// array column_[c][] holds Elements of column 'c'
    Element ** column_;
    
    /// colsiz_[c] is the number of Elements in column 'c'
    index_t * colsiz_;
    
    /// colmax_[c] is the number of Elements allocated in column 'c'
    index_t * colmax_;
    
    /// allocate column to hold specified number of values
    void allocateColumn(index_t col, index_t nb);
    
public:
    
    /// return the size of the matrix
    index_t size() const { return size_; }
    
    /// change the size of the matrix
    void resize(index_t s) { allocate(s); size_=s; }

    /// base for destructor
    void deallocate();
    
    /// default constructor
    SparMatSym();
    
    /// default destructor
    ~SparMatSym()  { deallocate(); }
    
    /// set to zero
    void reset();
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(index_t sz);
    
    /// returns the address of element at (x, y), no allocation is done
    real* address(index_t x, index_t y) const;
    
    /// returns a modifiable reference to the diagonal term at given index
    real& diagonal(index_t i);
    
    /// returns the element at (i, j), allocating if necessary
    real& element(index_t i, index_t j);
    
    /// returns the element at (i, j), allocating if necessary
    real& operator()(index_t i, index_t j)
    {
        return element(std::max(i, j), std::min(i, j));
    }

    /// scale the matrix by a scalar factor
    void scale(real);
    
    /// add terms with `i` and `j` in [start, start+cnt[ to `mat`
    void addDiagonalBlock(real* mat, index_t ldd, index_t start, index_t cnt, index_t mul) const;
    
    /// prepare matrix for multiplications by a vector (must be called)
    bool prepareForMultiply(int dim);
    
    /// multiplication of a vector: Y = Y + M * X with dim(X) = dim(M)
    void vecMulAdd(const real* X, real* Y) const;
    
    /// multiplication of a vector: Y = Y + M * X with dim(X) = dim(M)
    void vecMulAdd_ALT(const real* X, real* Y) const { vecMulAdd(X, Y); }

    /// 2D isotropic multiplication of a vector: Y = Y + M * X with dim(X) = 2 * dim(M)
    void vecMulAddIso2D(const real* X, real* Y) const;
    
    /// 3D isotropic multiplication of a vector: Y = Y + M * X with dim(X) = 3 * dim(M)
    void vecMulAddIso3D(const real* X, real* Y) const;
    
    /// true if matrix is non-zero
    bool notZero() const;
    
    /// number of elements in columns [start, stop[
    size_t nbElements(index_t start, index_t stop) const;
    
    /// number of elements in matrix
    size_t nbElements() const { return nbElements(0, size_); }

    /// returns a string which a description of the type of matrix
    std::string what() const;
    
    /// print matrix columns in sparse mode: ( i, j : value ) if |value| >= inf
    void printSparse(std::ostream&, real inf, index_t start, index_t stop) const;
    
    /// print matrix in sparse mode: ( i, j : value ) if |value| >= inf
    void printSparse(std::ostream& os, real inf) const { printSparse(os, inf, 0, size_); }

    /// print content of one column
    void printColumn(std::ostream&, index_t);
    
    /// print content of one column
    void printSummary(std::ostream&, index_t start, index_t stop);

    /// debug function
    int bad() const;
};


#endif

