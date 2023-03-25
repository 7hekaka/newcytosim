// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MATRIX_H
#define MATRIX_H

#include <string>
#include <iostream>
#include "real.h"


/// Matrix is a model for all the large matrices
/**
 This interface is a model for all the large matrices, that is not enforced
 by the compiler, as Matrices are their own class and do not derive from this.
 */
class Matrix
{
private:
    
    /// Disabled copy constructor (@todo: write copy constructor)
    Matrix(Matrix const&);
    
    /// Disabled copy assignment (@todo: write copy assignement)
    Matrix& operator = (Matrix const&);

protected:

    /// size of matrix
    size_t   size_;

public:
    
    /// constructor
    Matrix() { size_ = 0; }

    /// constructor
    Matrix(size_t s) { size_ = s; }
    
    /// return the size of the matrix
    size_t size() const { return size_; }
    
    /// change the size of the matrix
    void resize(size_t s) { allocate(s); size_=s; }

    //----------------------------------------------------------------------
    
    /// allocate the matrix to hold ( sz * sz ), all values may be lost
    virtual void allocate(size_t alc) = 0;
        
    /// returns the address of element at (x, y), no allocation is done
    virtual real* addr(size_t x, size_t y) const = 0;
    
    /// returns the address of element at (x, y), allocating if necessary
    virtual real& operator()(size_t x, size_t y) = 0;
    
    /// returns the value of element at (x, y) or zero if not allocated
    real value(size_t x, size_t y) const;
    
    //----------------------------------------------------------------------
    
    /// set all the elements to zero
    virtual void reset() = 0;
    
    /// scale the matrix by a scalar factor
    virtual void scale(real) = 0;
    
    /// copy the block ( x, y, x+sx, y+sy ) into `mat`
    void copyBlock(real* mat, size_t ldd, size_t sx, size_t nx, size_t sy, size_t ny) const;
    
    /// add terms with `i` and `j` in [start, start+cnt[ to `mat`
    virtual void addDiagonalBlock(real* mat, size_t ldd, size_t start, size_t cnt, size_t mul, size_t amp=1) const;
    
    //----------------------------------------------------------------------
    
    /// Optional optimization to accelerate multiplications below
    virtual bool prepareForMultiply(int dim) { return true; }
    
    /// Vector multiplication: Y <- Y + M * X, size(X) = size(Y) = size(M)
    virtual void vecMulAdd(const real* X, real* Y) const = 0;
    
    /// Vector multiplication: Y <- M * X, size(X) = size(Y) = size(M)
    virtual void vecMul(const real* X, real* Y) const;
    
    //----------------------------------------------------------------------
    
    /// maximum absolute value among all the elements
    virtual real norm_inf() const;
    
    /// true if the matrix is non-zero
    virtual bool notZero() const;
    
    /// number of element which are not null
    virtual size_t nbElements(size_t start, size_t stop) const;
    
    /// number of blocks which are not null
    size_t nbElements() const { return nbElements(0, size_); }

    /// returns a string which a description of the type of matrix
    virtual std::string what() const = 0;
    
    /// print matrix columns in sparse mode: ( i, j : value ) if |value| >= inf
    void printSparse(std::ostream&, real inf, size_t start, size_t stop) const;
    
    /// print matrix in sparse mode: ( i, j : value ) if |value| >= inf
    void printSparse(std::ostream& os, real inf) const { printSparse(os, inf, 0, size_); }
    
    /// printf debug function in full lines, all columns
    virtual void printFull(std::ostream&) const;
    
};


#endif
