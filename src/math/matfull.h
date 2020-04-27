// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MATFULL_H
#define MATFULL_H

#include "real.h"
#include "assert_macro.h"
#include <iostream>
#include <string>

/// A non-symmetric square real Matrix
/**
 
 This matrix uses line-major storage and size is padded to a multiple of 4
 to enable optimal use of SIMD operations for matrix-vector multiplication.
 
 */
class MatrixFull
{    
private:
    
    /// size of matrix
    size_t     size_;

    /// number of block on a line
    size_t     nblk_;
    
    /// size of memory which has been allocated
    size_t     allo_;
    
    /// array of pointers to the blocks
    real*      mat_;
    
    /// index of block
    size_t block(size_t i, size_t j) const
    {
        assert_true( i < size_ );
        assert_true( j < size_ );
        return ( j  >> 2 ) + nblk_ * ( i >> 2 );
    }

public:
    
    /// return the size of the matrix
    size_t size() const { return size_; }
    
    /// change the size of the matrix
    void resize(size_t s) { allocate(s); size_=s; nblk_= ~3&(s+3); }
    
    /// default constructor
    MatrixFull();
    
    /// the deallocation
    void deallocate();
    
    /// allocate the matrix to be able to hold nb_block (arg 1) blocks
    void allocate(size_t nb_block);
    
    /// default destructor
    virtual ~MatrixFull() { deallocate(); }
    
    /// the address holding element (i, j)
    real* addr(size_t i, size_t j) const;
    
    /// returns the address of element at (i, j), allocating if necessary
    real value(size_t i, size_t j) const { return *addr(i, j); }

    /// returns the address of element at (i, j), allocating if necessary
    real& operator()(size_t i, size_t j) { return *addr(i, j); }

    
    /// reset values to zero
    void reset();
    
    /// scale all values
    void scale(real a);
    
    /// transpose
    void transpose();

    /// total number of elements allocated
    size_t nbElements() const { return size_ * size_; }

    /// vector multiplication: Y <- M * X
    void vecMulAdd(const real* X, real* Y) const;
    
    /// vector multiplication: Y <- M * X
    void vecMul(const real* X, real* Y) const;
    
    /// vector multiplication: Y <- M * X
    void vecMulAVX(const real* X, real* Y) const;

/*
    /// isotropic vector multiplication: Y = Y + M * X, size(X) = size(Y) = 2 * size(M)
    void vecMulAddIso2D(const real* X, real* Y) const { }
    
    /// isotropic vector multiplication: Y = Y + M * X, size(X) = size(Y) = 3 * size(M)
    void vecMulAddIso3D(const real* X, real* Y) const { }
*/
    /// maximum of the absolute value of all elements
    real norm_inf() const;
    
    /// output
    void print(std::ostream&) const;
    
    /// returns a string which a description of the type of matrix
    std::string what() const;
};


/// output operator
inline std::ostream& operator << (std::ostream& os, MatrixFull const& arg)
{
    std::ios::fmtflags fgs = os.flags();
    arg.print(os);
    os.setf(fgs);
    return os;
}

#endif
