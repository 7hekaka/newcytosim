// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

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
 
 FJN 27.04.2020
 */
class MatrixFull final
{
    /// number of real in a matrix block
    static constexpr index_t SB = 16U;

private:
    
    /// size of matrix
    index_t size_;

    /// number of block on a line
    index_t nblk_;
    
    /// size of memory which has been allocated
    index_t allo_;
    
    /// array of pointers to the blocks
    real*   mat_;
    
    /// index of block
    index_t block(index_t i, index_t j) const
    {
        assert_true( i < size_ );
        assert_true( j < size_ );
        return ( j >> 2 ) + nblk_ * ( i >> 2 );
    }

public:
    
    /// return the size of the matrix
    index_t size() const { return size_; }
    
    /// change the size of the matrix
    void resize(index_t s) { allocate(s); size_=s; nblk_= (~3U&(s+3U))>>2; }
    
    /// default constructor
    MatrixFull();
    
    /// the deallocation
    void deallocate();
    
    /// allocate the matrix to be able to hold `nb` blocks
    void allocate(index_t nb);
    
    /// default destructor
    ~MatrixFull() { deallocate(); }
    
    /// return memory used to store values
    real* data() const { return mat_; }
    
    /// the address holding element (i, j)
    real* address(index_t i, index_t j) const;
    
    /// returns the address of element at (i, j), allocating if necessary
    real value(index_t i, index_t j) const { return *address(i, j); }

    /// returns the address of element at (i, j), allocating if necessary
    real& operator()(index_t i, index_t j) { return *address(i, j); }
    
    /// returns the address of element at (i, j), allocating if necessary
    real  operator()(index_t i, index_t j) const { return *address(i, j); }

    
    /// reset with 'dia' on diagonal and 'off' elsewhere
    void reset(real off, real dia);
    
    /// reset terms 'kl' below the diagonal or 'ku' above
    void truncate(index_t kl, index_t ku);

    /// import column-major matrix
    void importMatrix(index_t size, real const*, index_t lld);
    
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
    void vecMul0(const real* X, real* Y) const;
    
    /// vector multiplication: Y <- transposed(M) * X
    void transVecMulAdd(const real* X, real* Y) const;
    
    /// vector multiplication: Y <- transposed(M) * X
    void transVecMul(const real* X, real* Y) const
    { zero_real(size_, Y); transVecMulAdd(X, Y); }

    /// maximum of the absolute value of all elements
    real norm_inf() const;
    
    /// output part of matrix
    void print(std::ostream&, index_t imin, index_t imax, index_t jmin, index_t jmax) const;

    /// output
    void print(std::ostream& os) const { print(os, 0, size_, 0, size_); }
    
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
