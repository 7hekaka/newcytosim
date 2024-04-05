// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#ifndef MATSYM_H
#define MATSYM_H

#include "real.h"
#include <cstdio>
#include <string>

/// A real symmetric Matrix
/**
the full lower triangular is stored
 (not used in Cytosim)
 */
class MatrixSymmetric final
{
private:
    
    /// leading dimension of array
    size_t ldim_;
    
    /// size of matrix
    size_t size_;

    /// size of memory which has been allocated
    size_t allocated_;
    
    // full upper triangle:
    real* sym_;
    
    // if 'false', the destructor will not delete the memory;
    bool in_charge_;
    
public:
    
    /// return the size of the matrix
    size_t size() const { return size_; }
    
    /// change the size of the matrix
    void resize(size_t s) { allocate(s); size_=s; }

    /// base for destructor
    void deallocate();
    
    /// default constructor
    MatrixSymmetric();
    
    
    /// constructor from an existing array
    MatrixSymmetric(size_t s)
    {
        sym_ = nullptr;
        resize(s);
        ldim_ = s;
        sym_ = new_real(s*s);
        zero_real(s*s, sym_);
        in_charge_ = true;
    }

    /// constructor from an existing array
    MatrixSymmetric(size_t s, real* array, size_t ldd)
    {
        free_real(sym_);
        size_ = s;
        ldim_ = ldd;
        sym_ = array;
        in_charge_ = false;
    }
    
    /// default destructor
    ~MatrixSymmetric()  { deallocate(); }
    
    /// set to zero
    void reset();
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(size_t alc);
    
    /// returns address of data array
    real* data() const { return sym_; }

    /// returns the address of element at (x, y), no allocation is done
    real* address(size_t x, size_t y) const;
    
    /// returns the address of element at (x, y), allocating if necessary
    real& operator()(size_t i, size_t j);
    
    /// scale the matrix by a scalar factor
    void scale(real a);
    
    /// multiplication of a vector: Y = Y + M * X, dim(X) = dim(M)
    void vecMulAdd(const real* X, real* Y) const;
    
    /// 2D isotropic multiplication of a vector: Y = Y + M * X
    void vecMulAddIso2D(const real* X, real* Y) const;
    
    /// 3D isotropic multiplication of a vector: Y = Y + M * X
    void vecMulAddIso3D(const real* X, real* Y) const;
    
    /// true if matrix is non-zero
    bool notZero() const;
    
    /// number of element which are non-zero
    size_t nbElements(size_t start, size_t stop) const;
    
    /// returns a string which a description of the type of matrix
    std::string what() const;
};

#endif
