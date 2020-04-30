// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

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
class MatrixSymmetric
{
private:
    
    /// leading dimension of array
    size_t   msLDD;
    
    /// size of matrix
    size_t   size_;

    /// size of memory which has been allocated
    size_t   allocated_;
    
    // full upper triangle:
    real* val;
    
    // if 'false', destructor will not call delete[] val;
    bool in_charge;
    
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
        val = nullptr;
        resize(s);
        msLDD = s;
        val = new_real(s*s);
        zero_real(s*s, val);
        in_charge = true;
    }

    /// constructor from an existing array
    MatrixSymmetric(size_t s, real* array, size_t ldd)
    {
        free_real(val);
        size_ = s;
        msLDD = ldd;
        val = array;
        in_charge = false;
    }
    
    /// default destructor
    virtual ~MatrixSymmetric()  { deallocate(); }
    
    /// set to zero
    void reset();
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(size_t alc);
    
    /// returns address of data array
    real* data() const { return val; }

    /// returns the address of element at (x, y), no allocation is done
    real* addr(size_t x, size_t y) const;
    
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
    bool isNotZero() const;
    
    /// number of element which are non-zero
    size_t nbElements(size_t start, size_t stop) const;
    
    /// returns a string which a description of the type of matrix
    std::string what() const;
};

#endif
