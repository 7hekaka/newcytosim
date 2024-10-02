// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#ifndef SPARMAT_H
#define SPARMAT_H

#include "real.h"
#include <cstdio>
#include <string>

/// a real (non-symmetric) sparse Matrix
/**
 This class is not used currently in Cytosim
 */
class SparMat final
{
private:
    
    /// size of matrix
    index_t size_;

    /// size of memory which has been allocated
    index_t alloc_;

    // array [ size ][ ? ] holding the values for each column
    real ** mxCol;
    
    // array [ size ][ ? ] holding the line index for each column
    index_t ** mxRow;
    
    // allocate column to hold nb values
    void allocateColumn(index_t column_index, index_t num_values);
    
public:
    
    /// return the size of the matrix
    index_t size() const { return size_; }
    
    /// change the size of the matrix
    void resize(index_t s) { allocate(s); size_=s; }

    /// base for destructor
    void deallocate();
    
    /// default constructor
    SparMat();
    
    /// default destructor
    ~SparMat()  { deallocate(); }
    
    /// set to zero
    void reset();
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(index_t sz);
        
    /// returns the address of element at (x, y), no allocation is done
    real* address(index_t x, index_t y ) const;
    
    /// returns the address of element at (x, y), allocating if necessary
    real& operator()(index_t x, index_t y );
    
    /// scale the matrix by a scalar factor
    void scale( real a );
    
    /// add terms with `i` and `j` in [start, start+cnt[ to `mat`
    void addDiagonalBlock(real* mat, index_t ldd, index_t start, index_t cnt, index_t mul) const;
    
    /// multiplication of a vector: Y = Y + M * X, dim(X) = dim(M)
    void vecMulAdd(const real* X, real* Y) const;
    
    /// 2D isotropic multiplication of a vector: Y = Y + M * X
    void vecMulAddIso2D(const real* X, real* Y) const;
    
    /// 3D isotropic multiplication of a vector: Y = Y + M * X
    void vecMulAddIso3D(const real* X, real* Y) const;
    
    /// true if matrix is non-zero
    bool notZero() const;
    
    /// number of element which are non-zero
    size_t nbElements(index_t start, index_t stop) const;
    
    /// number of blocks which are not null
    size_t nbElements() const { return nbElements(0, size_); }

    /// returns a string which a description of the type of matrix
    std::string what() const;
    
    /// print matrix columns in sparse mode: ( i, j : value ) if |value| >= inf
    void printSparse(std::ostream&, real inf, index_t start, index_t stop) const;
    
    /// print matrix in sparse mode: ( i, j : value ) if |value| >= inf
    void printSparse(std::ostream& os, real inf) const { printSparse(os, inf, 0, size_); }

    /// debug function
    int bad() const;
};


#endif
