// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#ifndef SPARMATSYM1_H
#define SPARMATSYM1_H

#include "real.h"
#include <cstdio>
#include <string>

#define SPARMAT1_OPTIMIZED_MULTIPLY 1
#define SPARMAT1_USES_COLNEXT 1

///real symmetric sparse Matrix, with optimized multiplication
/**
 SparMatSym1 is similar to SparMatSym2 and uses the same sparse storage scheme,
 with independent arrays of elements for each column with sorted elements.
 Only the lower triangle of the matrix is stored.

 For multiplication, it uses a another format, from Numerical Recipes.
 The conversion is done by prepareForMultiply()
*/
class SparMatSym1 final
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

    /// diagonal elements
    real * diagon_;
    
#if SPARMAT1_USES_COLNEXT
    /// colidx_[i] is the index of the first non-empty column of index >= i
    size_t * colidx_;
#endif

    /// allocate column to hold specified number of values
    void allocateColumn(size_t jj, unsigned nb);
    
    /// allocate column to hold specified number of values
    void deallocateColumn(size_t jj);

    /// update colidx_[], a pointer to the next non-empty column
    void setColumnIndex();
    
#if SPARMAT1_OPTIMIZED_MULTIPLY

    /// allocated size of compact sparse storage arrays
    size_t nmax_;
    /// compact sparse storage indices
    unsigned * ija_;
    /// compact sparse storage values
    real * sa_;

#endif
    
    /// One column multiplication of a vector
    static void vecMulAddCol(const real* X, real* Y, real dia, size_t jj, Element col[], size_t cnt);
    
    /// One column multiplication of a vector, isotropic 2D version
    static void vecMulAddColIso2D(const real* X, real* Y, real dia, size_t jj, Element col[], size_t cnt);
    
    /// One column multiplication of a vector, isotropic 3D version
    static void vecMulAddColIso3D(const real* X, real* Y, real dia, size_t jj, Element col[], size_t cnt);


    /// One column multiplication of a vector
    void vecMulAddCol(const real* X, real* Y, size_t jj, real dia, size_t start, size_t stop) const;
    
    /// One column 2D isotropic multiplication of a vector
    void vecMulAddColIso2D(const real* X, real* Y, size_t jj, real dia, size_t start, size_t stop) const;
    
    /// One column 3D isotropic multiplication of a vector
    void vecMulAddColIso3D(const real* X, real* Y, size_t jj, real dia, size_t start, size_t stop) const;

    
    /// One column 2D isotropic multiplication of a vector
    void vecMulAddColIso2D_SSE(const double* X, double* Y, size_t jj, const double* dia, size_t start, size_t stop) const;
    
    /// One column 2D isotropic multiplication of a vector
    void vecMulAddColIso2D_SSEU(const double* X, double* Y, size_t jj, const double* dia, size_t start, size_t stop) const;

    
    /// One column 2D isotropic multiplication of a vector
    void vecMulAddColIso2D_AVX(const double* X, double* Y, size_t jj, const double* dia, size_t start, size_t stop) const;
    
    /// One column 2D isotropic multiplication of a vector
    void vecMulAddColIso2D_AVXU(const double* X, double* Y, size_t jj, const double* dia, size_t start, size_t stop) const;
    
    /// One column 3D isotropic multiplication of a vector
    void vecMulAddColIso3D_AVX(const double* X, double* Y, size_t jj, const double* dia, size_t start, size_t stop) const;

public:
    
    /// default constructor
    SparMatSym1();
    
    /// default destructor
    ~SparMatSym1()  { deallocate(); }

    /// return the size of the matrix
    size_t size() const { return size_; }
    
    /// change the size of the matrix
    void resize(size_t s) { allocate(s); size_=s; }
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(size_t sz);
    
    /// return total allocated memory
    size_t allocated() const;

    /// release memory
    void deallocate();
    
    /// set to zero
    void reset();
    
    /// return column at index j
    Element const* column(size_t j) const { return column_[j]; }
    
    /// number of elements in j-th column
    size_t column_size(size_t j) const { return colsiz_[j]; }

    /// returns the address of element at (x, y), no allocation is done
    real* addr(size_t x, size_t y) const;

    /// returns a modifiable diagonal element
    real& diagonal(size_t i) { return diagon_[i]; }
    
    /// returns the element at (i, j), allocating if necessary
    real& element(size_t i, size_t j);
    
    /// returns the element at (i, j), allocating if necessary
    real& operator()(size_t i, size_t j)
    {
        if ( i == j )
            return diagon_[i];
        return element(std::max(i, j), std::min(i, j));
    }

    /// scale the matrix by a scalar factor
    void scale(real);
    
    /// add terms with `i` and `j` in [start, start+cnt[ to `mat`
    void addDiagonalBlock(real* mat, size_t ldd, size_t start, size_t cnt, size_t amp=1) const;
    
    /// add scaled terms with `i` in [start, start+cnt[ if ( j > i ) and ( j <= i + rank ) to `mat`
    void addLowerBand(real alpha, real* mat, size_t ldd, size_t start, size_t cnt, size_t rank) const;
    
    /// add `alpha*trace()` for blocks within [start, start+cnt[ if ( j <= i + rank ) to `mat`
    void addDiagonalTrace(real alpha, real* mat, size_t ldd, size_t start, size_t cnt, size_t rank, bool sym) const;

    /// create compressed storage from column-based data
    bool prepareForMultiply(int);


    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M)
    void vecMulAdd(const real* X, real* Y, size_t start, size_t stop) const;
    
    /// 2D isotropic multiplication of a vector: Y <- Y + M * X with dim(X) = 2 * dim(M)
    void vecMulAddIso2D(const real* X, real* Y, size_t start, size_t stop) const;
    
    /// 3D isotropic multiplication of a vector: Y <- Y + M * X with dim(X) = 3 * dim(M)
    void vecMulAddIso3D(const real* X, real* Y, size_t start, size_t stop) const;

    /// multiplication of a vector, for columns within [start, stop[
    void vecMul(const real* X, real* Y, size_t start, size_t stop) const;

    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M)
    void vecMulAdd(const real* X, real* Y)      const { vecMulAdd(X, Y, 0, size_); }
    
    /// 2D isotropic multiplication of a vector: Y <- Y + M * X with dim(X) = 2 * dim(M)
    void vecMulAddIso2D(const real* X, real* Y) const { vecMulAddIso2D(X, Y, 0, size_); }
    
    /// 3D isotropic multiplication of a vector: Y <- Y + M * X with dim(X) = 3 * dim(M)
    void vecMulAddIso3D(const real* X, real* Y) const { vecMulAddIso3D(X, Y, 0, size_); }

    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M)
    void vecMulAdd_ALT(const real* X, real* Y)  const { vecMulAdd(X, Y, 0, size_); }
    
    /// multiplication of a vector: Y <- M * X with dim(X) = dim(M)
    void vecMul(const real* X, real* Y)         const { vecMul(X, Y, 0, size_); }

    /// true if matrix is non-zero
    bool isNotZero() const;
    
    /// number of elements in columns [start, stop[
    size_t nbElements(size_t start, size_t stop, size_t& alc) const;
    
    /// number of diagonal elements in columns [start, stop[
    size_t nbDiagonalElements(size_t start, size_t stop) const;

    /// number of elements in matrix
    size_t nbElements() const { size_t alc; return nbElements(0, size_, alc); }

    /// returns a string which a description of the type of matrix
    std::string what() const;
    
    /// print matrix columns in sparse mode: ( i, j : value ) if |value| >= inf
    void printSparse(std::ostream&, real inf, size_t start, size_t stop) const;
    
    /// print matrix in sparse mode: ( i, j : value ) if |value| >= inf
    void printSparse(std::ostream& os, real inf) const { printSparse(os, inf, 0, size_); }

    /// print content of one column
    void printColumn(std::ostream&, size_t);
    
    /// print content of one column
    void printSummary(std::ostream&, size_t start, size_t stop);

    /// printf debug function in sparse mode: i, j : value
    void printSparseArray(std::ostream&) const;
    
    /// debug function
    int bad() const;
};


#endif

