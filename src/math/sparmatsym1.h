// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#ifndef SPARMATSYM1_H
#define SPARMATSYM1_H

#include "assert_macro.h"
#include "real.h"
#include <cstdio>
#include <string>

#define SPARMAT1_COMPACTED 0
#define SPARMAT1_USES_COLNEXT 1

///real symmetric sparse Matrix, with optimized multiplication
/**
 SparMatSym1 is equivalent to SparMatSym2 and uses a similar storage scheme.
 It stores diagonal elements in 'diagon_' and off-diagonal elements in sorted
 array, with independent array for each column.
 The matrix is symmetric and only the lower triangle of the matrix is stored.

 if SPARMAT1_COMPACTED==1, it uses a compact storage for multiplication.
 This incurrs additional memory and performance may not be improved.
 The conversion is done by prepareForMultiply()
*/
class SparMatSym1 final
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
    
    /// number of columns allocated
    index_t alloc_;
    
    /// array column_[c][] holds Elements of column 'c'
    Element ** column_;
    
    /// colsiz_[c] is the number of Elements in column 'c'
    index_t * colsiz_;
    
    /// colmax_[c] is the number of Elements allocated in column 'c'
    index_t * colmax_;

    /// diagonal elements
    real * diagon_;
    
#if SPARMAT1_USES_COLNEXT
    /// colidx_[i] is the index of the first non-empty column of index >= i
    index_t * colidx_;
#endif

    /// allocate column to hold specified number of values
    void allocateColumn(index_t jj, index_t nb);
    
    /// allocate column to hold specified number of values
    void deallocateColumn(index_t jj);

    /// update colidx_[], a pointer to the next non-empty column
    void setColumnIndex();
    
#if SPARMAT1_COMPACTED

    /// allocated size of compact sparse storage arrays
    index_t nmax_;
    /// index to ija_[] and elm_[] where each column starts
    index_t * off_;
    /// compact sparse storage indices
    index_t * ija_;
    /// compact sparse storage values
    real * elm_;

#endif
    
    /// One column multiplication of a vector
    static void vecMulAddCol(const real* X, real* Y, index_t jj, real dia, Element col[], index_t cnt);
    
    /// One column multiplication of a vector, isotropic 2D version
    static void vecMulAddColIso2D(const real* X, real* Y, index_t jj, real dia, Element col[], index_t cnt);
    
    /// One column multiplication of a vector, isotropic 3D version
    static void vecMulAddColIso3D(const real* X, real* Y, index_t jj, real dia, Element col[], index_t cnt);


    /// One column multiplication of a vector
    void vecMulAddCol(const real* X, real* Y, index_t jj, real dia, index_t start, index_t stop) const;
    
    /// One column 2D isotropic multiplication of a vector
    void vecMulAddColIso2D(const real* X, real* Y, index_t jj, real dia, index_t start, index_t stop) const;
    
    /// One column 3D isotropic multiplication of a vector
    void vecMulAddColIso3D(const real* X, real* Y, index_t jj, real dia, index_t start, index_t stop) const;

    
    /// One column 2D isotropic multiplication of a vector
    void vecMulAddColIso2D_SSE(const double* X, double* Y, index_t jj, const double* dia, index_t start, index_t stop) const;
    
    /// One column 2D isotropic multiplication of a vector
    void vecMulAddColIso2D_SSEU(const double* X, double* Y, index_t jj, const double* dia, index_t start, index_t stop) const;

    
    /// One column 2D isotropic multiplication of a vector
    void vecMulAddColIso2D_AVX(const double* X, double* Y, index_t jj, const double* dia, index_t start, index_t stop) const;
    
    /// One column 2D isotropic multiplication of a vector
    void vecMulAddColIso2D_AVXU(const double* X, double* Y, index_t jj, const double* dia, index_t start, index_t stop) const;
    
    /// One column 3D isotropic multiplication of a vector
    void vecMulAddColIso3D_AVX(const double* X, double* Y, index_t jj, const double* dia, index_t start, index_t stop) const;

public:
    
    /// default constructor
    SparMatSym1();
    
    /// default destructor
    ~SparMatSym1()  { deallocate(); }

    /// return the size of the matrix
    index_t size() const { return size_; }
    
    /// change the size of the matrix
    void resize(index_t s) { allocate(s); size_=s; }
    
    /// allocate the matrix to hold ( sz * sz )
    void allocate(index_t sz);
    
    /// return total allocated memory
    index_t allocated() const;

    /// release memory
    void deallocate();
    
    /// set to zero
    void reset();
    
    /// number of columns
    index_t num_columns() const { return size_; }
    
    /// return column at index j
    Element const* column(index_t j) const { assert_true(j<size_); return column_[j]; }
    
    /// number of elements in j-th column
    index_t column_size(index_t j) const { assert_true(j<size_); return colsiz_[j]; }
    
    /// line index of n-th element in j-th column
    index_t column_index(index_t j, index_t n) const { assert_true(j<size_); return column_[j][n].inx; }

    /// returns the address of element at (x, y), no allocation is done
    real* address(index_t x, index_t y) const;

    /// returns a modifiable diagonal element
    real& diagonal(index_t i) { return diagon_[i]; }
    
    /// returns the element at (i, j), allocating if necessary
    real& element(index_t i, index_t j);
    
    /// returns the element at (i, j), allocating if necessary
    real& operator()(index_t i, index_t j)
    {
        if ( i == j )
            return diagon_[i];
        return element(std::max(i, j), std::min(i, j));
    }

    /// scale the matrix by a scalar factor
    void scale(real);
    
    /// add terms with `i` and `j` in [start, start+cnt[ to `mat`
    void addDiagonalBlock(real* mat, index_t ldd, index_t start, index_t cnt, index_t mul) const;
    
    /// add scaled terms with `i` in [start, start+cnt[ if ( j > i ) and ( j <= i + rank ) to `mat`
    void addLowerBand(real alpha, real* mat, index_t ldd, index_t start, index_t cnt, index_t mul, index_t rank) const;
    
    /// add `alpha*trace()` for blocks within [start, start+cnt[ if ( j <= i + rank ) to `mat`
    void addDiagonalTrace(real alpha, real* mat, index_t ldd, index_t start, index_t cnt, index_t mul, index_t rank, bool sym) const;

    /// create compressed storage from column-based data
    bool prepareForMultiply(int);


    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M)
    void vecMulAdd(const real* X, real* Y, index_t start, index_t stop) const;
    
    /// 2D isotropic multiplication of a vector: Y <- Y + M * X with dim(X) = 2 * dim(M)
    void vecMulAddIso2D(const real* X, real* Y, index_t start, index_t stop) const;
    
    /// 3D isotropic multiplication of a vector: Y <- Y + M * X with dim(X) = 3 * dim(M)
    void vecMulAddIso3D(const real* X, real* Y, index_t start, index_t stop) const;

    /// multiplication of a vector, for columns within [start, stop[
    void vecMul(const real* X, real* Y, index_t start, index_t stop) const;

    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M)
    void vecMulAdd(const real* X, real* Y)      const { vecMulAdd(X, Y, 0, size_); }
    
    /// 2D isotropic multiplication of a vector: Y <- Y + M * X with dim(X) = 2 * dim(M)
    void vecMulAddIso2D(const real* X, real* Y) const { vecMulAddIso2D(X, Y, 0, size_); }
    
    /// 3D isotropic multiplication of a vector: Y <- Y + M * X with dim(X) = 3 * dim(M)
    void vecMulAddIso3D(const real* X, real* Y) const { vecMulAddIso3D(X, Y, 0, size_); }

    /// multiplication of a vector: Y <- M * X with dim(X) = dim(M)
    void vecMul(const real* X, real* Y)         const { vecMul(X, Y, 0, size_); }
    
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M)
    void vecMulAdd_ALT(const real* X, real* Y, index_t start, index_t stop) const;
    
    /// multiplication of a vector: Y <- Y + M * X with dim(X) = dim(M)
    void vecMulAdd_ALT(const real* X, real* Y)  const { vecMulAdd_ALT(X, Y, 0, size_); }

    /// true if matrix is non-zero
    bool notZero() const;
    
    /// number of elements in columns [start, stop[
    size_t nbElements(index_t start, index_t stop, size_t& alc) const;
    
    /// number of diagonal elements in columns [start, stop[
    size_t nbDiagonalElements(index_t start, index_t stop) const;

    /// number of elements in matrix
    size_t nbElements() const { size_t alc; return nbElements(0, size_, alc); }

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

    /// printf debug function in sparse mode: i, j : value
    void printSparseArray(std::ostream&) const;
    
    /// debug function
    int bad() const;
};


#endif

