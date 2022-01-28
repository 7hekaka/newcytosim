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

    /*
     The row-indexed scheme sets up two one-dimensional arrays: sa and ija.
     sa[] stores matrix element values as real
     ija[] stores integer values.
     The storage rules are:
     * The first N locations of sa[] store A’s diagonal matrix elements, in order.
       (Diagonal elements are stored even if they are zero; this is at most a slight storage
       inefficiency, since diagonal elements are nonzero in most realistic applications.)
     * Each of the first N locations of ija[] stores the index of the array sa[] that
       contains the first off-diagonal element of the corresponding row of the matrix.
       (If there are no off-diagonal elements for that row, it is one greater than
       the index in sa[] of the most recently stored element of a previous row.)
     * ija[0] is always equal to N+2. (It can be read to determine N.)
     * ija[N] is one greater than the index in sa[] of the last off-diagonal
      element of the last row. (It can be read to determine the number of nonzero
     elements in the matrix, or the number of elements in the arrays sa and ija.)
     * sa[N] of is not used and can be set arbitrarily.
     * Entries in sa[] at locations N+2 contain A’s off-diagonal values,
       ordered by rows and, within each row, ordered by columns.
      Entries in ija[] at locations N+2 contain the column number of the corresponding
     element in sa[].

    3: 0: 1: 0: 0:
    0: 4: 0: 0: 0:
    0: 7: 5: 9: 0:
    0: 0: 0: 0: 2:
    0: 0: 0: 6: 5:

    In row-indexed compact storage, this 5x5 matrix is represented as follows:
    ija[k]  6  7  7  9 10 11  2  1  3  4  3
    sa[k]   3. 4. 5. 0. 5. X  1. 7. 9. 2. 6.
     
    The two arrays are of size 5 + nnz + 1, since there are nnz=5 off-diagonal non-zero elements.
     
     Here X is an arbitrary value. Notice that, according to the storage rules, the value of N
     (namely 5) is N = ija[0]-1, and the length of each array is ija[N], namely 11.
     The diagonal element in row i is sa[i], and the off-diagonal elements in that row are in
     sa[k] where k loops from ija[i] to ija[i+1]-1, if the upper limit is greater or equal to
     the lower one (as in C’s for loops).
    */
    size_t nmax_;
    unsigned * ija_;
    real * sa_;

#endif
    
    /// One column multiplication of a vector
    void vecMulAddCol(const real* X, real* Y, Element col[], size_t cnt) const;
    
    /// One column multiplication of a vector, isotropic 2D version
    void vecMulAddColIso2D(const real* X, real* Y, Element col[], size_t cnt) const;
    
    /// One column multiplication of a vector, isotropic 3D version
    void vecMulAddColIso3D(const real* X, real* Y, Element col[], size_t cnt) const;


    /// One column multiplication of a vector
    void vecMulAddCol(const real* X, real* Y, size_t jj, real const* dia, size_t start, size_t stop) const;
    
    /// One column 2D isotropic multiplication of a vector
    void vecMulAddColIso2D(const real* X, real* Y, size_t jj, real const* dia, size_t start, size_t stop) const;
    
    /// One column 3D isotropic multiplication of a vector
    void vecMulAddColIso3D(const real* X, real* Y, size_t jj, real const* dia, size_t start, size_t stop) const;

    
    /// One column 2D isotropic multiplication of a vector
    void vecMulAddColIso2D_SSE(const double* X, double* Y, size_t jj, double const* dia, size_t start, size_t stop) const;
    
    /// One column 2D isotropic multiplication of a vector
    void vecMulAddColIso2D_SSEU(const double* X, double* Y, size_t jj, double const* dia, size_t start, size_t stop) const;

    
    /// One column 2D isotropic multiplication of a vector
    void vecMulAddColIso2D_AVX(const double* X, double* Y, size_t jj, double const* dia, size_t start, size_t stop) const;
    
    /// One column 2D isotropic multiplication of a vector
    void vecMulAddColIso2D_AVXU(const double* X, double* Y, size_t jj, double const* dia, size_t start, size_t stop) const;
    
    /// One column 3D isotropic multiplication of a vector
    void vecMulAddColIso3D_AVX(const double* X, double* Y, size_t jj, double const* dia, size_t start, size_t stop) const;

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
    real& diagonal(size_t i);
    
    /// returns the address of element at (x, y), allocating if necessary
    real& operator()(size_t x, size_t y);
    
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

