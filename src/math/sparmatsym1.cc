// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#include <cmath>
#include <iomanip>
#include <sstream>
#include <iostream>

#include "sparmatsym1.h"
#include "assert_macro.h"
#include "blas.h"
#include "simd.h"

#ifdef __AVX__
#  define SPARMAT1_USES_AVX REAL_IS_DOUBLE
#  define SPARMAT1_USES_SSE REAL_IS_DOUBLE
#elif USE_SIMD
#  define SPARMAT1_USES_AVX 0
#  define SPARMAT1_USES_SSE REAL_IS_DOUBLE
#else
#  define SPARMAT1_USES_AVX 0
#  define SPARMAT1_USES_SSE 0
#endif


SparMatSym1::SparMatSym1()
{
    size_   = 0;
    alloc_  = 0;
    column_ = nullptr;
    colsiz_ = nullptr;
    colmax_ = nullptr;
    diagon_ = nullptr;
#if SPARMAT1_COMPACTED
    nmax_ = 0;
    off_  = new index_t[2]();
    ija_  = nullptr;
    elm_  = nullptr;
#endif
#if SPARMAT1_USES_COLNEXT
    colidx_ = new index_t[2]();
#endif
}


void SparMatSym1::allocate(index_t alc)
{
    if ( alc > alloc_ )
    {
        /*
         'chunk' can be increased to gain performance:
         more memory will be used, but reallocation will be less frequent
         */
        constexpr index_t chunk = 32;
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );

        //fprintf(stderr, "SMS1 allocate matrix %lu\n", alc);
        Element ** column_new = new Element*[alc];
        index_t * colsiz_new = new index_t[alc];
        index_t * colmax_new = new index_t[alc];
        real * diagon_new = new_real(alc);
        
        index_t ii = 0;
        if ( column_ )
        {
            for ( ; ii < alloc_; ++ii )
            {
                column_new[ii] = column_[ii];
                colsiz_new[ii] = colsiz_[ii];
                colmax_new[ii] = colmax_[ii];
                diagon_new[ii] = diagon_[ii];
            }
            delete[] column_;
            delete[] colsiz_;
            delete[] colmax_;
            free_real(diagon_);
        }
        
        column_ = column_new;
        colsiz_ = colsiz_new;
        colmax_ = colmax_new;
        diagon_ = diagon_new;
        alloc_ = alc;

        for ( ; ii < alc; ++ii )
        {
            column_[ii] = nullptr;
            colsiz_[ii] = 0;
            colmax_[ii] = 0;
            diagon_[ii] = 0;
        }
#if SPARMAT1_COMPACTED
        delete [] off_;
        off_ = new index_t[alc+2];
#endif
#if SPARMAT1_USES_COLNEXT
        delete[] colidx_;
        colidx_ = new index_t[alc+2];
        for ( index_t n = 0; n <= alc; ++n )
            colidx_[n] = n;
#endif
    }
}


index_t SparMatSym1::allocated() const
{
    index_t res = 0;
    if ( column_ )
        for ( index_t i = 0; i < alloc_; ++i )
            res += colmax_[i];
    return res;
}
    

void SparMatSym1::deallocate()
{
    if ( column_ )
    {
        for ( index_t ii = 0; ii < alloc_; ++ii )
            delete[] column_[ii];
        delete[] column_; column_ = nullptr;
        delete[] colsiz_; colsiz_ = nullptr;
        delete[] colmax_; colmax_ = nullptr;
        free_real(diagon_); diagon_ = nullptr;
#if SPARMAT1_USES_COLNEXT
        delete[] colidx_; colidx_ = nullptr;
#endif
#if SPARMAT1_COMPACTED
        delete[] ija_;
        free_real(elm_);
        ija_ = nullptr;
        elm_ = nullptr;
#endif
    }
    alloc_ = 0;
}


/// copy `cnt` elements from `src` to `dst`
static void copy(index_t cnt, SparMatSym1::Element * src, SparMatSym1::Element * dst)
{
    for ( index_t ii = 0; ii < cnt; ++ii )
        dst[ii] = src[ii];
}


void SparMatSym1::allocateColumn(const index_t j, index_t alc)
{
    assert_true( j < size_ );
    if ( alc > colmax_[j] )
    {
        constexpr index_t chunk = 4;
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );
        //fprintf(stderr, "SMS1 column %lu: %u --> %u\n", j, colmax_[j], alc);

        Element * ptr = new Element[alc];
        if ( column_[j] )
        {
            //copy over previous column elements
            copy(colsiz_[j], column_[j], ptr);
            
            //release old memory
            delete[] column_[j];
        }
        column_[j] = ptr;
        colmax_[j] = alc;
        assert_true( alc == colmax_[j] );
    }
}


void SparMatSym1::deallocateColumn(const index_t j)
{
    //fprintf(stderr, "SMS1 col %lu: %u --> null\n", j, colmax_[j]);
    delete[] column_[j];
    column_[j] = nullptr;
    colmax_[j] = 0;
    colsiz_[j] = 0;
}


/**
 This should work for any value of (i, j)
 This allocates to be able to hold the matrix element if necessary
 The diagonal element is always first in each column
*/
real& SparMatSym1::element(index_t ii, index_t jj)
{
    assert_true( ii > jj );
    assert_true( jj < size_ );
    //fprintf(stderr, "SMS1( %6i %6i )\n", i, j);

    assert_true( colsiz_[jj] <= colmax_[jj] );
    // make space for one additional element:
    if ( 1 + colsiz_[jj] > colmax_[jj] )
        allocateColumn(jj, 1 + colsiz_[jj]);
    
    // the column pointer should not change anymore now:
    Element * col = column_[jj];
    // place value, beyond range, acting as a fence for linear search below:
    Element * fence = col + colsiz_[jj];
    fence->reset(ii);
    // search column, given that elements are kept ordered in the column:
    Element * e = col;
    while ( e->inx < ii )
        ++e;
    
    if ( e->inx == ii )
    {
        colsiz_[jj] += ( e == fence );
        //printColumn(std::clog, jj);
        return e->val;
    }
    
    assert_true( colsiz_[jj]+1 <= colmax_[jj] );
    // shift all values above 'e', including 'e':
    col += colsiz_[jj];
    while ( --col >= e )
        col[1] = col[0];

    // add the requested term
    e->reset(ii);
    ++colsiz_[jj];
    //printColumn(std::clog, jj);
    return e->val;
}


real* SparMatSym1::address(index_t i, index_t j) const
{
    if ( i == j )
        return diagon_ + i;
    
    // swap to get ii <= jj (address lower triangle)
    index_t ii = std::max(i, j);
    index_t jj = std::min(i, j);

    for ( index_t kk = 0; kk < colsiz_[jj]; ++kk )
        if ( column_[jj][kk].inx == ii )
            return &( column_[jj][kk].val );
    
    return nullptr;
}


//------------------------------------------------------------------------------
#pragma mark -

void SparMatSym1::reset()
{
    for ( index_t j = 0; j < size_; ++j )
    {
#if 0
        Element * col = column_[j];
        for ( index_t i = 0; i < colsiz_[j]; ++i )
            col[i].val = 0;
#endif
        colsiz_[j] = 0;
        diagon_[j] = 0;
    }
}


bool SparMatSym1::notZero() const
{
    //check for any non-zero sparse term:
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        if ( diagon_[jj] != 0 )
            return true;
        for ( index_t kk = 0 ; kk < colsiz_[jj] ; ++kk )
        {
            assert_true(column_[jj][kk].val > jj);
            if ( column_[jj][kk].val != 0 )
                return true;
        }
    }
    //if here, the matrix is empty
    return false;
}


void SparMatSym1::scale(const real alpha)
{
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        diagon_[jj] *= alpha;
        for ( index_t n = 0; n < colsiz_[jj]; ++n )
            column_[jj][n].val *= alpha;
    }
}


void SparMatSym1::addDiagonalBlock(real* mat, const index_t ldd, index_t start, index_t cnt,
                                   const index_t mul) const
{
    assert_true( start + cnt <= size_ );
    
    for ( index_t j = 0; j < cnt; ++j )
    {
        Element * col = column_[j+start];
        real * dst = mat + ( j + ldd * j ) * mul;
        dst[0] += diagon_[j+start];
        for ( index_t n = 0; n < colsiz_[j+start]; ++n )
        {
            // assuming lower triangle is stored:
            assert_true( col[n].inx > j + start );
            index_t ij = col[n].inx - ( j + start );
            if ( j + ij < cnt )
            {
                //printf("SMS2 %4i %4i % .4f\n", j, ij, a);
                dst[mul*ij] += col[n].val;
                dst[mul*ldd*ij] += col[n].val;
            }
        }
    }
}


/*
addresses `mat' using lower banded storage for a symmetric matrix
mat(i, j) is stored in mat[i-j+ldd*j]
*/
void SparMatSym1::addLowerBand(real alpha, real* mat, const index_t ldd, index_t start, index_t cnt,
                               const index_t mul, const index_t rank) const
{
    start *= mul;
    cnt *= mul;
    assert_true( start + cnt <= size_ );

    for ( index_t j = 0; j < cnt; ++j )
    {
        Element * col = column_[j+start];
        index_t sup = std::min(cnt-j, rank+1);
        real * dst = mat + ( j + ldd * j ) * mul;
        dst[0] += alpha * diagon_[j+start];
        for ( index_t n = 0; n < colsiz_[j+start]; ++n )
        {
            // assuming lower triangle is stored:
            assert_true( col[n].inx >= j + start );
            index_t ij = col[n].inx - ( j + start );
            if ( ij < sup )
            {
                //printf("SMS2 %4i %4i % .4f\n", ii, jj, a);
                dst[ij] += alpha * col[n].val;
            }
        }
    }
}


void SparMatSym1::addDiagonalTrace(real alpha, real* mat, const index_t ldd, index_t start, index_t cnt,
                                   const index_t mul, const index_t rank, const bool sym) const
{
    start *= mul;
    index_t end = start + mul * cnt;
    assert_true( cnt <= size_ );

    for ( index_t jj = start; jj < end; ++jj )
    {
        index_t j = ( jj - start ) / mul;
        // with banded storage, mat(i, j) is stored in mat[i-j+ldd*j]
        mat[j+ldd*j] += alpha * diagon_[jj];
        for ( index_t n = 0; n < colsiz_[jj]; ++n )
        {
            assert_true( column_[jj][n].inx > jj );
            // assuming lower triangle is stored:
            index_t ii = column_[jj][n].inx - start;
            index_t i = ii / mul;
            if (( ii-jj == (i-j)*mul ) & ( i < cnt ) & ( i <= j + rank ))
            {
                real a = alpha * column_[jj][n].val;
                //fprintf(stderr, "SMS1 %4lu %4lu : %.4f\n", i, j, a);
                // with banded storage, mat(i, j) is stored in mat[i-j+ldd*j]
                mat[i+ldd*j] += a;
                if ( sym ) mat[j+ldd*i] += a;
            }
        }
    }
}


int SparMatSym1::bad() const
{
    if ( size_ <= 0 ) return 1;
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        if ( diagon_[jj] != diagon_[jj] )
            return 7;
        for ( index_t kk = 0 ; kk < colsiz_[jj] ; ++kk )
        {
            index_t ii = column_[jj][kk].inx;
            if ( ii >= size_ ) return 2;
            if ( ii <= jj ) return 3;
        }
    }
    return 0;
}


//all allocated elements are counted, even if zero
size_t SparMatSym1::nbElements(index_t start, index_t stop, size_t& alc) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_);
#if SPARMAT1_COMPACTED
    alc = nmax_ + alloc_;
#else
    alc = alloc_;
#endif
    size_t cnt = 0; // diagonal elements
    for ( index_t j = start; j < stop; ++j )
    {
        alc += colmax_[j];
        cnt += colsiz_[j];
    }
    return cnt;
}


size_t SparMatSym1::nbDiagonalElements(index_t start, index_t stop) const
{
    return stop - start;
}


std::string SparMatSym1::what() const
{
    size_t alc;
    size_t cnt = nbElements(0, size_, alc);
    std::ostringstream msg;
#if SPARMAT1_COMPACTED
    msg << "c";
#endif
#if SPARMAT1_USES_AVX
    msg << "SMS1x ";
#elif SPARMAT1_USES_SSE
    msg << "SMS1e ";
#else
    msg << "SMS1 ";
#endif
    msg << cnt << " (" << alc << ")";
    return msg.str();
}


void SparMatSym1::printSparse(std::ostream& os, real inf, index_t start, index_t stop) const
{
    stop = std::min(stop, size_);
    os << "% SparMatSym1 size " << size_ << ":\n";
    for ( index_t jj = start; jj < stop; ++jj )
    {
        if ( diagon_[jj] != 0 || colsiz_[jj] > 0 )
        {
            os << "% column " << jj << " dia " << diagon_[jj] << "\n";
            for ( index_t n = 0 ; n < colsiz_[jj] ; ++n )
            {
                index_t ii = column_[jj][n].inx;
                real v = column_[jj][n].val;
                if ( abs_real(v) >= inf )
                    os << std::setw(6) << ii << " " << std::setw(6) << jj << " " << std::setw(16) << v << "\n";
            }
        }
    }
}


void SparMatSym1::printSummary(std::ostream& os, index_t start, index_t stop)
{
    stop = std::min(stop, size_);
    os << "% SparMatSym1 size " << size_ << ":";
    for ( index_t jj = start; jj < stop; ++jj )
    {
        os << "\n   " << jj << "   " << colsiz_[jj];
#if SPARMAT1_USES_COLNEXT
        os << " index " << colidx_[jj];
#endif
    }
    std::endl(os);
}


void SparMatSym1::printColumn(std::ostream& os, const index_t jj)
{
    std::streamsize p = os.precision(1);
    Element const* col = column_[jj];
    os << "SMS1 col " << jj << " (" << diagon_[jj] << ") : ";
    for ( index_t n = 0; n < colsiz_[jj]; ++n )
    {
        real v = col[n].val;
        if ( v == 0 )
            os << " !";  // this is a waste
        else if ( abs_real(v) < REAL_EPSILON )
            os << " *";  // this element could be removed
        else
            os << " ";
        os << col[n].inx << " (" << std::fixed << v << ") ";
    }
    std::endl(os);
    os.precision(p);
}


void SparMatSym1::printSparseArray(std::ostream& os) const
{
#if SPARMAT1_COMPACTED
    if ( ija_ )
    {
        std::ios::fmtflags fgs = os.flags();
        
        os << "off ";
        for ( index_t j = 0; j <= size_; ++j )
            os << " " << std::setw(4) << off_[j];
        os << "\n";

        index_t end = off_[size_];
        os << "ija ";
        for ( index_t n = 0; n < end; ++n )
            os << " " << std::setw(4) << ija_[n];
        os << "\n";
        
        std::streamsize p = os.precision(2);
        os << "elm ";
        for ( index_t n = 0; n < end; ++n )
            os << " " << std::setw(4) << elm_[n];
        os << "\n";
        os.precision(p);
        os.setf(fgs);
        return;
    }
#endif
    os << "SMS1 has no compact storage\n";
}


//------------------------------------------------------------------------------
#pragma mark - Column-Vector Multiplication

/**
Multiply by column `jj` with off-diagonal elements provided in `col` of size `cnt`
*/
void SparMatSym1::vecMulAddCol(const real* X, real* Y, index_t jj,
                               real dia, Element col[], index_t cnt)
{
    const real X0 = X[jj];
    real Y0 = Y[jj] + dia * X0;
    for ( index_t n = 0; n < cnt; ++n )
    {
        const index_t ii = col[n].inx;
        const real a = col[n].val;
        assert_true( ii > jj );
        Y[ii] += a * X0;
        Y0 += a * X[ii];
    }
    Y[jj] = Y0;
}

/**
 Multiply by column `jj` with off-diagonal elements provided in `col` of size `cnt`
 2 x multiplexed
 */
void SparMatSym1::vecMulAddColIso2D(const real* X, real* Y, index_t jj,
                                    real dia, Element col[], index_t cnt)
{
    const real X0 = X[jj  ];
    const real X1 = X[jj+1];
    real Y0 = Y[jj  ] + dia * X0;
    real Y1 = Y[jj+1] + dia * X1;
    for ( index_t n = 0; n < cnt; ++n )
    {
        const index_t ii = 2 * col[n].inx;
        const real a = col[n].val;
        assert_true( ii > jj );
        Y[ii  ] += a * X0;
        Y[ii+1] += a * X1;
        Y0 += a * X[ii  ];
        Y1 += a * X[ii+1];
    }
    Y[jj  ] = Y0;
    Y[jj+1] = Y1;
}

/**
 Multiply by column `jj` with off-diagonal elements provided in `col` of size `cnt`
 3 x multiplexed
*/
void SparMatSym1::vecMulAddColIso3D(const real* X, real* Y, index_t jj,
                                    real dia, Element col[], index_t cnt)
{
    const real X0 = X[jj  ];
    const real X1 = X[jj+1];
    const real X2 = X[jj+2];
    real Y0 = Y[jj  ] + dia * X0;
    real Y1 = Y[jj+1] + dia * X1;
    real Y2 = Y[jj+2] + dia * X2;
    for ( index_t n = 0; n < cnt; ++n )
    {
        const index_t ii = 3 * col[n].inx;
        const real a = col[n].val;
        assert_true( ii > jj );
        Y[ii  ] += a * X0;
        Y[ii+1] += a * X1;
        Y[ii+2] += a * X2;
        Y0 += a * X[ii  ];
        Y1 += a * X[ii+1];
        Y2 += a * X[ii+2];
    }
    Y[jj  ] = Y0;
    Y[jj+1] = Y1;
    Y[jj+2] = Y2;
}


//------------------------------------------------------------------------------
#pragma mark - Prepare Multiplication

void SparMatSym1::setColumnIndex()
{
#if SPARMAT1_USES_COLNEXT
    colidx_[size_] = size_;
    if ( size_ > 0 )
    {
        index_t inx = size_;
        index_t nxt = size_;
        while ( inx-- > 0 )
        {
            if ( diagon_[inx] != 0 || colsiz_[inx] > 0 )
                nxt = inx;
            else if ( colsiz_[inx] == 0 )
                deallocateColumn(inx);
            colidx_[inx] = nxt;
        }
    }
#endif
#if ( 0 )
    std::clog << std::endl;
    index_t cnt = 0;
    for ( index_t j = 0; j < size_; ++j )
    {
        cnt += ( colsiz_[j] == 0 );
        //printColumn(std::clog, j);
    }
    std::clog << "SMS1 has " << cnt << " / " << size_ << " empty columns\n";
#endif
}


#if !SPARMAT1_COMPACTED

bool SparMatSym1::prepareForMultiply(int)
{
    setColumnIndex();
    return true;
}

#else


/**
 Create the Compressed Sparse Row (CSR) format described in
     "SPARSKIT: a basic tool kit for sparse matrix computations"
 Youcef Saad, Version 2, 1994
 
 Attention: indices start here at zero, and many things are thus shifted by one,
 compared to FORTRAN code and many related documentation.

 * dia[] stores the diagonal elements of the matrix.
 Off-diagonal element are stored following the Compressed Sparse Row format,
 using three one-dimensional arrays: off[], elm[] and ija[].
 * elm[] stores matrix element values as real
 * ija[] stores the line index corresponding to elm[] at the same position.
 * off[j] is the index in elm[] and ija[] of the first off-element of column j

 Thus:
 * Array dia[] is of size N and off[] is of size N+1.
 * off[size] = number of off-diagonal elements.
 * Arrays ija[] and elm[] hold all off-diagonal non-zero elements.
 the off-diagonal elements in row i are in elm[k] where ija[i] <= k < ija[i+1]
 
 Example:
 
    3: 0: 1: 0: 0:
    0: 4: 0: 0: 0:
    0: 7: 5: 9: 0:
    0: 0: 0: 0: 2:
    0: 0: 0: 6: 5:

    In row-indexed compact storage, this 5x5 matrix is represented as follows:
    dia[]  3. 4. 5. 0. 5.
    off[]  0  1  1  3  4  5
    ija[]  2  1  3  4  3
    elm[]  1. 7. 9. 2. 6.

    Since the number of off-diagonal non-zero elements is 5.
 */

bool SparMatSym1::prepareForMultiply(int dim)
{
    assert_true( size_ <= alloc_ );
    
    setColumnIndex();

    // count number of non-zero elements
    index_t nnz = 0;
    for ( index_t jj = 0; jj < size_; ++jj )
        nnz += colsiz_[jj];
    
    //allocate classical sparse matrix storage (Numerical Recipes)
    if ( nnz > nmax_ )
    {
        constexpr index_t chunk = 1;
        nmax_ = ( nnz + chunk - 1 ) & ~( chunk - 1 );
        delete[] ija_;
        free_real(elm_);
        ija_ = new unsigned[nmax_];
        elm_ = new_real(nmax_);
    }
    
    // Create the Compressed Sparse Storage
    unsigned cnt = 0;
    off_[0] = 0;
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        Element * col = column_[jj];
        for ( index_t k = 0; k < colsiz_[jj]; ++k )
        {
            assert_true( col[k].inx > jj );
            if ( col[k].val != 0 )
            {
                // non-zero non-diagonal element
                assert_true( cnt < nnz );
                elm_[cnt] = col[k].val;
                ija_[cnt] = dim * col[k].inx;
                ++cnt;
            }
        }
        off_[jj+1] = cnt;
    }
    if ( cnt > nnz )
        ABORT_NOW("internal error");
    //std::clog << "SMS1 has " << cnt << " (" << nnz << ") non-zero elements\n";
    
    //printSparse(std::clog);
    //printSparseArray(std::cout);
    return true;
}

//------------------------------------------------------------------------------
#pragma mark - Optimized Column-Vector multiplication


void SparMatSym1::vecMulAddCol(const real* X, real* Y, index_t jj,
                               const real dia, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= nmax_ );
    real X0 = X[jj];
    real Y0 = Y[jj] + dia * X0;
    for ( index_t n = start; n < stop; ++n )
    {
        real a = elm_[n];
        auto ii = ija_[n];
        Y[ii] += a * X0;
        Y0 += a * X[ii];
    }
    Y[jj] = Y0;
}

void SparMatSym1::vecMulAddColIso2D(const real* X, real* Y, index_t jj,
                                    const real dia, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= nmax_ );
    real X0 = X[jj  ];
    real X1 = X[jj+1];
    real Y0 = Y[jj  ] + dia * X0;
    real Y1 = Y[jj+1] + dia * X1;
    for ( index_t n = start; n < stop; ++n )
    {
        real a = elm_[n];
        auto ii = ija_[n];
        assert_true( ii > jj );
        Y0 += a * X[ii  ];
        Y1 += a * X[ii+1];
        Y[ii  ] += a * X0;
        Y[ii+1] += a * X1;
    }
    Y[jj  ] = Y0;
    Y[jj+1] = Y1;
}


void SparMatSym1::vecMulAddColIso3D(const real* X, real* Y, index_t jj,
                                    const real dia, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= nmax_ );
    real X0 = X[jj  ];
    real X1 = X[jj+1];
    real X2 = X[jj+2];
    real Y0 = Y[jj  ] + dia * X0;
    real Y1 = Y[jj+1] + dia * X1;
    real Y2 = Y[jj+2] + dia * X2;
    for ( index_t n = start; n < stop; ++n )
    {
        auto ii = ija_[n];
        assert_true( ii > jj );
        real a = elm_[n];
        Y0 += a * X[ii  ];
        Y1 += a * X[ii+1];
        Y2 += a * X[ii+2];
        Y[ii  ] += a * X0;
        Y[ii+1] += a * X1;
        Y[ii+2] += a * X2;
    }
    Y[jj  ] = Y0;
    Y[jj+1] = Y1;
    Y[jj+2] = Y2;
}

//------------------------------------------------------------------------------
#pragma mark - 2D SIMD

#if SPARMAT1_USES_SSE

static inline void multiply2(const double* X, double* Y, index_t ii,
                             const double* val, vec2 const& xx, vec2& ss)
{
    vec2 aa = loaddup2(val);
    ss = fmadd2(load2(X+ii), aa, ss);
    store2(Y+ii, fmadd2(xx, aa, load2(Y+ii)));
}


void SparMatSym1::vecMulAddColIso2D_SSE(const double* X, double* Y, index_t jj,
                                        const double* dia, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= nmax_ );
    const vec2 xx = load2(X+jj);
    vec2 ss = fmadd2(loaddup2(dia), xx, load2(Y+jj));
    // there is a dependence here for 'ss'
    for ( index_t n = start; n < stop; ++n )
        multiply2(X, Y, ija_[n], elm_+n, xx, ss);
    store2(Y+jj, ss);
}


void SparMatSym1::vecMulAddColIso2D_SSEU(const double* X, double* Y, index_t jj,
                                         const double* dia, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= nmax_ );
    const vec2 xx = load2(X+jj);
    vec2 s0 = mul2(loaddup2(dia), xx);
    vec2 s1 = load2(Y+jj);
    vec2 s2 = setzero2();
    vec2 s3 = setzero2();
    
    index_t n = start;
#if ( 0 )
    // unrolling by 8 may exceed the number of registers in the CPU
#pragma nounroll
    if ( end >= n + 8 )
    {
        vec2 s4 = setzero2();
        vec2 s5 = setzero2();
        vec2 s6 = setzero2();
        vec2 s7 = setzero2();
        index_t end = n + 8 * ( ( stop - n ) / 8 );
        // process 8 by 8:
        for ( ; n < end; n += 8 )
        {
            const auto i0 = ija_[n  ];
            const auto i1 = ija_[n+1];
            const auto i2 = ija_[n+2];
            const auto i3 = ija_[n+3];
            const auto i4 = ija_[n+4];
            const auto i5 = ija_[n+5];
            const auto i6 = ija_[n+6];
            const auto i7 = ija_[n+7];
            vec2 y0 = load2(Y+i0);
            vec2 y1 = load2(Y+i1);
            vec2 y2 = load2(Y+i2);
            vec2 y3 = load2(Y+i3);
            vec2 y4 = load2(Y+i4);
            vec2 y5 = load2(Y+i5);
            vec2 y6 = load2(Y+i6);
            vec2 y7 = load2(Y+i7);
            vec2 a0 = loaddup2(elm_+n);
            vec2 a1 = loaddup2(elm_+n+1);
            vec2 a2 = loaddup2(elm_+n+2);
            vec2 a3 = loaddup2(elm_+n+3);
            vec2 a4 = loaddup2(elm_+n+4);
            vec2 a5 = loaddup2(elm_+n+5);
            vec2 a6 = loaddup2(elm_+n+6);
            vec2 a7 = loaddup2(elm_+n+7);
            s0 = fmadd2(load2(X+i0), a0, s0);
            s1 = fmadd2(load2(X+i1), a1, s1);
            s2 = fmadd2(load2(X+i2), a2, s2);
            s3 = fmadd2(load2(X+i3), a3, s3);
            s4 = fmadd2(load2(X+i4), a4, s4);
            s5 = fmadd2(load2(X+i5), a5, s5);
            s6 = fmadd2(load2(X+i6), a6, s6);
            s7 = fmadd2(load2(X+i7), a7, s7);
            store2(Y+i0, fmadd2(xx, a0, y0));
            store2(Y+i1, fmadd2(xx, a1, y1));
            store2(Y+i2, fmadd2(xx, a2, y2));
            store2(Y+i3, fmadd2(xx, a3, y3));
            store2(Y+i4, fmadd2(xx, a4, y4));
            store2(Y+i5, fmadd2(xx, a5, y5));
            store2(Y+i6, fmadd2(xx, a6, y6));
            store2(Y+i7, fmadd2(xx, a7, y7));
        }
        // collapse into lower summation registers:
        s0 = add2(s0, s4);
        s1 = add2(s1, s5);
        s2 = add2(s2, s6);
        s3 = add2(s3, s7);
    }
#endif
    
    index_t end = n + 4 * ( ( stop - n ) / 4 );
    // process 4 by 4:
#pragma nounroll
    for ( ; n < end; n += 4 )
    {
#if ( 0 )
        /*
         Since all the indices are different, the blocks can be processed in
         parallel, and micro-operations can be interleaved to avoid latency.
         The compiler however cannot assume this, because the indices of the
         blocks are not known at compile time.
         */
        multiply2(X, Y, ija_[n  ], elm_+n  , xx, s0);
        multiply2(X, Y, ija_[n+1], elm_+n+1, xx, s1);
        multiply2(X, Y, ija_[n+2], elm_+n+2, xx, s2);
        multiply2(X, Y, ija_[n+3], elm_+n+3, xx, s3);
#else
        /* we remove here the apparent dependency on the values of Y[],
         which are read and written, but at different indices.
         The compiler can reorder instructions to avoid lattencies */
        const auto i0 = ija_[n  ];
        const auto i1 = ija_[n+1];
        const auto i2 = ija_[n+2];
        const auto i3 = ija_[n+3];
        vec2 y0 = load2(Y+i0);
        vec2 y1 = load2(Y+i1);
        vec2 y2 = load2(Y+i2);
        vec2 y3 = load2(Y+i3);
        vec2 a0 = loaddup2(elm_+n);
        vec2 a1 = loaddup2(elm_+n+1);
        vec2 a2 = loaddup2(elm_+n+2);
        vec2 a3 = loaddup2(elm_+n+3);
        s0 = fmadd2(load2(X+i0), a0, s0);
        s1 = fmadd2(load2(X+i1), a1, s1);
        s2 = fmadd2(load2(X+i2), a2, s2);
        s3 = fmadd2(load2(X+i3), a3, s3);
        store2(Y+i0, fmadd2(xx, a0, y0));
        store2(Y+i1, fmadd2(xx, a1, y1));
        store2(Y+i2, fmadd2(xx, a2, y2));
        store2(Y+i3, fmadd2(xx, a3, y3));
#endif
    }
    // collapse 's0'
    s0 = add2(add2(s0,s1), add2(s2,s3));
    // process remaining blocks:
#pragma nounroll
    for ( ; n < stop; ++n )
        multiply2(X, Y, ija_[n], elm_+n, xx, s0);
    store2(Y+jj, s0);
}

#endif

#if SPARMAT1_USES_AVX && SPARMAT1_COMPACTED

/*
Accumulation is done here in the higher part of 'ss'
The high position of 'xx' is not used
The low position of 'ss' is used locally
*/
static inline void multiply4(const double* X, double* Y, index_t ii,
                             const double* val, vec4 const& xx, vec4& ss)
{
    vec4 x = blend22(xx, broadcast2(X+ii));  // hi <- X , lo <- xx
    ss = blend22(load2crap(Y+ii), ss);    // hi <- ss, lo <- Y
    ss = fmadd4(broadcast1(val), x, ss);
    store2(Y+ii, getlo(ss));
}


void SparMatSym1::vecMulAddColIso2D_AVX(const double* X, double* Y, index_t jj,
                                        const double* dia, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= nmax_ );
    const vec4 xx = broadcast2(X+jj);  // hi position
    vec4 ss = fmadd4(broadcast1(dia), xx, broadcast2(Y+jj));
    // there is a dependence here for 'ss'
    for ( index_t n = start; n < stop; ++n )
        multiply4(X, Y, ija_[n], elm_+n, xx, ss);
    store2(Y+jj, gethi(ss));
}


void SparMatSym1::vecMulAddColIso2D_AVXU(const double* X, double* Y, index_t jj,
                                         const double* dia, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= nmax_ );
    const vec4 xx = broadcast2(X+jj);  // hi and lo position
    vec4 s0 = mul4(broadcast1(dia), xx);
    vec4 s1 = broadcast2(Y+jj);
    vec4 s2 = setzero4();
    vec4 s3 = setzero4();
    
    unsigned * inx = ija_ + start;
    const double * val = elm_ + start;
    const double * end = elm_ + stop;
    const double * halt = end - 3;  // val+3 <= end-1  is  val < end-3;
        // process 4 by 4:
    #pragma nounroll
    for ( ; val < halt; val += 4, inx += 4 )
    {
#if ( 0 )
        /*
         Since all the indices are different, the blocks can be processed in
         parallel, and micro-operations can be interleaved to avoid latency.
         The compiler however cannot assume this, because the indices of the
         blocks are not known at compile time.
         */
        multiply4(X, Y, inx[0], val  , xx, s0);
        multiply4(X, Y, inx[1], val+1, xx, s1);
        multiply4(X, Y, inx[2], val+2, xx, s2);
        multiply4(X, Y, inx[3], val+3, xx, s3);
#else
        /* we remove here the apparent dependency on the values of Y[],
         which are read and written, but at different indices.
         The compiler can reorder instructions to avoid lattencies */
        //__m128i ii = _mm_slli_epi32(_mm_loadu_si128((__m128i*)(ija_+n)), 0x1);
        //printi(ii, "indx");
        const auto i0 = inx[0];
        const auto i1 = inx[1];
        const auto i2 = inx[2];
        const auto i3 = inx[3];
        s0 = blend22(load2crap(Y+i0),s0);    // lo = Y
        s1 = blend22(load2crap(Y+i1),s1);    // lo = Y
        s2 = blend22(load2crap(Y+i2),s2);    // lo = Y
        s3 = blend22(load2crap(Y+i3),s3);    // lo = Y
        vec4 x0 = blend22(xx,broadcast2(X+i0));   // hi = X , lo <- xx
        vec4 x1 = blend22(xx,broadcast2(X+i1));   // hi = X , lo <- xx
        vec4 x2 = blend22(xx,broadcast2(X+i2));   // hi = X , lo <- xx
        vec4 x3 = blend22(xx,broadcast2(X+i3));   // hi = X , lo <- xx
        s0 = fmadd4(broadcast1(val  ), x0, s0);
        s1 = fmadd4(broadcast1(val+1), x1, s1);
        s2 = fmadd4(broadcast1(val+2), x2, s2);
        s3 = fmadd4(broadcast1(val+3), x3, s3);
        store2(Y+i0, getlo(s0));
        store2(Y+i1, getlo(s1));
        store2(Y+i2, getlo(s2));
        store2(Y+i3, getlo(s3));
#endif
    }
    // collapse into 's0'
    s0 = add4(add4(s0,s1), add4(s2,s3));
    // process remaining values:
#pragma nounroll
    for ( ; val < end; ++val, ++inx )
        multiply4(X, Y, inx[0], val, xx, s0);
    store2(Y+jj, gethi(s0));
}

#endif
#endif

//------------------------------------------------------------------------------
#pragma mark - 3D SIMD

#if SPARMAT1_USES_SSE
#endif


#if SPARMAT1_USES_AVX && SPARMAT1_COMPACTED
void SparMatSym1::vecMulAddColIso3D_AVX(const double* X, double* Y, index_t jj,
                                        const double* dia, index_t start, index_t stop) const
{
    //fprintf(stderr, "col %lu: %lu %lu\n", jj, start, stop);
    assert_true( start <= stop );
    assert_true( stop <= nmax_ );
    vec4 zz = setzero4();
    vec4 xx = blend31(loadu4(X+jj), zz);
    vec4 yy = fmadd4(broadcast1(dia), xx, loadu4(Y+jj));
    real * val = elm_ + start + 1;
    real const*end = elm_ + stop;
    unsigned *inx = ija_ + start;
    while ( val < end )
    {
        index_t ii = inx[0];
        index_t kk = inx[1];
        assert_true( kk > ii );
        inx += 2;
        vec4 aa = broadcast1(val-1);
        vec4 bb = broadcast1(val);
        vec4 nn = loadu4(Y+ii);
        vec4 mm = loadu4(Y+kk);
        val += 2;
        yy = fmadd4(aa, loadu4(X+ii), yy);
        zz = fmadd4(bb, loadu4(X+kk), zz);
        storeu4(Y+ii, fmadd4(aa, xx, nn));
        storeu4(Y+kk, fmadd4(bb, xx, mm));
    }
    yy = add4(yy, zz);
    while ( val <= end )
    {
        index_t ii = *inx++;
        assert_true( ii > jj );
        vec4 aa = broadcast1(val-1);
        ++val;
        yy = fmadd4(aa, loadu4(X+ii), yy);
        storeu4(Y+ii, fmadd4(aa, xx, loadu4(Y+ii)));
    }
    store3(Y+jj, yy);
}
#endif


//------------------------------------------------------------------------------
#pragma mark - Matrix-Vector Add-multiply


void SparMatSym1::vecMulAdd_ALT(const real* X, real* Y, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_);

    for ( index_t jj = start; jj < stop; ++jj )
        vecMulAddCol(X, Y, jj, diagon_[jj], column_[jj], colsiz_[jj]);
}


void SparMatSym1::vecMulAdd(const real* X, real* Y, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_);

#if SPARMAT1_USES_COLNEXT
    for ( index_t jj = colidx_[start]; jj < stop; jj = colidx_[jj+1] )
#else
    for ( index_t jj = start; jj < stop; ++jj )
#endif
    {
#if SPARMAT1_COMPACTED
        vecMulAddCol(X, Y, jj, diagon_[jj], off_[jj], off_[jj+1]);
#else
        vecMulAddCol(X, Y, jj, diagon_[jj], column_[jj], colsiz_[jj]);
#endif
    }
}


void SparMatSym1::vecMulAddIso2D(const real* X, real* Y, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_);

#if SPARMAT1_USES_COLNEXT
    for ( index_t jj = colidx_[start]; jj < stop; jj = colidx_[jj+1] )
#else
    for ( index_t jj = start; jj < stop; ++jj )
#endif
    {
#if SPARMAT1_COMPACTED
#  if SPARMAT1_USES_AVX
        vecMulAddColIso2D_AVXU(X, Y, 2*jj, diagon_+jj, off_[jj], off_[jj+1]);
#  elif SPARMAT1_USES_SSE
        vecMulAddColIso2D_SSEU(X, Y, 2*jj, diagon_+jj, off_[jj], off_[jj+1]);
#  else
        vecMulAddColIso2D(X, Y, 2*jj, diagon_[jj], off_[jj], off_[jj+1]);
#  endif
#else
        vecMulAddColIso2D(X, Y, 2*jj, diagon_[jj], column_[jj], colsiz_[jj]);
#endif

    }
}


void SparMatSym1::vecMulAddIso3D(const real* X, real* Y, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_);

#if SPARMAT1_USES_COLNEXT
    for ( index_t jj = colidx_[start]; jj < stop; jj = colidx_[jj+1] )
#else
    for ( index_t jj = start; jj < stop; ++jj )
#endif
    {
#if SPARMAT1_COMPACTED
#  if SPARMAT1_USES_AVX
        vecMulAddColIso3D_AVX(X, Y, 3*jj, diagon_+jj, off_[jj], off_[jj+1]);
#  else
        vecMulAddColIso3D(X, Y, 3*jj, diagon_[jj], off_[jj], off_[jj+1]);
#  endif
#else
        vecMulAddColIso3D(X, Y, 3*jj, diagon_[jj], column_[jj], colsiz_[jj]);
#endif
    }
}


//------------------------------------------------------------------------------
#pragma mark - Matrix-Vector multiplication

void SparMatSym1::vecMul(const real* X, real* Y, index_t start, index_t stop) const
{
    zero_real(stop-start, Y+start);
    vecMulAdd(X, Y, start, stop);
}
