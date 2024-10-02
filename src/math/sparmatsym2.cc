// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#include <cmath>
#include <iomanip>
#include <sstream>
#include <iostream>

#include "sparmatsym2.h"
#include "assert_macro.h"
#include "blas.h"
#include "simd.h"

#ifdef __AVX__
#  define SPARMAT2_USES_AVX 1
#  define SPARMAT2_USES_SSE 1
#  include "simd_float.h"
#elif USE_SIMD
#  define SPARMAT2_USES_AVX 0
#  define SPARMAT2_USES_SSE 1
#  include "simd_float.h"
#else
#  define SPARMAT2_USES_AVX 0
#  define SPARMAT2_USES_SSE 0
#endif


SparMatSym2::SparMatSym2()
{
    size_   = 0;
    alloc_  = 0;
    column_ = nullptr;
    colsiz_ = nullptr;
    colmax_ = nullptr;
#if SPARMAT2_COMPACTED
    alcDSS_ = 0;
    colDSS_ = nullptr;
    rowDSS_ = nullptr;
    valDSS_ = nullptr;
#endif
#if SPARMAT2_USES_COLNEXT
    colidx_ = new unsigned[2]();
#endif
}


void SparMatSym2::allocate(index_t alc)
{
    if ( alc > alloc_ )
    {
        /*
         'chunk' can be increased to gain performance:
         more memory will be used, but reallocation will be less frequent
         */
        constexpr index_t chunk = 64;
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );

        //fprintf(stderr, "SMS2 allocate matrix %lu\n", alc);
        Element ** column_new = new Element*[alc];
        unsigned * colsiz_new = new unsigned[alc];
        unsigned * colmax_new = new unsigned[alc];
        
        index_t ii = 0;
        if ( column_ )
        {
            for ( ; ii < alloc_; ++ii )
            {
                column_new[ii] = column_[ii];
                colsiz_new[ii] = colsiz_[ii];
                colmax_new[ii] = colmax_[ii];
            }
            delete[] column_;
            delete[] colsiz_;
            delete[] colmax_;
        }
        
        column_ = column_new;
        colsiz_ = colsiz_new;
        colmax_ = colmax_new;
        alloc_ = alc;

        for ( ; ii < alc; ++ii )
        {
            column_[ii] = nullptr;
            colsiz_[ii] = 0;
            colmax_[ii] = 0;
        }
#if SPARMAT2_COMPACTED
        delete[] rowDSS_;
        rowDSS_ = new unsigned[alc+2];
#endif
#if SPARMAT2_USES_COLNEXT
        delete[] colidx_;
        colidx_ = new unsigned[alc+2];
        for ( unsigned n = 0; n <= alc; ++n )
            colidx_[n] = n;
#endif
    }
}


index_t SparMatSym2::allocated() const
{
    index_t res = 0;
    if ( column_ )
        for ( index_t i = 0; i < alloc_; ++i )
            res += colmax_[i];
    return res;
}


void SparMatSym2::deallocate()
{
    if ( column_ )
    {
        for ( index_t i = 0; i < alloc_; ++i )
            delete[] column_[i];
        delete[] column_; column_ = nullptr;
        delete[] colsiz_; colsiz_ = nullptr;
        delete[] colmax_; colmax_ = nullptr;
#if SPARMAT2_USES_COLNEXT
        delete[] colidx_; colidx_ = nullptr;
#endif
#if SPARMAT2_COMPACTED
        delete[] colDSS_; colDSS_ = nullptr;
        delete[] rowDSS_; rowDSS_ = nullptr;
        free_real(valDSS_); valDSS_ = nullptr;
#endif
    }
    alloc_ = 0;
}


/// copy `cnt` elements from `src` to `dst`
static void copy(index_t cnt, SparMatSym2::Element * src, SparMatSym2::Element * dst)
{
    for ( index_t ii = 0; ii < cnt; ++ii )
        dst[ii] = src[ii];
}


void SparMatSym2::allocateColumn(const index_t j, index_t alc)
{
    assert_true( j < size_ );
    if ( alc > colmax_[j] )
    {
        constexpr index_t chunk = 8;
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );
        //fprintf(stderr, "SMS2 column %lu: %u --> %u\n", j, colmax_[j], alc);

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


/**
This will allocate for 5 additional values in this column,
corresponding to the maximal elements in a column from Meca::addLink4()
In this way we minimize the number of allocations in `operator()`
*/
real& SparMatSym2::diagonal(index_t i)
{
    assert_true( i < size_ );
    
    if ( 8 + colsiz_[i] > colmax_[i] )
        allocateColumn(i, 8 + colsiz_[i]);

    Element * col = column_[i];
    if ( colsiz_[i] == 0 )
    {
        col->reset(i);
        colsiz_[i] = 1;
    }

    assert_true( col->inx == i );
    return col->val;
}

/**
 This should work for any value of (i, j)
 This allocates to be able to hold the matrix element if necessary
 The diagonal element is always first in each column
*/
real& SparMatSym2::element(index_t ii, index_t jj)
{
    assert_true( ii >= jj );
    assert_true( jj < size_ );
    //fprintf(stderr, "SMS2( %6i %6i )\n", i, j);

    assert_true( colsiz_[jj] <= colmax_[jj] );
    // make space always for two additional elements
    if ( 2 + colsiz_[jj] > colmax_[jj] )
        allocateColumn(jj, 2 + colsiz_[jj]);
    
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
        if ( fence == col )
        {
            // the column is empty, insert diagonal term first:
            e->reset(jj);
            e += ( ii != jj );
            e->reset(ii);
            colsiz_[jj] = 1 + ( ii != jj );
            //printColumn(std::clog, jj);
            return e->val;
        }
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


real* SparMatSym2::address(index_t i, index_t j) const
{
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

void SparMatSym2::reset()
{
    for ( index_t jj = 0; jj < size_; ++jj )
        colsiz_[jj] = 0;
}


bool SparMatSym2::notZero() const
{
    //check for any non-zero sparse term:
    for ( index_t jj = 0; jj < size_; ++jj )
        for ( index_t kk = 0 ; kk < colsiz_[jj] ; ++kk )
            if ( column_[jj][kk].val != 0 )
                return true;
    
    //if here, the matrix is empty
    return false;
}


void SparMatSym2::scale(const real alpha)
{
    for ( index_t jj = 0; jj < size_; ++jj )
        for ( index_t n = 0; n < colsiz_[jj]; ++n )
            column_[jj][n].val *= alpha;
}


void SparMatSym2::addDiagonalBlock(real* mat, const index_t ldd, index_t start, index_t cnt,
                                   const index_t mul) const
{
    assert_true( start + cnt <= size_ );

    for ( index_t j = 0; j < cnt; ++j )
    {
        Element * col = column_[j+start];
        real * dst = mat + ( j + ldd * j ) * mul;
        for ( index_t n = 0; n < colsiz_[j+start]; ++n )
        {
            // assuming lower triangle is stored:
            assert_true( col[n].inx >= j + start );
            index_t ij = col[n].inx - ( j + start );
            if ( j + ij < cnt )
            {
                //printf("SMS2 %4i %4i % .4f\n", j, ij, a);
                dst[mul*ij] += col[n].val;
                if ( ij )
                    dst[mul*ldd*ij] += col[n].val;
            }
        }
    }
}


void SparMatSym2::addLowerBand(real alpha, real* mat, const index_t ldd, index_t start, index_t cnt,
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


void SparMatSym2::addDiagonalTrace(real alpha, real* mat, const index_t ldd,
                                   const index_t start, const index_t cnt,
                                   const index_t mul, const index_t rank, const bool sym) const
{
    ABORT_NOW("unfinished SparMatSym2::addDiagonalTrace()");
}


int SparMatSym2::bad() const
{
    if ( size_ <= 0 ) return 1;
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        for ( index_t kk = 0 ; kk < colsiz_[jj] ; ++kk )
        {
            index_t inx = column_[jj][kk].inx;
            if ( inx >= size_ ) return 2;
            if ( inx <= jj ) return 3;
        }
    }
    return 0;
}


size_t SparMatSym2::nbElements(index_t start, index_t stop, size_t& alc) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_);
#if SPARMAT2_COMPACTED
    alc = alcDSS_;
#else
    alc = 0;
#endif
    size_t cnt = 0;
    for ( index_t j = start; j < stop; ++j )
    {
        alc += colmax_[j];
        cnt += colsiz_[j];
    }
    return cnt;
}


size_t SparMatSym2::nbDiagonalElements(index_t start, index_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_);
    size_t cnt = 0;
    for ( index_t jj = start; jj < stop; ++jj )
        cnt += ( colsiz_[jj] > 0 ) && ( column_[jj][0].val != 0.0 );
    return cnt;
}


std::string SparMatSym2::what() const
{
    size_t alc;
    size_t cnt = nbElements(0, size_, alc);
    std::ostringstream msg;
#if SPARMAT2_COMPACTED
    msg << "c";
#endif
#if SPARMAT2_USES_AVX
    msg << "SMS2x ";
#elif SPARMAT2_USES_SSE
    msg << "SMS2e ";
#else
    msg << "SMS2 ";
#endif
    msg << cnt << " (" << alc << ")";
    return msg.str();
}


void SparMatSym2::printSparse(std::ostream& os, real inf, index_t start, index_t stop) const
{
    stop = std::min(stop, size_);
    os << "% SparMatSym2 size " << size_ << ":\n";
    for ( index_t jj = start; jj < stop; ++jj )
    {
        if ( colsiz_[jj] > 0 )
        {
            os << "% column " << jj << "\n";
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


void SparMatSym2::printSummary(std::ostream& os, index_t start, index_t stop)
{
    stop = std::min(stop, size_);
    os << "SparMatSym2 size " << size_ << ":";
    for ( index_t jj = start; jj < stop; ++jj )
    {
        os << "\n   " << jj << "   " << colsiz_[jj];
#if SPARMAT2_USES_COLNEXT
        os << " index " << colidx_[jj];
#endif
    }
    std::endl(os);
}


void SparMatSym2::printColumn(std::ostream& os, const index_t jj)
{
    std::streamsize p = os.precision(1);
    Element const* col = column_[jj];
    os << "SparMatSym2 col " << jj << ":";
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


void SparMatSym2::printSparseArray(std::ostream& os) const
{
#if SPARMAT2_COMPACTED
    if ( rowDSS_ )
    {
        std::ios::fmtflags fgs = os.flags();
        std::streamsize p = os.precision(2);
        
        index_t cnt = rowDSS_[size_];
        os << "\n% SparMatSym2 size " << size_ << " DSS storage:";
        os << "\nvalues   ";
        for ( index_t i = 0; i < cnt; ++i )
            os << " " << std::setw(5) << valDSS_[i];
        
        os << "\ncolumns  ";
        for ( index_t i = 0; i < cnt; ++i )
            os << " " << std::setw(5) << colDSS_[i];
        
        os << "\nrowIndex ";
        for ( index_t i = 0; i <= size_; ++i )
            os << " " << std::setw(5) << rowDSS_[i];
        
        os.precision(p);
        os.setf(fgs);
        std::endl(os);
        return;
    }
#endif
    os << "SMS2 has no compact storage\n";
}


//------------------------------------------------------------------------------
#pragma mark - Column-Vector Multiplication

/**
Multiply by column `jj` provided in `col` of size `cnt`
Attention: this assumes that col[0] is the diagonal element
*/
void SparMatSym2::vecMulAddCol(const real* X, real* Y, Element col[], index_t cnt)
{
    assert_true( cnt > 0 );
    const auto jj = col[0].inx;
    const real X0 = X[jj];
    real Y0 = Y[jj] + col[0].val * X0;
    for ( index_t n = 1; n < cnt; ++n )
    {
        const auto ii = col[n].inx;
        const real a = col[n].val;
        assert_true( ii > jj );
        Y[ii] += a * X0;
        Y0 += a * X[ii];
    }
    Y[jj] = Y0;
}

/**
 Multiply by column `jj` provided in `col` of size `cnt`
 Attention: this assumes that col[0] is the diagonal element
 */
void SparMatSym2::vecMulAddColIso2D(const real* X, real* Y, Element col[], index_t cnt)
{
    assert_true( cnt > 0 );
    const auto jj = 2 * col[0].inx;
    const real X0 = X[jj  ];
    const real X1 = X[jj+1];
    real Y0 = Y[jj  ] + col[0].val * X0;
    real Y1 = Y[jj+1] + col[0].val * X1;
    for ( index_t n = 1; n < cnt; ++n )
    {
        const auto ii = 2 * col[n].inx;
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
Multiply by column `jj` provided in `col` of size `cnt`
Attention: this assumes that col[0] is the diagonal element
*/
void SparMatSym2::vecMulAddColIso3D(const real* X, real* Y, Element col[], index_t cnt)
{
    assert_true( cnt > 0 );
    const auto jj = 3 * col[0].inx;
    const real X0 = X[jj  ];
    const real X1 = X[jj+1];
    const real X2 = X[jj+2];
    real Y0 = Y[jj  ] + col[0].val * X0;
    real Y1 = Y[jj+1] + col[0].val * X1;
    real Y2 = Y[jj+2] + col[0].val * X2;
    for ( index_t n = 1; n < cnt; ++n )
    {
        const auto ii = 3 * col[n].inx;
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

void SparMatSym2::setColumnIndex()
{
#if SPARMAT2_USES_COLNEXT
    colidx_[size_] = size_;
    if ( size_ > 0 )
    {
        index_t inx = size_;
        index_t nxt = size_;
        while ( inx-- > 0 )
        {
            if ( colsiz_[inx] > 0 )
                nxt = inx;
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
        printColumn(std::clog, j);
    }
    std::clog << "SMS2 has " << cnt << " / " << size_ << " empty columns\n";
#endif
}

#if !SPARMAT2_COMPACTED

bool SparMatSym2::prepareForMultiply(int)
{
    setColumnIndex();
    return true;
}

#else

/**
 This builds the lower half DSS Symmetric Matrix Storage.
 This is also called the Compact Column Storage
 In this class the matrix is symmetric and storing the upper half would be equivalent.
 
 The storage uses three arrays: values, columns, and rows
 
 * valDSS_[] contains the non-zero elements of the sparse matrix.
 * colDSS_[i] is the index of the column corresponding to valDSS_[i].
 * rowDSS_[j] gives the index of the first non-zero element in row j.
 
 The length of valDSS_[] and colDSS_[i] is equal to the number of non-zero elements in the matrix.
 The length of rowDSS_[] is the number of rows in the matrix plus one.
 
 As rowDSS_[] gives the location of the first non-zero element within a row,
 and the non-zero elements are stored consecutively, the number of non-zero
 elements in the i-th row is equal to ( rowDSS_[i+1] - rowDSS_[i] ).
 
 To have this relationship hold for the last row of the matrix, an additional entry
 is added to the end of rowDSS_, equal to the number of non-zero elements plus one.
 Thus the size of rowDSS_[] is one plus the number of row in the matrix, and
 rowDSS_[size] = number_of_non_zero_elements
 */
bool SparMatSym2::prepareForMultiply(int dimension)
{
    assert_true( size_ <= alloc_ );
    
    setColumnIndex();

    //count number of non-zero elements, always including the diagonal term
    index_t nnz = 0;
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        if ( colsiz_[jj] > 0 )
            nnz += colsiz_[jj];
        else
            ++nnz;
    }
    
    // allocate DSS sparse matrix storage
    if ( nnz > alcDSS_ )
    {
        constexpr index_t chunk = 16;
        alcDSS_ = ( nnz + chunk - 1 ) & ~( chunk -1 );
        delete[] colDSS_;
        free_real(valDSS_);
        colDSS_ = new unsigned[alcDSS_];
        valDSS_ = new_real(alcDSS_);
    }
    
    /*
     Create the DSS sparse matrix storage.
     */
    index_t cnt = 0;
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        rowDSS_[jj] = cnt;
        // always put diagonal term:
        valDSS_[cnt] = 0.0;
        colDSS_[cnt] = jj * dimension;
        real & dia = valDSS_[cnt];
        ++cnt;
        if ( colsiz_[jj] > 0 )
        {
            Element * col = column_[jj];
            assert_true( col[0].inx == jj );
            dia = col[0].val;
            for ( index_t k = 1; k < colsiz_[jj]; ++k )
            if ( col[k].val != 0 )
            {
                assert_true( cnt < alcDSS_ );
                valDSS_[cnt] = col[k].val;
                colDSS_[cnt] = col[k].inx * dimension;
                ++cnt;
            }
        }
    }
    if ( cnt > nnz )
        ABORT_NOW("internal error");
    if ( size_ > 0 )
        rowDSS_[size_] = cnt;
    
    //printSparse(std::clog, 0, 0, size_);
    //printSparseArray(std::clog);
    return true;
}

//------------------------------------------------------------------------------
#pragma mark - Optimized Column-Vector multiplication


void SparMatSym2::vecMulAddCol(const real* X, real* Y,
                               index_t start, index_t stop) const
{
    assert_true( start < stop );
    assert_true( stop <= alcDSS_ );
    index_t jj = colDSS_[start];
    real X0 = X[jj];
    real Y0 = Y[jj] + valDSS_[start] * X0;
    for ( index_t n = start+1; n < stop; ++n )
    {
        real a = valDSS_[n];
        index_t ii = colDSS_[n];
        Y[ii] += a * X0;
        Y0 += a * X[ii];
    }
    Y[jj] = Y0;
}

void SparMatSym2::vecMulAddColIso2D(const real* X, real* Y,
                                    index_t start, index_t stop) const
{
    assert_true( start < stop );
    assert_true( stop <= alcDSS_ );
    index_t jj = colDSS_[start];
    real X0 = X[jj  ];
    real X1 = X[jj+1];
    real Y0 = Y[jj  ] + valDSS_[start] * X0;
    real Y1 = Y[jj+1] + valDSS_[start] * X1;
    for ( index_t n = start+1; n < stop; ++n )
    {
        index_t ii = colDSS_[n];
        assert_true( ii > jj );
        real a = valDSS_[n];
        Y0 += a * X[ii  ];
        Y1 += a * X[ii+1];
        Y[ii  ] += a * X0;
        Y[ii+1] += a * X1;
    }
    Y[jj  ] = Y0;
    Y[jj+1] = Y1;
}


void SparMatSym2::vecMulAddColIso3D(const real* X, real* Y,
                                    index_t start, index_t stop) const
{
    assert_true( start < stop );
    assert_true( stop <= alcDSS_ );
    index_t jj = colDSS_[start];
    real X0 = X[jj  ];
    real X1 = X[jj+1];
    real X2 = X[jj+2];
    real Y0 = Y[jj  ] + valDSS_[start] * X0;
    real Y1 = Y[jj+1] + valDSS_[start] * X1;
    real Y2 = Y[jj+2] + valDSS_[start] * X2;
    for ( index_t n = start+1; n < stop; ++n )
    {
        index_t ii = colDSS_[n];
        assert_true( ii > jj );
        real a = valDSS_[n];
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

#if REAL_IS_DOUBLE && SPARMAT2_USES_SSE

static inline void multiply2(const double* X, double* Y, index_t ii,
                             const double* val, vec2 const& xx, vec2& ss)
{
    vec2 aa = loaddup2(val);
    ss = fmadd2(load2(X+ii), aa, ss);
    store2(Y+ii, fmadd2(xx, aa, load2(Y+ii)));
}


void SparMatSym2::vecMulAddColIso2D_SSE(const double* X, double* Y,
                                        index_t start, index_t stop) const
{
    assert_true( start < stop );
    assert_true( stop <= alcDSS_ );
    index_t jj = colDSS_[start];
    const vec2 xx = load2(X+jj);
    //process diagonal element:
    vec2 ss = fmadd2(loaddup2(valDSS_+start), xx, load2(Y+jj));
    // there is a dependence here for 'ss'
    for ( index_t n = start+1; n < stop; ++n )
        multiply2(X, Y, colDSS_[n], valDSS_+n, xx, ss);
    store2(Y+jj, ss);
}


void SparMatSym2::vecMulAddColIso2D_SSEU(const double* X, double* Y,
                                         index_t start, index_t stop) const
{
    assert_true( start < stop );
    assert_true( stop <= alcDSS_ );
    
    index_t jj = colDSS_[start];
    const vec2 xx = load2(X+jj);
    //process diagonal element:
    vec2 s0 = mul2(loaddup2(valDSS_+start), xx);
    vec2 s1 = load2(Y+jj);
    vec2 s2 = setzero2();
    vec2 s3 = setzero2();
    
    unsigned * inx = colDSS_ + start + 1;
    const double * val = valDSS_ + start + 1;
    const double * end = valDSS_ + stop;

#if ( 1 )
    // unrolling by 6 might not lead to much improvement
    const double * pause = end - 5;  // val+7 <= end-1  is  val < end-7;
    if ( val < pause )
    {
        // process 8 by 8:
        #pragma nounroll
        for ( ; val < pause; val += 6 )
        {
            const auto i0 = inx[0];
            const auto i1 = inx[1];
            const auto i2 = inx[2];
            const auto i3 = inx[3];
            const auto i4 = inx[4];
            const auto i5 = inx[5];
            inx += 6;
            vec2 y0 = load2(Y+i0);
            vec2 y1 = load2(Y+i1);
            vec2 y2 = load2(Y+i2);
            vec2 y3 = load2(Y+i3);
            vec2 a0 = loaddup2(val);
            vec2 a1 = loaddup2(val+1);
            vec2 a2 = loaddup2(val+2);
            vec2 a3 = loaddup2(val+3);
            s0 = fmadd2(load2(X+i0), a0, s0);
            s1 = fmadd2(load2(X+i1), a1, s1);
            s2 = fmadd2(load2(X+i2), a2, s2);
            s3 = fmadd2(load2(X+i3), a3, s3);
            vec2 y4 = load2(Y+i4);
            vec2 y5 = load2(Y+i5);
            vec2 a4 = loaddup2(val+4);
            vec2 a5 = loaddup2(val+5);
            store2(Y+i0, fmadd2(xx, a0, y0));
            store2(Y+i1, fmadd2(xx, a1, y1));
            store2(Y+i2, fmadd2(xx, a2, y2));
            store2(Y+i3, fmadd2(xx, a3, y3));
            s0 = fmadd2(load2(X+i4), a4, s0);
            s1 = fmadd2(load2(X+i5), a5, s1);
            store2(Y+i4, fmadd2(xx, a4, y4));
            store2(Y+i5, fmadd2(xx, a5, y5));
        }
    }
#endif
    
    // process 4 by 4:
    const double * halt = end - 3;  // val+3 <= end-1  is  val < end-3;
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
        multiply2(X, Y, inx[0], val  , xx, s0);
        multiply2(X, Y, inx[1], val+1, xx, s1);
        multiply2(X, Y, inx[2], val+2, xx, s2);
        multiply2(X, Y, inx[3], val+3, xx, s3);
#else
        /* we remove here the apparent dependency on the values of Y[],
         which are read and written, but at different indices.
         The compiler can reorder instructions to avoid lattencies */
        const auto i0 = inx[0];
        const auto i1 = inx[1];
        const auto i2 = inx[2];
        const auto i3 = inx[3];
        vec2 y0 = load2(Y+i0);
        vec2 y1 = load2(Y+i1);
        vec2 y2 = load2(Y+i2);
        vec2 y3 = load2(Y+i3);
        vec2 a0 = loaddup2(val);
        vec2 a1 = loaddup2(val+1);
        vec2 a2 = loaddup2(val+2);
        vec2 a3 = loaddup2(val+3);
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
    while ( val < end )
        multiply2(X, Y, *inx++, val++, xx, s0);
    store2(Y+jj, s0);
}

#endif

#if SPARMAT2_USES_AVX && REAL_IS_DOUBLE

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


void SparMatSym2::vecMulAddColIso2D_AVX(const double* X, double* Y,
                                        index_t start, index_t stop) const
{
    assert_true( start < stop );
    assert_true( stop <= alcDSS_ );
    index_t jj = colDSS_[start];
    const vec4 xx = broadcast2(X+jj);  // hi position
    //process diagonal element:
    vec4 ss = fmadd4(broadcast1(valDSS_+start), xx, broadcast2(Y+jj));
    // there is a dependence here for 'ss'
    for ( index_t n = start+1; n < stop; ++n )
        multiply4(X, Y, colDSS_[n], valDSS_+n, xx, ss);
    store2(Y+jj, gethi(ss));
}


void SparMatSym2::vecMulAddColIso2D_AVXU(const double* X, double* Y,
                                         index_t start, index_t stop) const
{
    assert_true( start < stop );
    assert_true( stop <= alcDSS_ );
    index_t jj = colDSS_[start];
    const vec4 xx = broadcast2(X+jj);  // hi and lo position
    //process diagonal element:
    vec4 s0 = mul4(broadcast1(valDSS_+start), xx);
    vec4 s1 = broadcast2(Y+jj);
    vec4 s2 = setzero4();
    vec4 s3 = setzero4();
    
    unsigned * inx = colDSS_ + start + 1;
    const double * val = valDSS_ + start + 1;
    const double * end = valDSS_ + stop;
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

#if SPARMAT2_USES_AVX && SPARMAT2_COMPACTED && REAL_IS_DOUBLE
void SparMatSym2::vecMulAddColIso3D_AVX(const double* X, double* Y,
                                        index_t start, index_t stop) const
{
    assert_true( start < stop );
    assert_true( stop <= alcDSS_ );
    //Attention: this assumes that jj is the diagonal element
    index_t jj = colDSS_[start];
    unsigned * inx = colDSS_ + start;
    double const* val = valDSS_ + start;
    double const* end = valDSS_ + stop;
    
    //printf("SparMatSym2 column %lu has %lu elements\n", jj, stop - start);
    vec4 zz = setzero4();
    vec4 xx = blend31(loadu4(X+jj), zz);
    vec4 yy = fmadd4(broadcast1(val), xx, loadu4(Y+jj));
    // process one element when the number of values is even
    if (( stop & 1 ) == ( start & 1 ))
    {
        index_t ii = *(++inx);
        vec4 aa = broadcast1(++val);
        zz = mul4(aa, loadu4(X+ii));
        storeu4(Y+ii, fmadd4(aa, xx, loadu4(Y+ii)));
    }
    ++val;
    ++inx;
    while ( val < end )
    {
        index_t ii = inx[0];
        index_t kk = inx[1];
        assert_true( kk > ii );
        inx += 2;
        vec4 aa = broadcast1(val);
        vec4 bb = broadcast1(val+1);
        vec4 nn = loadu4(Y+ii);
        vec4 mm = loadu4(Y+kk);
        val += 2;
        yy = fmadd4(aa, loadu4(X+ii), yy);
        zz = fmadd4(bb, loadu4(X+kk), zz);
        storeu4(Y+ii, fmadd4(aa, xx, nn));
        storeu4(Y+kk, fmadd4(bb, xx, mm));
    }
    store3(Y+jj, add4(yy, zz));
    assert_true( val == end );
}
#endif


#if SPARMAT2_USES_SSE && SPARMAT2_COMPACTED && !REAL_IS_DOUBLE
void SparMatSym2::vecMulAddColIso3D_SSE(const float* X, float* Y,
                                        index_t start, index_t stop) const
{
    assert_true( start < stop );
    assert_true( stop <= alcDSS_ );
    //Attention: this assumes that jj is the diagonal element
    index_t jj = colDSS_[start];
    unsigned * inx = colDSS_ + start;
    float const* val = valDSS_ + start;
    float const* end = valDSS_ + stop;
    
    //printf("SparMatSym2 column %lu has %lu elements\n", jj, stop - start);
    vec4f s0 = setzero4f();
    vec4f xx = blend31f(loadu4f(X+jj), s0);
    vec4f s1 = fmadd4f(broadcast1f(val), xx, loadu4f(Y+jj));
    // process one element when the number of values is even
    if (( stop & 1 ) == ( start & 1 ))
    {
        index_t ii = *(++inx);
        vec4f aa = broadcast1f(++val);
        s0 = mul4f(aa, loadu4f(X+ii));
        storeu4f(Y+ii, fmadd4f(aa, xx, loadu4f(Y+ii)));
    }
    ++val;
    ++inx;
#if ( 0 )
    vec4f s2 = setzero4f();
    vec4f s3 = setzero4f();
    float const* halt = end - 2;
    #pragma ivdep unroll (4)
    #pragma clang loop unroll(disable)
    while ( val < halt )
    {
        index_t i0 = inx[0];
        index_t i1 = inx[1];
        index_t i2 = *(inx+2);
        index_t i3 = *(inx+3);
        assert_true( i1 > i0 );
        inx += 4;
#if 0
        vec4f a0 = broadcast1f(val);
        vec4f a1 = broadcast1f(val+1);
        vec4f a2 = broadcast1f(val+2);
        vec4f a3 = broadcast1f(val+3);
#else
        vec4f a3 = loadu4f(val);
        vec4f a0 = broadcastXf(a3);
        vec4f a1 = broadcastYf(a3);
        vec4f a2 = broadcastZf(a3);
        a3 = broadcastTf(a3);
#endif
        val += 4;
        s0 = fmadd4f(a0, loadu4f(X+i0), s0);
        s1 = fmadd4f(a1, loadu4f(X+i1), s1);
        s2 = fmadd4f(a2, loadu4f(X+i2), s2);
        s3 = fmadd4f(a3, loadu4f(X+i3), s3);
        storeu4f(Y+i0, fmadd4f(a0, xx, loadu4f(Y+i0)));
        storeu4f(Y+i1, fmadd4f(a1, xx, loadu4f(Y+i1)));
        storeu4f(Y+i2, fmadd4f(a2, xx, loadu4f(Y+i2)));
        storeu4f(Y+i3, fmadd4f(a3, xx, loadu4f(Y+i3)));
    }
    s0 = add4f(s0, s2);
    s1 = add4f(s1, s3);
#endif
    while ( val < end )
    {
        index_t i0 = inx[0];
        index_t i1 = inx[1];
        assert_true( i1 > i0 );
        inx += 2;
#if 1
        vec4f a0 = broadcast1f(val);
        vec4f a1 = broadcast1f(val+1);
#else
        vec4f a1 = load2f(val);
        vec4f a0 = broadcastXf(a1);
        a1 = broadcastYf(a1);
#endif
        val += 2;
        s0 = fmadd4f(a0, loadu4f(X+i0), s0);
        s1 = fmadd4f(a1, loadu4f(X+i1), s1);
        storeu4f(Y+i0, fmadd4f(a0, xx, loadu4f(Y+i0)));
        storeu4f(Y+i1, fmadd4f(a1, xx, loadu4f(Y+i1)));
    }
    store3f(Y+jj, add4f(s1, s0));
    assert_true( val == end );
}
#endif

//------------------------------------------------------------------------------
#pragma mark - Matrix-Vector Add-multiply


void SparMatSym2::vecMulAdd_ALT(const real* X, real* Y, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_);

    for ( index_t jj = start; jj < stop; ++jj )
        if ( colsiz_[jj] > 0 )
        {
            assert_true(column_[jj][0].inx == jj);
            vecMulAddCol(X, Y, column_[jj], colsiz_[jj]);
        }
}


void SparMatSym2::vecMulAdd(const real* X, real* Y, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_);

#if SPARMAT2_USES_COLNEXT
    for ( index_t jj = colidx_[start]; jj < stop; jj = colidx_[jj+1] )
#else
    for ( index_t jj = start; jj < stop; ++jj )
    if ( colsiz_[jj] > 0 )
#endif
    {
#if SPARMAT2_COMPACTED
        // check this is diagonal term:
        assert_true(colDSS_[rowDSS_[jj]] == jj);
        vecMulAddCol(X, Y, rowDSS_[jj], rowDSS_[jj+1]);
#else
        assert_true(column_[jj][0].inx == jj);
        vecMulAddCol(X, Y, column_[jj], colsiz_[jj]);
#endif
    }
}


void SparMatSym2::vecMulAddIso2D(const real* X, real* Y, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_);

#if SPARMAT2_USES_COLNEXT
    for ( index_t jj = colidx_[start]; jj < stop; jj = colidx_[jj+1] )
#else
    for ( index_t jj = start; jj < stop; ++jj )
    if ( colsiz_[jj] > 0 )
#endif
    {
#if SPARMAT2_COMPACTED
        assert_true(colDSS_[rowDSS_[jj]] == 2*jj);
#  if SPARMAT2_USES_AVX && REAL_IS_DOUBLE
        //vecMulAddColIso2D_AVXU(X, Y, rowDSS_[jj], rowDSS_[jj+1]);
        vecMulAddColIso2D_SSEU(X, Y, rowDSS_[jj], rowDSS_[jj+1]);
#  else
        vecMulAddColIso2D(X, Y, rowDSS_[jj], rowDSS_[jj+1]);
#  endif
#else
        assert_true(column_[jj][0].inx == jj);
        vecMulAddColIso2D(X, Y, column_[jj], colsiz_[jj]);
#endif

    }
}


void SparMatSym2::vecMulAddIso3D(const real* X, real* Y, index_t start, index_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_);

#if SPARMAT2_USES_COLNEXT
    for ( index_t jj = colidx_[start]; jj < stop; jj = colidx_[jj+1] )
#else
    for ( index_t jj = start; jj < stop; ++jj )
    if ( colsiz_[jj] > 0 )
#endif
    {
#if SPARMAT2_COMPACTED
        assert_true(colDSS_[rowDSS_[jj]] == 3*jj);
#  if SPARMAT2_USES_AVX && REAL_IS_DOUBLE
        vecMulAddColIso3D_AVX(X, Y, rowDSS_[jj], rowDSS_[jj+1]);
#  elif SPARMAT2_USES_SSE && !REAL_IS_DOUBLE
        vecMulAddColIso3D_SSE(X, Y, rowDSS_[jj], rowDSS_[jj+1]);
#  else
        vecMulAddColIso3D(X, Y, rowDSS_[jj], rowDSS_[jj+1]);
#  endif
#else
        assert_true(column_[jj][0].inx == jj);
        vecMulAddColIso3D(X, Y, column_[jj], colsiz_[jj]);
#endif
    }
}


//------------------------------------------------------------------------------
#pragma mark - Matrix-Vector multiplication

void SparMatSym2::vecMul(const real* X, real* Y, index_t start, index_t stop) const
{
    zero_real(stop-start, Y+start);
    vecMulAdd(X, Y, start, stop);
}
