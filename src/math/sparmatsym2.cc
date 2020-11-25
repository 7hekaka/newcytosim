// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#include <cmath>
#include "sparmatsym2.h"
#include "assert_macro.h"
#include "blas.h"

#include <iomanip>
#include <sstream>
#include <iostream>


#ifdef __AVX__
#  define MATRIX2_USES_AVX REAL_IS_DOUBLE
#  define MATRIX2_USES_SSE 0
#  include "simd.h"
#elif defined(__SSE3__)
#  define MATRIX2_USES_AVX REAL_IS_DOUBLE
#  define MATRIX2_USES_SSE REAL_IS_DOUBLE
#  include "simd.h"
#else
#  define MATRIX2_USES_AVX 0
#  define MATRIX2_USES_SSE 0
#endif


SparMatSym2::SparMatSym2()
{
    size_   = 0;
    alloc_  = 0;
    column_ = nullptr;
    colsiz_ = nullptr;
    colmax_ = nullptr;
#if MATRIX2_OPTIMIZED_MULTIPLY
    alcDSS_ = 0;
    colDSS_ = nullptr;
    rowDSS_ = nullptr;
    valDSS_ = nullptr;
#endif
#if MATRIX2_USES_COLNEXT
    next_ = new unsigned[1];
    next_[0] = 0;
#endif
}


void SparMatSym2::allocate(size_t alc)
{
    if ( alc > alloc_ )
    {
        /*
         'chunk' can be increased to gain performance:
         more memory will be used, but reallocation will be less frequent
         */
        constexpr size_t chunk = 64;
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );

        //fprintf(stderr, "SMS2 allocate matrix %u\n", alc);
        Element ** column_new = new Element*[alc];
        unsigned * colsiz_new = new unsigned[alc];
        unsigned * colmax_new = new unsigned[alc];
        
        size_t ii = 0;
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
#if MATRIX2_OPTIMIZED_MULTIPLY
        delete[] rowDSS_;
        rowDSS_ = new unsigned[alc+1];
#endif
#if MATRIX2_USES_COLNEXT
        delete[] next_;
        next_ = new unsigned[alc+1];
        for ( size_t n = 0; n <= alc; ++n )
            next_[n] = n;
#endif
    }
}


void SparMatSym2::deallocate()
{
    if ( column_ )
    {
        for ( size_t ii = 0; ii < alloc_; ++ii )
            delete[] column_[ii];
        delete[] column_; column_ = nullptr;
        delete[] colsiz_; colsiz_ = nullptr;
        delete[] colmax_; colmax_ = nullptr;
#if MATRIX2_OPTIMIZED_MULTIPLY
        delete[] colDSS_; colDSS_ = nullptr;
        delete[] rowDSS_; rowDSS_ = nullptr;
        free_real(valDSS_); valDSS_ = nullptr;
#endif
#if MATRIX2_USES_COLNEXT
        delete[] next_;  next_ = nullptr;
#endif
    }
    alloc_ = 0;
}


/// copy `cnt` elements from `src` to `dst`
void copy(size_t cnt, SparMatSym2::Element * src, SparMatSym2::Element * dst)
{
    for ( size_t ii = 0; ii < cnt; ++ii )
        dst[ii] = src[ii];
}

/// move `cnt` elements to the next index, starting at vec[0]
void shift(size_t cnt, SparMatSym2::Element * vec)
{
    for ( size_t ii = cnt; ii > 0; --ii )
        vec[ii] = vec[ii-1];
}


void SparMatSym2::allocateColumn(const size_t jj, size_t alc)
{
    assert_true( jj < size_ );
    if ( alc > colmax_[jj] )
    {
        //fprintf(stderr, "SMS2 allocate column %i size %u\n", jj, alc);
        constexpr size_t chunk = 16;
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );
        Element * ptr = new Element[alc];
        
        if ( column_[jj] )
        {
            //copy over previous column elements
            copy(colsiz_[jj], column_[jj], ptr);
            
            //release old memory
            delete[] column_[jj];
        }
        column_[jj] = ptr;
        colmax_[jj] = alc;
    }
}


SparMatSym2::Element * SparMatSym2::insertElement(const size_t jj, size_t inx)
{
    assert_true( jj < size_ );
    // allocate space for new Element if necessary:
    if ( colsiz_[jj] >= colmax_[jj] )
    {
        constexpr size_t chunk = 16;
        size_t alc = ( colsiz_[jj] + chunk ) & ~( chunk -1 );
        Element * ptr = new Element[alc];
        if ( column_[jj] )
        {
            copy(inx, column_[jj], ptr);
            copy(colsiz_[jj]-inx, column_[jj]+inx, ptr+inx+1);
            delete[] column_[jj];
        }
        column_[jj] = ptr;
        colmax_[jj] = alc;
    }
    else
    {
        shift(colsiz_[jj]-inx, column_[jj]+inx);
    }
    column_[jj][inx].reset(0);
    ++colsiz_[jj];
    return column_[jj]+inx;
}


real& SparMatSym2::diagonal(size_t i)
{
    assert_true( i < size_ );
    
    Element * col;
    
    if ( colsiz_[i] == 0 )
    {
        allocateColumn(i, 1);
        col = column_[i];
        //diagonal term always first:
        col->reset(i);
        colsiz_[i] = 1;
    }
    else
    {
        col = column_[i];
        assert_true( col->inx == i );
    }
    
    return col->val;
}

/**
 This allocates to be able to hold the matrix element if necessary
*/
real& SparMatSym2::operator()(size_t i, size_t j)
{
    assert_true( i < size_ );
    assert_true( j < size_ );
    //fprintf(stderr, "SMS2( %6i %6i )\n", i, j);
        
    Element * col;
    
    if ( i == j )
    {
        // return diagonal element
        if ( colsiz_[j] <= 0 )
        {
            allocateColumn(j, 1);
            col = column_[j];
            // put diagonal term always first:
            col->reset(j);
            colsiz_[j] = 1;
        }
        else
        {
            col = column_[j];
            assert_true( col->inx == j );
        }
        return col->val;
    }
 
    // swap to get ii > jj (address lower triangle)
    size_t ii = std::max(i, j);
    size_t jj = std::min(i, j);

    //check if the column is empty:
    if ( colsiz_[jj] < 2 )
    {
        allocateColumn(jj, 2);
        col = column_[jj];
        if ( colsiz_[jj] == 0 )
        {
            // put diagonal term always first:
            col->reset(jj);
        }
        //add the requested term:
        col[1].reset(ii);
        colsiz_[jj] = 2;
        return col[1].val;
    }
    
    col = column_[jj];
    Element * e = col + 1;
    Element * lst = col + colsiz_[jj] - 1;
    
    //search, knowing that elements are kept ordered in the column:
    while ( e->inx < ii )
    {
        if ( ++e > lst )
        {
            // add one element last
            size_t n = colsiz_[jj];
            if ( n >= colmax_[jj] )
            {
                allocateColumn(jj, n+1);
                col = column_[jj];
            }
            ++colsiz_[jj];
            col[n].reset(ii);
            return col[n].val;
        }
    }
    
    if ( e->inx == ii )
        return e->val;
    
    size_t n = (size_t)( e - col );

    assert_true( col[n].inx > ii );
    col = insertElement(jj, n);
    assert_true( n < colmax_[jj] );

    // add the requested term
    col->reset(ii);

    //printColumn(std::clog, jj);
    return col->val;
}


real* SparMatSym2::addr(size_t i, size_t j) const
{
    // swap to get ii <= jj (address lower triangle)
    size_t ii = std::max(i, j);
    size_t jj = std::min(i, j);

    for ( size_t kk = 0; kk < colsiz_[jj]; ++kk )
        if ( column_[jj][kk].inx == ii )
            return &( column_[jj][kk].val );
    
    return nullptr;
}


//------------------------------------------------------------------------------
#pragma mark -

void SparMatSym2::reset()
{
    for ( size_t jj = 0; jj < size_; ++jj )
        colsiz_[jj] = 0;
}


bool SparMatSym2::isNotZero() const
{
    //check for any non-zero sparse term:
    for ( size_t jj = 0; jj < size_; ++jj )
        for ( size_t kk = 0 ; kk < colsiz_[jj] ; ++kk )
            if ( column_[jj][kk].val != 0 )
                return true;
    
    //if here, the matrix is empty
    return false;
}


void SparMatSym2::scale(const real alpha)
{
    for ( size_t jj = 0; jj < size_; ++jj )
        for ( size_t n = 0; n < colsiz_[jj]; ++n )
            column_[jj][n].val *= alpha;
}


void SparMatSym2::addDiagonalBlock(real* mat, const size_t ldd,
                                   const size_t start, const size_t cnt,
                                   const size_t amp) const
{
    size_t end = start + cnt;
    assert_true( end <= size_ );
    
    for ( size_t jj = start; jj < end; ++jj )
    {
        size_t j = amp * ( jj - start );
        for ( size_t n = 0; n < colsiz_[jj]; ++n )
        {
            size_t ii = column_[jj][n].inx;
            // assuming lower triangle is stored:
            if ( ii < end )
            {
                size_t i = amp * ( ii - start );
                //printf("SMS2 %4i %4i % .4f\n", ii, jj, a);
                mat[i+ldd*j] += column_[jj][n].val;
                if ( j != i )
                    mat[j+ldd*i] += column_[jj][n].val;
            }
        }
    }
}


/*
addresses `mat' using lower banded storage for a symmetric matrix
mat(i, j) is stored in mat[i-j+ldd*j]
*/
void SparMatSym2::addTriangularBlockBanded(real alpha, real* mat, const size_t ldd,
                                           const size_t start, const size_t cnt,
                                           const size_t rank) const
{
    size_t end = start + cnt;
    assert_true( end <= size_ );
    
    for ( size_t jj = start; jj < end; ++jj )
    {
        size_t j = jj - start;
        for ( size_t n = 0; n < colsiz_[jj]; ++n )
        {
            size_t ii = column_[jj][n].inx;
            // assuming lower triangle is stored:
            if ( ii < end )
            {
                size_t i = ii - start;
                //printf("SMS2 %4i %4i % .4f\n", ii, jj, a);
                assert_true( i > j );
                // with banded storage, mat(i, j) is stored in mat[i-j+ldd*j]
                if ( i <= j + rank )
                    mat[i-j+ldd*j] += alpha * column_[jj][n].val;
            }
        }
    }
}



void SparMatSym2::addDiagonalTrace(real alpha, real* mat, const size_t ldd,
                                   const size_t start, const size_t cnt) const
{
    fprintf(stderr, "unfinished SparMatSym2::addDiagonalTrace()\n");
    exit(1);
}


void SparMatSym2::addDiagonalTraceBanded(real alpha, real* mat, const size_t ldd,
                                         const size_t start, const size_t cnt,
                                         const size_t rank) const
{
    fprintf(stderr, "unfinished SparMatSym2::addDiagonalTraceBanded()\n");
    exit(1);
}


int SparMatSym2::bad() const
{
    if ( size_ <= 0 ) return 1;
    for ( size_t jj = 0; jj < size_; ++jj )
    {
        for ( size_t kk = 0 ; kk < colsiz_[jj] ; ++kk )
        {
            if ( column_[jj][kk].inx >= size_ ) return 2;
            if ( column_[jj][kk].inx <= jj )   return 3;
        }
    }
    return 0;
}


size_t SparMatSym2::nbElements(size_t start, size_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= size_ );
    //all allocated elements are counted, even if zero
    size_t cnt = 0;
    for ( size_t jj = start; jj < stop; ++jj )
        cnt += colsiz_[jj];
    return cnt;
}


size_t SparMatSym2::nbDiagonalElements(size_t start, size_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= size_ );
    //all allocated elements are counted, even if zero
    size_t cnt = 0;
    for ( size_t jj = start; jj < stop; ++jj )
        cnt += ( colsiz_[jj] > 0 ) && ( column_[jj][0].val != 0.0 );
    return cnt;
}


std::string SparMatSym2::what() const
{
    std::ostringstream msg;
#if MATRIX2_USES_AVX
    msg << "SMS2x " << nbElements();
#elif MATRIX2_USES_SSE
    msg << "SMS2e " << nbElements();
#else
    msg << "SMS2 " << nbElements();
#endif
    return msg.str();
}


void SparMatSym2::printSparse(std::ostream& os, real inf, size_t start, size_t stop) const
{
    stop = std::min(stop, size_);
    char str[256];
    os << "\n% SparseMatSym2 size " << size_;
    for ( size_t jj = start; jj < stop; ++jj )
    {
        if ( colsiz_[jj] > 0 )
        {
            os << "\n% column " << jj;
            for ( size_t n = 0 ; n < colsiz_[jj] ; ++n )
            {
                real v = column_[jj][n].val;
                if ( abs_real(v) >= inf )
                {
                    snprintf(str, sizeof(str), "\n%6u %6lu %16.6f", column_[jj][n].inx, jj, v);
                    os << str;
                }
            }
        }
    }
    std::endl(os);
}


void SparMatSym2::printColumns(std::ostream& os, size_t start, size_t stop)
{
    stop = std::min(stop, size_);
    os << "SMS2 size " << size_ << ":";
    for ( size_t jj = start; jj < stop; ++jj )
    {
        os << "\n   " << jj << "   " << colsiz_[jj];
#if MATRIX2_USES_COLNEXT
        os << " next " << next_[jj];
#endif
    }
    std::endl(os);
}


void SparMatSym2::printColumn(std::ostream& os, const size_t jj)
{
    Element const* col = column_[jj];
    os << "SMS2 col " << jj << ":";
    for ( size_t n = 0; n < colsiz_[jj]; ++n )
    {
        os << "\n" << col[n].inx << " :";
        os << " " << col[n].val;
    }
    std::endl(os);
}


void SparMatSym2::printSparseArray(std::ostream& os) const
{
#if MATRIX2_OPTIMIZED_MULTIPLY
    std::ios::fmtflags fgs = os.flags();
    std::streamsize p = os.precision();

    size_t cnt = rowDSS_[size_];
    os << "\n% SparseMatSym2 size " << size_ << " DSS storage:";
    os << "\nvalues   ";
    os.precision(2);
    for ( size_t i = 0; i < cnt; ++i )
        os << " " << std::setw(5) << valDSS_[i];
    
    os << "\ncolumns  ";
    for ( size_t i = 0; i < cnt; ++i )
        os << " " << std::setw(5) << colDSS_[i];
                
    os.precision(2);
    os << "\nrowIndex ";
    for ( size_t i = 0; i <= size_; ++i )
        os << " " << std::setw(5) << rowDSS_[i];

    os.precision(p);
    os.setf(fgs);
    std::endl(os);
#else
    os << "no alternative sparse matrix storage\n";
#endif
}


//------------------------------------------------------------------------------
#pragma mark - Column-Vector Multiplication

/**
Multiply by column `jj` provided in `col` of size `cnt`
*/
void SparMatSym2::vecMulAddCol(const real* X, real* Y, size_t jj, Element col[], size_t cnt) const
{
    assert_true( cnt > 0 );
    const real X0 = X[jj];
    real Y0 = Y[jj] + col[0].val * X0;
    for ( size_t n = 1 ; n < cnt ; ++n )
    {
        const size_t ii = col[n].inx;
        const real a = col[n].val;
        Y[ii] += a * X0;
        assert_true( ii > jj );
        Y0 += a * X[ii];
    }
    Y[jj] = Y0;
}

/**
 Multiply by column `jj` provided in `col` of size `cnt`
 */
void SparMatSym2::vecMulAddColIso2D(const real* X, real* Y, size_t jj, Element col[], size_t cnt) const
{
    assert_true( cnt > 0 );
    const real X0 = X[jj  ];
    const real X1 = X[jj+1];
    real Y0 = Y[jj  ] + col[0].val * X0;
    real Y1 = Y[jj+1] + col[0].val * X1;
    for ( size_t n = 1 ; n < cnt ; ++n )
    {
        const size_t ii = 2 * col[n].inx;
        const real a = col[n].val;
        Y[ii  ] += a * X0;
        Y[ii+1] += a * X1;
        assert_true( ii > jj );
        Y0 += a * X[ii  ];
        Y1 += a * X[ii+1];
    }
    Y[jj  ] = Y0;
    Y[jj+1] = Y1;
}

/**
Multiply by column `jj` provided in `col` of size `cnt`
*/
void SparMatSym2::vecMulAddColIso3D(const real* X, real* Y, size_t jj, Element col[], size_t cnt) const
{
    assert_true( cnt > 0 );
    const real X0 = X[jj  ];
    const real X1 = X[jj+1];
    const real X2 = X[jj+2];
    real Y0 = Y[jj  ] + col[0].val * X0;
    real Y1 = Y[jj+1] + col[0].val * X1;
    real Y2 = Y[jj+2] + col[0].val * X2;
    for ( size_t n = 1 ; n < cnt ; ++n )
    {
        const size_t ii = 3 * col[n].inx;
        const real a = col[n].val;
        Y[ii  ] += a * X0;
        Y[ii+1] += a * X1;
        Y[ii+2] += a * X2;
        assert_true( ii > jj );
        Y0 += a * X[ii  ];
        Y1 += a * X[ii+1];
        Y2 += a * X[ii+2];
    }
    Y[jj  ] = Y0;
    Y[jj+1] = Y1;
    Y[jj+1] = Y2;
}


//------------------------------------------------------------------------------
#pragma mark - Prepare Multiplication


#if !MATRIX2_OPTIMIZED_MULTIPLY

bool SparMatSym2::prepareForMultiply(int)
{
    return true;
}

#else


#if MATRIX2_USES_COLNEXT
void SparMatSym2::setNextColumn()
{
    next_[size_] = size_;

    if ( size_ > 0 )
    {
        size_t inx = size_;
        size_t nxt = size_;
        while ( --inx > 0 )
        {
            if ( colsiz_[inx] > 0 )
                nxt = inx;
            next_[inx] = nxt;
        }
        if ( colsiz_[0] > 0 )
            next_[0] = 0;
        else
            next_[0] = nxt;
    }
}
#endif


bool SparMatSym2::prepareForMultiply(int dim)
{
    assert_true( size_ <= alloc_ );
    
#if MATRIX2_USES_COLNEXT
    setNextColumn();
#endif
    
#if ( 0 )
    size_t cnt = 0;
    for ( size_t jj = 0; jj < size_; ++jj )
        cnt += ( colsiz_[jj] == 0 );
    std::clog << "SparMatSym2 has " << cnt << " / " << size_ << " empty columns\n";
#endif

    //count number of non-zero elements, always including the diagonal term
    size_t nbe = 0;
    for ( size_t jj = 0; jj < size_; ++jj )
    {
        if ( colsiz_[jj] > 0 )
            nbe += colsiz_[jj];
        else
            nbe ++;
    }
    
    // allocate DSS sparse matrix storage
    if ( nbe > alcDSS_ )
    {
        constexpr size_t chunk = 16;
        alcDSS_ = ( nbe + chunk - 1 ) & ~( chunk -1 );
        delete[] colDSS_;
        free_real(valDSS_);
        colDSS_ = new unsigned[alcDSS_];
        valDSS_ = new_real(alcDSS_);
    }
    
    /*
     Create the DSS sparse matrix storage.
     */
    unsigned inx = 0;
    for ( size_t jj = 0; jj < size_; ++jj )
    {
        rowDSS_[jj] = inx;
        if ( colsiz_[jj] > 0 )
        {
            Element * col = column_[jj];
            assert_true( col[0].inx == jj );
            for ( size_t n = 0; n < colsiz_[jj]; ++n )
            {
                assert_true( inx < alcDSS_ );
                valDSS_[inx] = col[n].val;
                colDSS_[inx] = col[n].inx * dim;
                ++inx;
            }
        }
        else {
            valDSS_[inx] = 0.0;
            colDSS_[inx] = jj * dim;
            ++inx;
        }
    }
    if ( inx != nbe ) ABORT_NOW("internal error");
    rowDSS_[size_] = inx;
    
    //printSparse(std::clog, 0, 0, 4);
    //printSparseArray(std::clog);
    return true;
}

//------------------------------------------------------------------------------
#pragma mark - Optimized Column-Vector multiplication


void SparMatSym2::vecMulAddCol(const real* X, real* Y,
                               size_t start, size_t stop) const
{
    assert_true( start <= stop );
    size_t jj = colDSS_[start];
    real X0 = X[jj];
    real Y0 = Y[jj] + valDSS_[start] * X0;
    for ( size_t n = start+1; n < stop; ++n )
    {
        real a = valDSS_[n];
        size_t ii = colDSS_[n];
        Y[ii] += a * X0;
        Y0    += a * X[ii];
    }
    Y[jj] = Y0;
}

void SparMatSym2::vecMulAddColIso2D(const real* X, real* Y,
                                    size_t start, size_t stop) const
{
    assert_true( start <= stop );
    size_t jj = colDSS_[start];
    assert_true( stop <= 2*size_ );
    real X0 = X[jj  ];
    real X1 = X[jj+1];
    real Y0 = Y[jj  ] + valDSS_[start] * X0;
    real Y1 = Y[jj+1] + valDSS_[start] * X1;
    for ( size_t n = start+1; n < stop; ++n )
    {
        size_t ii = valDSS_[n];
        assert_true( ii > jj );
        real a = colDSS_[n];
        Y0      += a * X[ii  ];
        Y1      += a * X[ii+1];
        Y[ii  ] += a * X0;
        Y[ii+1] += a * X1;
    }
    Y[jj  ] = Y0;
    Y[jj+1] = Y1;
}


void SparMatSym2::vecMulAddColIso3D(const real* X, real* Y,
                                    size_t start, size_t stop) const
{
    assert_true( start <= stop );
    size_t jj = colDSS_[start];
    real X0 = X[jj  ];
    real X1 = X[jj+1];
    real X2 = X[jj+2];
    real Y0 = Y[jj  ] + valDSS_[start] * X0;
    real Y1 = Y[jj+1] + valDSS_[start] * X1;
    real Y2 = Y[jj+2] + valDSS_[start] * X2;
    for ( size_t n = start+1; n < stop; ++n )
    {
        size_t ii = colDSS_[n];
        assert_true( ii > jj );
        real a = valDSS_[n];
        Y0      += a * X[ii  ];
        Y1      += a * X[ii+1];
        Y2      += a * X[ii+2];
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

#if MATRIX2_USES_SSE

inline void multiply2(const real* X, real* Y, size_t ii,
                      const real* val, vec2 const& xx, vec2& ss)
{
    vec2 aa = loaddup2(val);
    ss = fmadd2(load2(X+ii), aa, ss);
    store2(Y+ii, fmadd2(xx, aa, load2(Y+ii)));
}


void SparMatSym2::vecMulAddColIso2D_SSE(const real* X, real* Y,
                                        size_t start, size_t stop) const
{
    assert_true( start <= stop );
    size_t jj = colDSS_[start];
    const vec2 xx = load2(X+jj);
    vec2 ss = fmadd2(loaddup2(valDSS_+start), xx, load2(Y+jj));
    // there is a dependence here for 'ss'
    for ( size_t n = start+1; n < stop; ++n )
        multiply2(X, Y, colDSS_[n], valDSS_+n, xx, ss);
    store2(Y+jj, ss);
}


void SparMatSym2::vecMulAddColIso2D_SSEU(const real* X, real* Y,
                                         size_t start, size_t stop) const
{
    assert_true( start <= stop );
    size_t jj = colDSS_[start];
    const vec2 xx = load2(X+jj);
    vec2 s0 = mul2(loaddup2(valDSS_+start), xx);
    vec2 s1 = load2(Y+jj);
    vec2 s2 = setzero2();
    vec2 s3 = setzero2();
    
    size_t n = start+1;
#if ( 0 )
    // unrolling by 8 may exceed the number of registers in the CPU
#pragma nounroll
    if ( end >= n + 8 )
    {
        vec2 s4 = setzero2();
        vec2 s5 = setzero2();
        vec2 s6 = setzero2();
        vec2 s7 = setzero2();
        size_t end = n + 8 * ( ( stop - n ) / 8 );
        // process 8 by 8:
        for ( ; n < end; n += 8 )
        {
            const size_t i0 = colDSS_[n  ];
            const size_t i1 = colDSS_[n+1];
            const size_t i2 = colDSS_[n+2];
            const size_t i3 = colDSS_[n+3];
            const size_t i4 = colDSS_[n+4];
            const size_t i5 = colDSS_[n+5];
            const size_t i6 = colDSS_[n+6];
            const size_t i7 = colDSS_[n+7];
            vec2 y0 = load2(Y+i0);
            vec2 y1 = load2(Y+i1);
            vec2 y2 = load2(Y+i2);
            vec2 y3 = load2(Y+i3);
            vec2 y4 = load2(Y+i4);
            vec2 y5 = load2(Y+i5);
            vec2 y6 = load2(Y+i6);
            vec2 y7 = load2(Y+i7);
            vec2 a0 = loaddup2(valDSS_+n);
            vec2 a1 = loaddup2(valDSS_+n+1);
            vec2 a2 = loaddup2(valDSS_+n+2);
            vec2 a3 = loaddup2(valDSS_+n+3);
            vec2 a4 = loaddup2(valDSS_+n+4);
            vec2 a5 = loaddup2(valDSS_+n+5);
            vec2 a6 = loaddup2(valDSS_+n+6);
            vec2 a7 = loaddup2(valDSS_+n+7);
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
    
    size_t end = n + 4 * ( ( stop - n ) / 4 );
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
        multiply2(X, Y, colDSS_[n  ], valDSS_+n  , xx, s0);
        multiply2(X, Y, colDSS_[n+1], valDSS_+n+1, xx, s1);
        multiply2(X, Y, colDSS_[n+2], valDSS_+n+2, xx, s2);
        multiply2(X, Y, colDSS_[n+3], valDSS_+n+3, xx, s3);
#else
        /* we remove here the apparent dependency on the values of Y[],
         which are read and written, but at different indices.
         The compiler can reorder instructions to avoid lattencies */
        const size_t i0 = colDSS_[n  ];
        const size_t i1 = colDSS_[n+1];
        const size_t i2 = colDSS_[n+2];
        const size_t i3 = colDSS_[n+3];
        vec2 y0 = load2(Y+i0);
        vec2 y1 = load2(Y+i1);
        vec2 y2 = load2(Y+i2);
        vec2 y3 = load2(Y+i3);
        vec2 a0 = loaddup2(valDSS_+n);
        vec2 a1 = loaddup2(valDSS_+n+1);
        vec2 a2 = loaddup2(valDSS_+n+2);
        vec2 a3 = loaddup2(valDSS_+n+3);
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
        multiply2(X, Y, colDSS_[n], valDSS_+n, xx, s0);
    store2(Y+jj, s0);
}

#endif

#if MATRIX2_USES_AVX

/*
Accumulation is done here in the higher part of 'ss'
The high position of 'xx' is not used
The low position of 'ss' is used locally
*/
inline void multiply4(const real* X, real* Y, size_t ii,
                      const real* val, vec4 const& xx, vec4& ss)
{
    vec4 x = blend4(xx, broadcast2(X+ii), 0b1100);  // hi <- X , lo <- xx
    ss = blend4(cast4(load2(Y+ii)), ss, 0b1100);    // hi <- ss, lo <- Y
    ss = fmadd4(broadcast1(val), x, ss);
    store2(Y+ii, getlo(ss));
}


void SparMatSym2::vecMulAddColIso2D_AVX(const real* X, real* Y,
                                        size_t start, size_t stop) const
{
    assert_true( start <= stop );
    size_t jj = colDSS_[start];
    const vec4 xx = broadcast2(X+jj);  // hi position
    vec4 ss = fmadd4(broadcast1(valDSS_+start), xx, broadcast2(Y+jj));
    // there is a dependence here for 'ss'
    for ( size_t n = start+1; n < stop; ++n )
        multiply4(X, Y, colDSS_[n], valDSS_+n, xx, ss);
    store2(Y+jj, gethi(ss));
}


void SparMatSym2::vecMulAddColIso2D_AVXU(const real* X, real* Y,
                                         size_t start, size_t stop) const
{
    assert_true( start <= stop );
    size_t jj = colDSS_[start];
    const vec4 xx = broadcast2(X+jj);  // hi and lo position
    vec4 s0 = mul4(broadcast1(valDSS_+start), xx);
    vec4 s1 = broadcast2(Y+jj);
    vec4 s2 = setzero4();
    vec4 s3 = setzero4();
    
    unsigned * inx = colDSS_ + start;
    const real * val = valDSS_ + start + 1;
    const real * end = val + 4 * ((stop-start)/4);
    // process 4 by 4:
#pragma nounroll
    for ( ; val < end; val += 4 )
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
        const size_t i0 = inx[0];
        const size_t i1 = inx[1];
        const size_t i2 = inx[2];
        const size_t i3 = inx[3];
        s0 = blend4(cast4(load2(Y+i0)),s0,0b1100);    // lo = Y
        s1 = blend4(cast4(load2(Y+i1)),s1,0b1100);    // lo = Y
        s2 = blend4(cast4(load2(Y+i2)),s2,0b1100);    // lo = Y
        s3 = blend4(cast4(load2(Y+i3)),s3,0b1100);    // lo = Y
        vec4 x0 = blend4(xx,broadcast2(X+i0),0b1100);   // hi = X , lo <- xx
        vec4 x1 = blend4(xx,broadcast2(X+i1),0b1100);   // hi = X , lo <- xx
        vec4 x2 = blend4(xx,broadcast2(X+i2),0b1100);   // hi = X , lo <- xx
        vec4 x3 = blend4(xx,broadcast2(X+i3),0b1100);   // hi = X , lo <- xx
        s0 = fmadd4(broadcast1(val  ), x0, s0);
        s1 = fmadd4(broadcast1(val+1), x1, s1);
        s2 = fmadd4(broadcast1(val+2), x2, s2);
        s3 = fmadd4(broadcast1(val+3), x3, s3);
        store2(Y+i0, getlo(s0));
        store2(Y+i1, getlo(s1));
        store2(Y+i2, getlo(s2));
        store2(Y+i3, getlo(s3));
#endif
        inx += 4;
    }
    // collapse into 's0'
    s0 = add4(add4(s0,s1), add4(s2,s3));
    // process remaining values:
    end = valDSS_ + stop;
#pragma nounroll
    for ( ; val < end; ++val, ++inx )
        multiply4(X, Y, inx[0], val, xx, s0);
    store2(Y+jj, gethi(s0));
}

#endif
#endif

//------------------------------------------------------------------------------
#pragma mark - 3D SIMD

#if MATRIX2_USES_AVX && MATRIX2_OPTIMIZED_MULTIPLY
void SparMatSym2::vecMulAddColIso3D_AVX(const real* X, real* Y,
                                        size_t start, size_t stop) const
{
    assert_true( start <= stop );
    size_t jj = colDSS_[start];
    unsigned * inx = colDSS_ + start;
    real const* val = valDSS_ + start;
    real const* end = valDSS_ + stop;
    
    //printf("SparMatSym2 column %lu has %lu elements\n", jj, stop - start);
    vec4 zz = setzero4();
    vec4 xx = blend4(loadu4(X+jj), zz, 0b1000);
    vec4 yy = fmadd4(broadcast1(val), xx, loadu4(Y+jj));
    // process one element when the number of values is even
    if ( 0 == (( stop - start ) & 1) )
    {
        size_t ii = *(++inx);
        vec4 aa = broadcast1(++val);
        zz = mul4(aa, loadu4(X+ii));
        storeu4(Y+ii, fmadd4(aa, xx, loadu4(Y+ii)));
    }
    ++val;
    ++inx;
    while ( val < end )
    {
        size_t ii = *(inx  );
        size_t kk = *(inx+1);
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


//------------------------------------------------------------------------------
#pragma mark - Matrix-Vector Add-multiply


void SparMatSym2::vecMulAdd(const real* X, real* Y, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= size_ );

#if MATRIX2_USES_COLNEXT
    for ( size_t jj = next_[start]; jj < stop; jj = next_[jj+1] )
#else
    for ( size_t jj = start; jj < stop; ++jj )
#endif
    {
#if MATRIX2_OPTIMIZED_MULTIPLY
        vecMulAddCol(X, Y, rowDSS_[jj], rowDSS_[jj+1]);
#else
        if ( colsiz_[jj] > 0 )
        {
            assert_true(column_[jj][0].inx == jj);
            vecMulAddCol(X, Y, jj, column_[jj], colsiz_[jj]);
        }
#endif
    }
}


void SparMatSym2::vecMulAddIso2D(const real* X, real* Y, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= size_ );

#if MATRIX2_USES_COLNEXT
    for ( size_t jj = next_[start]; jj < stop; jj = next_[jj+1] )
#else
    for ( size_t jj = start; jj < stop; ++jj )
#endif
    {
#if MATRIX2_OPTIMIZED_MULTIPLY
#  if MATRIX2_USES_AVX
        vecMulAddColIso2D_AVXU(X, Y, rowDSS_[jj], rowDSS_[jj+1]);
#  elif MATRIX2_USES_SSE
        vecMulAddColIso2D_SSEU(X, Y, rowDSS_[jj], rowDSS_[jj+1]);
#  else
        vecMulAddColIso2D(X, Y, rowDSS_[jj], rowDSS_[jj+1]);
#  endif
#else
        if ( colsiz_[jj] > 0 )
        {
            assert_true(column_[jj][0].inx == jj);
            vecMulAddColIso2D(X, Y, 2*jj, column_[jj], colsiz_[jj]);
        }
#endif

    }
}


void SparMatSym2::vecMulAddIso3D(const real* X, real* Y, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= size_ );

#if MATRIX2_USES_COLNEXT
    for ( size_t jj = next_[start]; jj < stop; jj = next_[jj+1] )
#else
    for ( size_t jj = start; jj < stop; ++jj )
#endif
    {
#if MATRIX2_OPTIMIZED_MULTIPLY
#  if MATRIX2_USES_AVX
        vecMulAddColIso3D_AVX(X, Y, rowDSS_[jj], rowDSS_[jj+1]);
#  else
        vecMulAddColIso3D(X, Y, rowDSS_[jj], rowDSS_[jj+1]);
#  endif
#else
        if ( colsiz_[jj] > 0 )
        {
            assert_true(column_[jj][0].inx == jj);
            vecMulAddColIso3D(X, Y, 3*jj, column_[jj], colsiz_[jj]);
        }
#endif
    }
}


//------------------------------------------------------------------------------
#pragma mark - Matrix-Vector multiplication

void SparMatSym2::vecMul(const real* X, real* Y, size_t start, size_t stop) const
{
    zero_real(stop-start, Y+start);
    vecMulAdd(X, Y, start, stop);
}
