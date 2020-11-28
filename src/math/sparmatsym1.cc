// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#include <cmath>
#include "sparmatsym1.h"
#include "assert_macro.h"
#include "blas.h"

#include <iomanip>
#include <sstream>


#ifdef __AVX__
#  define MATRIX1_USES_AVX REAL_IS_DOUBLE
#  define MATRIX1_USES_SSE REAL_IS_DOUBLE
#  include "simd.h"
#elif defined(__SSE3__)
#  define MATRIX1_USES_AVX 0
#  define MATRIX1_USES_SSE REAL_IS_DOUBLE
#  include "simd.h"
#else
#  define MATRIX1_USES_AVX 0
#  define MATRIX1_USES_SSE 0
#endif


SparMatSym1::SparMatSym1()
{
    size_   = 0;
    alloc_  = 0;
    column_ = nullptr;
    colsiz_ = nullptr;
    colmax_ = nullptr;
#if MATRIX1_OPTIMIZED_MULTIPLY
    nmax_ = 0;
    ija_  = nullptr;
    sa_   = nullptr;
#endif
#if MATRIX1_USES_COLNEXT
    colidx_ = new unsigned[2];
    colidx_[0] = 0;
#endif
}


void SparMatSym1::allocate(size_t alc)
{
    if ( alc > alloc_ )
    {
        /*
         'chunk' can be increased to gain performance:
         more memory will be used, but reallocation will be less frequent
         */
        constexpr size_t chunk = 64;
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );

        //fprintf(stderr, "SMS1 allocate matrix %u\n", alc);
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
        
#if MATRIX1_USES_COLNEXT
        delete[] colidx_;
        colidx_ = new unsigned[alc+1];
        for ( size_t n = 0; n <= alc; ++n )
            colidx_[n] = n;
#endif
    }
}


void SparMatSym1::deallocate()
{
    if ( column_ )
    {
        for ( size_t ii = 0; ii < alloc_; ++ii )
            delete[] column_[ii];
        delete[] column_; column_ = nullptr;
        delete[] colsiz_; colsiz_ = nullptr;
        delete[] colmax_; colmax_ = nullptr;
#if MATRIX1_USES_COLNEXT
        delete[] colidx_; colidx_ = nullptr;
#endif
#if MATRIX1_OPTIMIZED_MULTIPLY
        delete[] ija_;
        free_real(sa_);
        ija_ = nullptr;
        sa_ = nullptr;
#endif
    }
    alloc_ = 0;
}


/// copy `cnt` elements from `src` to `dst`
void copy(size_t cnt, SparMatSym1::Element * src, SparMatSym1::Element * dst)
{
    for ( size_t ii = 0; ii < cnt; ++ii )
        dst[ii] = src[ii];
}

/// move `cnt` elements to the next index, starting at vec[0]
void shift(size_t cnt, SparMatSym1::Element * vec)
{
    for ( size_t ii = cnt; ii > 0; --ii )
        vec[ii] = vec[ii-1];
}


void SparMatSym1::allocateColumn(const size_t jj, size_t alc)
{
    assert_true( jj < size_ );
    if ( alc > colmax_[jj] )
    {
        //fprintf(stderr, "SMS1 allocate column %i size %u\n", jj, alc);
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


SparMatSym1::Element * SparMatSym1::insertElement(const size_t jj, size_t inx)
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


real& SparMatSym1::diagonal(size_t i)
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
real& SparMatSym1::operator()(size_t i, size_t j)
{
    assert_true( i < size_ );
    assert_true( j < size_ );
    //fprintf(stderr, "SMS1( %6i %6i )\n", i, j);
    
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


real* SparMatSym1::addr(size_t i, size_t j) const
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

void SparMatSym1::reset()
{
    for ( size_t jj = 0; jj < size_; ++jj )
        colsiz_[jj] = 0;
}


bool SparMatSym1::isNotZero() const
{
    //check for any non-zero sparse term:
    for ( size_t jj = 0; jj < size_; ++jj )
        for ( size_t kk = 0 ; kk < colsiz_[jj] ; ++kk )
            if ( column_[jj][kk].val != 0 )
                return true;
    
    //if here, the matrix is empty
    return false;
}


void SparMatSym1::scale(const real alpha)
{
    for ( size_t jj = 0; jj < size_; ++jj )
        for ( size_t n = 0; n < colsiz_[jj]; ++n )
            column_[jj][n].val *= alpha;
}


void SparMatSym1::addDiagonalBlock(real* mat, const size_t ldd,
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
                //printf("SMS1 %4i %4i % .4f\n", ii, jj, a);
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
void SparMatSym1::addTriangularBlockBanded(real alpha, real* mat, const size_t ldd,
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
            assert_true( ii >= jj );
            // assuming lower triangle is stored:
            if ( ii < end )
            {
                size_t i = ii - start;
                //printf("SMS1 %4i %4i % .4f\n", ii, jj, a);
                assert_true( i >= j );
                // with banded storage, mat(i, j) is stored in mat[i-j+ldd*j]
                if ( i <= j + rank )
                    mat[i-j+ldd*j] += alpha * column_[jj][n].val;
            }
        }
    }
}



void SparMatSym1::addDiagonalTrace(real alpha, real* mat, const size_t ldd,
                                   const size_t start, const size_t cnt) const
{
    fprintf(stderr, "unfinished SparMatSym1::addDiagonalTrace()\n");
    exit(1);
}


void SparMatSym1::addDiagonalTraceBanded(real alpha, real* mat, const size_t ldd,
                                         const size_t start, const size_t cnt,
                                         const size_t rank) const
{
    fprintf(stderr, "unfinished SparMatSym1::addDiagonalTraceBanded()\n");
    exit(1);
}


int SparMatSym1::bad() const
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


size_t SparMatSym1::nbElements(size_t start, size_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= size_ );
    //all allocated elements are counted, even if zero
    size_t cnt = 0;
    for ( size_t jj = start; jj < stop; ++jj )
        cnt += colsiz_[jj];
    return cnt;
}


size_t SparMatSym1::nbDiagonalElements(size_t start, size_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= size_ );
    //all allocated elements are counted, even if zero
    size_t cnt = 0;
    for ( size_t jj = start; jj < stop; ++jj )
        cnt += ( colsiz_[jj] > 0 ) && ( column_[jj][0].val != 0.0 );
    return cnt;
}


std::string SparMatSym1::what() const
{
    std::ostringstream msg;
#if MATRIX1_USES_AVX
    msg << "SMS1x " << nbElements();
#elif MATRIX1_USES_SSE
    msg << "SMS1e " << nbElements();
#else
    msg << "SMS1 " << nbElements();
#endif
    msg << " (" << nbDiagonalElements(0, size_) << ")";
    return msg.str();
}


void SparMatSym1::printSparse(std::ostream& os, real inf, size_t start, size_t stop) const
{
    os << "% SparMatSym1 size " << size_ << ":";
    stop = std::min(stop, size_);
    char str[256];
    for ( size_t jj = start; jj < stop; ++jj )
    {
        if ( colsiz_[jj] > 0 )
        {
            os << "% column " << jj << "\n";
            for ( size_t n = 0 ; n < colsiz_[jj] ; ++n )
            {
                real v = column_[jj][n].val;
                if ( abs_real(v) >= inf )
                {
                    snprintf(str, sizeof(str), "%6u %6lu %16.6f\n", column_[jj][n].inx, jj, v);
                    os << str;
                }
            }
        }
    }
}


void SparMatSym1::printColumns(std::ostream& os, size_t start, size_t stop)
{
    stop = std::min(stop, size_);
    os << "% SparMatSym1 size " << size_ << ":";
    for ( size_t jj = start; jj < stop; ++jj )
    {
        os << "\n   " << jj << "   " << colsiz_[jj];
#if MATRIX1_USES_COLNEXT
        os << " index " << colidx_[jj];
#endif
    }
    std::endl(os);
}


void SparMatSym1::printColumn(std::ostream& os, const size_t jj)
{
    Element const* col = column_[jj];
    os << "SMS1 col " << jj << ":";
    for ( size_t n = 0; n < colsiz_[jj]; ++n )
    {
        os << "\n" << col[n].inx << " :";
        os << " " << col[n].val;
    }
    std::endl(os);
}


void SparMatSym1::printSparseArray(std::ostream& os) const
{
#if MATRIX1_OPTIMIZED_MULTIPLY
    std::ios::fmtflags fgs = os.flags();
    size_t end = ija_[size_];
    
    os << "ija ";
    for ( size_t n = 0; n < end; ++n )
        os << " " << std::setw(6) << ija_[n];
    os << "\n";
    
    std::streamsize p = os.precision();
    os.precision(2);
    os << "sa  ";
    for ( size_t n = 0; n < end; ++n )
        os << " " << std::setw(6) << sa_[n];
    os << "\n";
    os.precision(p);
    os.setf(fgs);
#else
    os << "optimized sparse matrix storage unavailable\n";
#endif
}


//------------------------------------------------------------------------------
#pragma mark - Column-Vector Multiplication

/**
Multiply by column `jj` provided in `col` of size `cnt`
*/
void SparMatSym1::vecMulAddCol(const real* X, real* Y, size_t jj, Element col[], size_t cnt) const
{
    assert_true( cnt > 0 );
    assert_true( col[0].inx == jj );
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
void SparMatSym1::vecMulAddColIso2D(const real* X, real* Y, size_t jj, Element col[], size_t cnt) const
{
    assert_true( cnt > 0 );
    assert_true( 2*col[0].inx == jj );
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
void SparMatSym1::vecMulAddColIso3D(const real* X, real* Y, size_t jj, Element col[], size_t cnt) const
{
    assert_true( cnt > 0 );
    assert_true( 3*col[0].inx == jj );
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
    Y[jj+2] = Y2;
}


//------------------------------------------------------------------------------
#pragma mark - Prepare Multiplication

#if MATRIX1_USES_COLNEXT
void SparMatSym1::setColumnIndex()
{
    if ( size_ > 0 )
    {
        size_t inx = size_;
        size_t nxt = size_;
        while ( inx-- > 0 )
        {
            if ( colsiz_[inx] > 0 )
                nxt = inx;
            colidx_[inx] = nxt;
        }
    }
    colidx_[size_] = size_;
}
#else
void SparMatSym1::setColumnIndex()
{
}
#endif


#if !MATRIX1_OPTIMIZED_MULTIPLY

bool SparMatSym1::prepareForMultiply(int)
{
    setColumnIndex();
    return true;
}

#else

bool SparMatSym1::prepareForMultiply(int dim)
{
    assert_true( size_ <= alloc_ );
    
    setColumnIndex();
    
#if ( 0 )
    size_t cnt = 0;
    for ( size_t jj = 0; jj < size_; ++jj )
        cnt += ( colsiz_[jj] == 0 );
    std::clog << "SparMatSym1 has " << cnt << " / " << size_ << " empty columns\n";
#endif

    //count number of non-zero elements, including diagonal
    size_t nbe = 1;
    for ( size_t jj = 0; jj < size_; ++jj )
    {
        if ( colsiz_[jj] > 0 )
            nbe += colsiz_[jj];
        else
            nbe ++;
    }
    
    //allocate classical sparse matrix storage (Numerical Recipes)
    if ( nbe > nmax_ )
    {
        delete[] ija_;
        free_real(sa_);

        nmax_ = nbe + size_;
        ija_  = new unsigned[nmax_];
        sa_   = new_real(nmax_);
    }
    
    /*
     Create the compressed sparse format described in Numerical Recipe,
     Chapter 2.7 Sparse Linear Systems - Indexed Storage of Sparse Matrices
     indices however start here at zero, and everything is shifted by one index,
     compared to numerical recipe's code.
     */
    ija_[0] = size_+1;
    sa_[size_] = 42; // this is the arbitrary value
    size_t inx = size_;
    for ( size_t jj = 0; jj < size_; ++jj )
    {
        if ( colsiz_[jj] > 0 )
        {
            // diagonal term first:
            assert_true( column_[jj][0].inx == jj );
            sa_[jj] = column_[jj][0].val;
            // other non-zero elements:
            for ( size_t cc = 1; cc < colsiz_[jj]; ++cc )
            {
                ++inx;
                assert_true( inx < nbe );
                sa_[inx]  = column_[jj][cc].val;
                ija_[inx] = dim * column_[jj][cc].inx;
            }
        }
        else {
            sa_[jj] = 0.0;
        }
        ija_[jj+1] = inx+1;
    }
    if ( inx+1 != nbe ) ABORT_NOW("internal error");

    //printSparse(std::clog);
    //printSparseArray(std::clog);
    return true;
}

//------------------------------------------------------------------------------
#pragma mark - Optimized Column-Vector multiplication


void SparMatSym1::vecMulAddCol(const real* X, real* Y, size_t jj,
                               real const* dia, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= nmax_ );
    real X0 = X[jj];
    real Y0 = Y[jj] + dia[0] * X0;
    for ( size_t n = start; n < stop; ++n )
    {
        real a = sa_[n];
        size_t ii = ija_[n];
        Y[ii] += a * X0;
        Y0    += a * X[ii];
    }
    Y[jj] = Y0;
}

void SparMatSym1::vecMulAddColIso2D(const real* X, real* Y, size_t jj,
                                    real const* dia, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= nmax_ );
    real X0 = X[jj  ];
    real X1 = X[jj+1];
    real Y0 = Y[jj  ] + dia[0] * X0;
    real Y1 = Y[jj+1] + dia[0] * X1;
    for ( size_t n = start; n < stop; ++n )
    {
        size_t ii = ija_[n];
        assert_true( ii > jj );
        real a = sa_[n];
        Y0      += a * X[ii  ];
        Y1      += a * X[ii+1];
        Y[ii  ] += a * X0;
        Y[ii+1] += a * X1;
    }
    Y[jj  ] = Y0;
    Y[jj+1] = Y1;
}


void SparMatSym1::vecMulAddColIso3D(const real* X, real* Y, size_t jj,
                                    real const* dia, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= nmax_ );
    real X0 = X[jj  ];
    real X1 = X[jj+1];
    real X2 = X[jj+2];
    real Y0 = Y[jj  ] + dia[0] * X0;
    real Y1 = Y[jj+1] + dia[0] * X1;
    real Y2 = Y[jj+2] + dia[0] * X2;
    for ( size_t n = start; n < stop; ++n )
    {
        size_t ii = ija_[n];
        assert_true( ii > jj );
        real a = sa_[n];
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

#if MATRIX1_USES_SSE

inline void multiply2(const real* X, real* Y, size_t ii,
                      const real* val, vec2 const& xx, vec2& ss)
{
    vec2 aa = loaddup2(val);
    ss = fmadd2(load2(X+ii), aa, ss);
    store2(Y+ii, fmadd2(xx, aa, load2(Y+ii)));
}


void SparMatSym1::vecMulAddColIso2D_SSE(const real* X, real* Y, size_t jj,
                                        real const* dia, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= nmax_ );
    const vec2 xx = load2(X+jj);
    vec2 ss = fmadd2(loaddup2(dia), xx, load2(Y+jj));
    // there is a dependence here for 'ss'
    for ( size_t n = start; n < stop; ++n )
        multiply2(X, Y, ija_[n], sa_+n, xx, ss);
    store2(Y+jj, ss);
}


void SparMatSym1::vecMulAddColIso2D_SSEU(const real* X, real* Y, size_t jj,
                                         real const* dia, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= nmax_ );
    const vec2 xx = load2(X+jj);
    vec2 s0 = mul2(loaddup2(dia), xx);
    vec2 s1 = load2(Y+jj);
    vec2 s2 = setzero2();
    vec2 s3 = setzero2();
    
    size_t n = start;
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
            const size_t i0 = ija_[n  ];
            const size_t i1 = ija_[n+1];
            const size_t i2 = ija_[n+2];
            const size_t i3 = ija_[n+3];
            const size_t i4 = ija_[n+4];
            const size_t i5 = ija_[n+5];
            const size_t i6 = ija_[n+6];
            const size_t i7 = ija_[n+7];
            vec2 y0 = load2(Y+i0);
            vec2 y1 = load2(Y+i1);
            vec2 y2 = load2(Y+i2);
            vec2 y3 = load2(Y+i3);
            vec2 y4 = load2(Y+i4);
            vec2 y5 = load2(Y+i5);
            vec2 y6 = load2(Y+i6);
            vec2 y7 = load2(Y+i7);
            vec2 a0 = loaddup2(sa_+n);
            vec2 a1 = loaddup2(sa_+n+1);
            vec2 a2 = loaddup2(sa_+n+2);
            vec2 a3 = loaddup2(sa_+n+3);
            vec2 a4 = loaddup2(sa_+n+4);
            vec2 a5 = loaddup2(sa_+n+5);
            vec2 a6 = loaddup2(sa_+n+6);
            vec2 a7 = loaddup2(sa_+n+7);
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
        multiply2(X, Y, ija_[n  ], sa_+n  , xx, s0);
        multiply2(X, Y, ija_[n+1], sa_+n+1, xx, s1);
        multiply2(X, Y, ija_[n+2], sa_+n+2, xx, s2);
        multiply2(X, Y, ija_[n+3], sa_+n+3, xx, s3);
#else
        /* we remove here the apparent dependency on the values of Y[],
         which are read and written, but at different indices.
         The compiler can reorder instructions to avoid lattencies */
        const size_t i0 = ija_[n  ];
        const size_t i1 = ija_[n+1];
        const size_t i2 = ija_[n+2];
        const size_t i3 = ija_[n+3];
        vec2 y0 = load2(Y+i0);
        vec2 y1 = load2(Y+i1);
        vec2 y2 = load2(Y+i2);
        vec2 y3 = load2(Y+i3);
        vec2 a0 = loaddup2(sa_+n);
        vec2 a1 = loaddup2(sa_+n+1);
        vec2 a2 = loaddup2(sa_+n+2);
        vec2 a3 = loaddup2(sa_+n+3);
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
        multiply2(X, Y, ija_[n], sa_+n, xx, s0);
    store2(Y+jj, s0);
}

#endif

#if MATRIX1_USES_AVX && MATRIX1_OPTIMIZED_MULTIPLY

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


void SparMatSym1::vecMulAddColIso2D_AVX(const real* X, real* Y, size_t jj,
                                        real const* dia, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= nmax_ );
    const vec4 xx = broadcast2(X+jj);  // hi position
    vec4 ss = fmadd4(broadcast1(dia), xx, broadcast2(Y+jj));
    // there is a dependence here for 'ss'
    for ( size_t n = start; n < stop; ++n )
        multiply4(X, Y, ija_[n], sa_+n, xx, ss);
    store2(Y+jj, gethi(ss));
}


void SparMatSym1::vecMulAddColIso2D_AVXU(const real* X, real* Y, size_t jj,
                                         real const* dia, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= nmax_ );
    const vec4 xx = broadcast2(X+jj);  // hi and lo position
    vec4 s0 = mul4(broadcast1(dia), xx);
    vec4 s1 = broadcast2(Y+jj);
    vec4 s2 = setzero4();
    vec4 s3 = setzero4();
    
    unsigned * inx = ija_ + start;
    const real * val = sa_ + start;
    const real * end = sa_ + stop;
    const real * halt = end - 3;  // val+3 <= end-1  is  val < end-3;
        // process 4 by 4:
    #pragma nounroll
    for ( ; val < halt; val += 4 )
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
#pragma nounroll
    for ( ; val < end; ++val, ++inx )
        multiply4(X, Y, inx[0], val, xx, s0);
    store2(Y+jj, gethi(s0));
}

#endif
#endif

//------------------------------------------------------------------------------
#pragma mark - 3D SIMD

#if MATRIX1_USES_SSE
#endif


#if MATRIX1_USES_AVX && MATRIX1_OPTIMIZED_MULTIPLY
void SparMatSym1::vecMulAddColIso3D_AVX(const real* X, real* Y, size_t jj,
                                        real const* dia, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= nmax_ );
    vec4 zz = setzero4();
    vec4 xx = blend4(loadu4(X+jj), zz, 0b1000);
    vec4 yy = fmadd4(broadcast1(dia), xx, loadu4(Y+jj));
    real * val = sa_ + start;
    real const*end = sa_ + stop - 1;
    unsigned *inx = ija_ + start;
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
    yy = add4(yy, zz);
    while ( val <= end )
    {
        size_t ii = *inx++;
        assert_true( ii > jj );
        vec4 aa = broadcast1(val++);
        yy = fmadd4(aa, loadu4(X+ii), yy);
        storeu4(Y+ii, fmadd4(aa, xx, loadu4(Y+ii)));
    }
    store3(Y+jj, yy);
}
#endif


//------------------------------------------------------------------------------
#pragma mark - Matrix-Vector Add-multiply


void SparMatSym1::vecMulAdd(const real* X, real* Y, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_);

#if MATRIX1_USES_COLNEXT
    for ( size_t jj = colidx_[start]; jj < stop; jj = colidx_[jj+1] )
#else
    for ( size_t jj = start; jj < stop; ++jj )
#endif
    {
#if MATRIX1_OPTIMIZED_MULTIPLY
        vecMulAddCol(X, Y, jj, sa_+jj, ija_[jj], ija_[jj+1]);
#else
        if ( colsiz_[jj] > 0 )
        {
            assert_true(column_[jj][0].inx == jj);
            vecMulAddCol(X, Y, jj, column_[jj], colsiz_[jj]);
        }
#endif
    }
}


void SparMatSym1::vecMulAddIso2D(const real* X, real* Y, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_);

#if MATRIX1_USES_COLNEXT
    for ( size_t jj = colidx_[start]; jj < stop; jj = colidx_[jj+1] )
#else
    for ( size_t jj = start; jj < stop; ++jj )
#endif
    {
#if MATRIX1_OPTIMIZED_MULTIPLY
#  if MATRIX1_USES_AVX
        vecMulAddColIso2D_AVXU(X, Y, 2*jj, sa_+jj, ija_[jj], ija_[jj+1]);
#  elif MATRIX1_USES_SSE
        vecMulAddColIso2D_SSEU(X, Y, 2*jj, sa_+jj, ija_[jj], ija_[jj+1]);
#  else
        vecMulAddColIso2D(X, Y, 2*jj, sa_+jj, ija_[jj], ija_[jj+1]);
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


void SparMatSym1::vecMulAddIso3D(const real* X, real* Y, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_);

#if MATRIX1_USES_COLNEXT
    for ( size_t jj = colidx_[start]; jj < stop; jj = colidx_[jj+1] )
#else
    for ( size_t jj = start; jj < stop; ++jj )
#endif
    {
#if MATRIX1_OPTIMIZED_MULTIPLY
#  if MATRIX1_USES_AVX
        vecMulAddColIso3D_AVX(X, Y, 3*jj, sa_+jj, ija_[jj], ija_[jj+1]);
#  else
        vecMulAddColIso3D(X, Y, 3*jj, sa_+jj, ija_[jj], ija_[jj+1]);
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

void SparMatSym1::vecMul(const real* X, real* Y, size_t start, size_t stop) const
{
    zero_real(stop-start, Y+start);
    vecMulAdd(X, Y, start, stop);
}
