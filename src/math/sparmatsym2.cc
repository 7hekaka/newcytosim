// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#include <cmath>
#include "sparmatsym2.h"
#include "assert_macro.h"
#include "blas.h"

#include <iomanip>
#include <sstream>


#ifdef __AVX__
#  define MATRIX2_USES_AVX REAL_IS_DOUBLE
#  define MATRIX2_USES_SSE 0
#  include "simd.h"
#elif defined(__SSE3__)
#  define MATRIX2_USES_AVX 0
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
    nmax_ = 0;
    ija_  = nullptr;
    sa_   = nullptr;
#endif
#if MATRIX2_USES_COLNEXT
    next_ = new size_t[1];
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

        //fprintf(stderr, "SMS1 allocate matrix %u\n", alc);
        Element ** col_new    = new Element*[alc];
        size_t   * colsiz_new = new size_t[alc];
        size_t   * colmax_new = new size_t[alc];
        
        size_t ii = 0;
        if ( column_ )
        {
            for ( ; ii < alloc_; ++ii )
            {
                col_new[ii]    = column_[ii];
                colsiz_new[ii] = colsiz_[ii];
                colmax_new[ii] = colmax_[ii];
            }
            delete[] column_;
            delete[] colsiz_;
            delete[] colmax_;
        }
        
        column_ = col_new;
        colsiz_ = colsiz_new;
        colmax_ = colmax_new;
        alloc_ = alc;

        for ( ; ii < alc; ++ii )
        {
            column_[ii] = nullptr;
            colsiz_[ii] = 0;
            colmax_[ii] = 0;
        }
        
#if MATRIX2_USES_COLNEXT
        delete[] next_;
        next_ = new size_t[alc+1];
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
        delete[] ija_;
        free_real(sa_);
        ija_ = nullptr;
        sa_ = nullptr;
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
        //fprintf(stderr, "SMS1 allocate column %i size %u\n", jj, alc);
        constexpr size_t chunk = 16;
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );
        Element * col_new = new Element[alc];
        
        if ( column_[jj] )
        {
            //copy over previous column elements
            copy(colsiz_[jj], column_[jj], col_new);
            
            //release old memory
            delete[] column_[jj];
        }
        column_[jj]  = col_new;
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
        Element * col_new = new Element[alc];
        if ( column_[jj] )
        {
            copy(inx, column_[jj], col_new);
            copy(colsiz_[jj]-inx, column_[jj]+inx, col_new+inx+1);
            delete[] column_[jj];
        }
        column_[jj]  = col_new;
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


real& SparMatSym2::diagonal(size_t ix)
{
    assert_true( ix < size_ );
    
    Element * col;
    
    if ( colsiz_[ix] == 0 )
    {
        allocateColumn(ix, 1);
        col = column_[ix];
        //diagonal term always first:
        col->reset(ix);
        colsiz_[ix] = 1;
    }
    else
    {
        col = column_[ix];
        assert_true( col->inx == ix );
    }
    
    return col->val;
}

/**
 This allocate to be able to hold the matrix element if necessary
*/
real& SparMatSym2::operator()(size_t i, size_t j)
{
    assert_true( i < size_ );
    assert_true( j < size_ );
    //fprintf(stderr, "SMS2( %6i %6i )\n", i, j);
    
    // swap to get ii > jj (address lower triangle)
    size_t ii = std::max(i, j);
    size_t jj = std::min(i, j);

    Element * col = column_[jj];
    
    if ( colsiz_[jj] > 0 )
    {
        Element * e = col;
        Element * lst = col + colsiz_[jj] - 1;
        
        //check all elements in the column:
        while ( e <= lst )
        {
            if ( e->inx == ii )
                return e->val;
            ++e;
        }
    }
    else
    {
        allocateColumn(jj, 2);
        col = column_[jj];
        // put diagonal term always first:
        col->reset(jj);
        if ( ii == jj )
        {
            colsiz_[jj] = 1;
            return col[0].val;
        }
        //add the requested term:
        col[1].reset(ii);
        colsiz_[jj] = 2;
        return col[1].val;
    }
    
    // add the requested term at the end:
    size_t n = colsiz_[jj];
    
    // allocate space for new Element if necessary:
    if ( n >= colmax_[jj] )
    {
        allocateColumn(jj, n+1);
        col = column_[jj];
    }
    
    assert_true( n < colmax_[jj] );
    col[n].reset(ii);
    ++colsiz_[jj];
    
    //printColumn(jj);
    return col[n].val;
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
                //printf("SMS1 %4i %4i % .4f\n", ii, jj, a);
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


std::string SparMatSym2::what() const
{
    std::ostringstream msg;
#if MATRIX2_USES_AVX
    msg << "SMS1x " << nbElements();
#elif MATRIX2_USES_SSE
    msg << "SMS1e " << nbElements();
#else
    msg << "SMS1 " << nbElements();
#endif
    return msg.str();
}


void SparMatSym2::printSparse(std::ostream& os, real inf) const
{
    char str[256];
    for ( size_t jj = 0; jj < size_; ++jj )
    {
        if ( colsiz_[jj] > 0 )
            os << "% column " << jj << "\n";
        for ( size_t n = 0 ; n < colsiz_[jj] ; ++n )
        {
            real v = column_[jj][n].val;
            if ( abs_real(v) >= inf )
            {
                snprintf(str, sizeof(str), "%6lu %6lu %16.6f\n", column_[jj][n].inx, jj, v);
                os << str;
            }
        }
    }
}


void SparMatSym2::printColumns(std::ostream& os)
{
    os << "SMS1 size " << size_ << ":";
    for ( size_t jj = 0; jj < size_; ++jj )
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
    os << "SMS1 col " << jj << ":";
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

        nmax_  = nbe + size_;
        ija_   = new size_t[nmax_];
        sa_    = new_real(nmax_);
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


void SparMatSym2::vecMulAddCol(const real* X, real* Y, size_t jj,
                               real const* dia, size_t start, size_t stop) const
{
    assert_true( start <= stop );
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

void SparMatSym2::vecMulAddColIso2D(const real* X, real* Y, size_t jj,
                                    real const* dia, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= 2*size_ );
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


void SparMatSym2::vecMulAddColIso3D(const real* X, real* Y, size_t jj,
                                    real const* dia, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= 3*size_ );
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

#if MATRIX2_USES_SSE

inline void multiply2(const real* X, real* Y, size_t ii,
                      const real* val, vec2 const& xx, vec2& ss)
{
    vec2 aa = loaddup2(val);
    ss = fmadd2(load2(X+ii), aa, ss);
    store2(Y+ii, fmadd2(xx, aa, load2(Y+ii)));
}


void SparMatSym2::vecMulAddColIso2D_SSE(const real* X, real* Y, size_t jj,
                                        real const* dia, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    const vec2 xx = load2(X+jj);
    vec2 ss = fmadd2(loaddup2(dia), xx, load2(Y+jj));
    // there is a dependence here for 'ss'
    for ( size_t n = start; n < stop; ++n )
        multiply2(X, Y, ija_[n], sa_+n, xx, ss);
    store2(Y+jj, ss);
}


void SparMatSym2::vecMulAddColIso2D_SSEU(const real* X, real* Y, size_t jj,
                                         real const* dia, size_t start, size_t stop) const
{
    assert_true( start <= stop );
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


void SparMatSym2::vecMulAddColIso2D_AVX(const real* X, real* Y, size_t jj,
                                        real const* dia, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    const vec4 xx = broadcast2(X+jj);  // hi position
    vec4 ss = fmadd4(broadcast1(dia), xx, broadcast2(Y+jj));
    // there is a dependence here for 'ss'
    for ( size_t n = start; n < stop; ++n )
        multiply4(X, Y, ija_[n], sa_+n, xx, ss);
    store2(Y+jj, gethi(ss));
}


void SparMatSym2::vecMulAddColIso2D_AVXU(const real* X, real* Y, size_t jj,
                                         real const* dia, size_t start, size_t stop) const
{
    assert_true( start <= stop );
    const vec4 xx = broadcast2(X+jj);  // hi and lo position
    vec4 s0 = mul4(broadcast1(dia), xx);
    vec4 s1 = broadcast2(Y+jj);
    vec4 s2 = setzero4();
    vec4 s3 = setzero4();
    
    size_t * inx = ija_ + start;
    const real * val = sa_ + start;
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
    end = sa_ + stop;
#pragma nounroll
    for ( ; val < end; ++val, ++inx )
        multiply4(X, Y, inx[0], val, xx, s0);
    store2(Y+jj, gethi(s0));
}

#endif
#endif

//------------------------------------------------------------------------------
#pragma mark - 3D SIMD

#if MATRIX2_USES_SSE
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
        vecMulAddColIso2D_AVXU(X, Y, 2*jj, sa_+jj, ija_[jj], ija_[jj+1]);
#  elif MATRIX2_USES_SSE
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
        vecMulAddColIso3D(X, Y, 3*jj, sa_+jj, ija_[jj], ija_[jj+1]);
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
