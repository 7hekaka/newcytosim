// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#include "sparmatsym.h"
#include "assert_macro.h"
#include "blas.h"

#include <iomanip>
#include <sstream>


SparMatSym::SparMatSym()
{
    alloc_ = 0;
    column_ = nullptr;
    colsiz_ = nullptr;
    colmax_ = nullptr;
}


void SparMatSym::allocate(index_t alc)
{
    if ( alc > alloc_ )
    {
        constexpr index_t chunk = 64;
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );

        //fprintf(stderr, "SMS allocate matrix %lu\n", alc);
        Element ** column_new = new Element*[alc];
        index_t * colsiz_new = new index_t[alc];
        index_t * colmax_new = new index_t[alc];
        
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
    }
}


void SparMatSym::deallocate()
{
    if ( column_ )
    {
        for ( index_t ii = 0; ii < alloc_; ++ii )
            delete[] column_[ii];
        delete[] column_; column_ = nullptr;
        delete[] colsiz_; colsiz_ = nullptr;
        delete[] colmax_; colmax_ = nullptr;
    }
    alloc_ = 0;
}


/// copy `cnt` elements from `src` to `dst`
static void copy(index_t cnt, SparMatSym::Element * src, SparMatSym::Element * dst)
{
    for ( index_t ii = 0; ii < cnt; ++ii )
        dst[ii] = src[ii];
}


void SparMatSym::allocateColumn(const index_t jj, index_t alc)
{
    assert_true( jj < size_ );

    if ( alc > colmax_[jj] )
    {
        //fprintf(stderr, "SMS allocate column %i size %u\n", jj, alc);
        constexpr index_t chunk = 16;
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
        assert_true( alc == colmax_[jj] );
    }
}



real& SparMatSym::diagonal(index_t ix)
{
    assert_true( ix < size_ );
    
    Element * col;
    
    if ( colsiz_[ix] == index_t(0) )
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
 This allocates to be able to hold the matrix element if necessary
 */
real& SparMatSym::element(index_t ii, index_t jj)
{
    assert_true( ii >= jj );
    assert_true( jj < size_ );
    //fprintf(stderr, "SMS( %6i %6i )\n", i, j);

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
    
    // add the requested term at the end:
    index_t n = colsiz_[jj];

    // allocate space for new Element if necessary:
    if ( n >= colmax_[jj] )
    {
        allocateColumn(jj, n+1);
        col = column_[jj];
    }
    
    col[n].reset(ii);
    ++colsiz_[jj];
    
    //printColumn(jj);
    return col[n].val;
}


real* SparMatSym::address(index_t i, index_t j) const
{
    // swap to get ii > jj (address lower triangle)
    index_t ii = std::max(i, j);
    index_t jj = std::min(i, j);

    for ( index_t kk = 0; kk < colsiz_[jj]; ++kk )
        if ( column_[jj][kk].inx == ii )
            return &( column_[jj][kk].val );
    
    return nullptr;
}


//------------------------------------------------------------------------------
#pragma mark -

void SparMatSym::reset()
{
    for ( index_t jj = 0; jj < size_; ++jj )
        colsiz_[jj] = 0;
}


bool SparMatSym::notZero() const
{
    //check for any non-zero sparse term:
    for ( index_t jj = 0; jj < size_; ++jj )
        for ( index_t kk = 0 ; kk < colsiz_[jj] ; ++kk )
            if ( column_[jj][kk].val != 0 )
                return true;
    
    //if here, the matrix is empty
    return false;
}


void SparMatSym::scale(const real alpha)
{
    for ( index_t jj = 0; jj < size_; ++jj )
        for ( index_t n = 0; n < colsiz_[jj]; ++n )
            column_[jj][n].val *= alpha;
}


void SparMatSym::addDiagonalBlock(real* mat, const index_t ldd, index_t start,
                                  index_t cnt, const index_t mul) const
{
    assert_true( start + cnt <= size_ );
    
    for ( index_t j = 0; j < cnt; ++j )
    {
        Element * col = column_[j+start];
        real * dst = mat + ( j + ldd * j ) * mul;
        for ( index_t n = 0; n < colsiz_[j+start]; ++n )
        {
            // assuming lower triangle is stored:
            assert_true( col[n].inx > j + start );
            index_t ii = col[n].inx;
            if ( start <= ii && ii < cnt )
            {
                index_t ij = ii - ( j + start );
                //printf("SMS2 %4i %4i % .4f\n", j, ij, a);
                dst[mul*ij] += col[n].val;
                if ( ij )
                    dst[mul*ldd*ij] += col[n].val;
            }
        }
    }
}


int SparMatSym::bad() const
{
    if ( size_ <= 0 ) return 1;
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        for ( index_t kk = 0 ; kk < colsiz_[jj] ; ++kk )
        {
            if ( column_[jj][kk].inx >= size_ ) return 2;
            if ( column_[jj][kk].inx <= jj )   return 3;
        }
    }
    return 0;
}


size_t SparMatSym::nbElements(index_t start, index_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_);
    //all allocated elements are counted, even if zero
    index_t cnt = 0;
    for ( index_t jj = start; jj < stop; ++jj )
        cnt += colsiz_[jj];
    return cnt;
}

//------------------------------------------------------------------------------
#pragma mark -

std::string SparMatSym::what() const
{
    std::ostringstream msg;
    msg << "SMS " << nbElements();
    return msg.str();
}


void SparMatSym::printSparse(std::ostream& os, real, index_t start, index_t stop) const
{
    stop = std::min(stop, size_);
    std::streamsize p = os.precision(8);
    os << "% SparMatSym size " << size_ << ":\n";
    for ( index_t jj = start; jj < stop; ++jj )
    {
        for ( index_t n = 0 ; n < colsiz_[jj] ; ++n )
        {
            os << column_[jj][n].inx << " " << jj << " ";
            os << column_[jj][n].val << "\n";
        }
    }
    os.precision(p);
}


void SparMatSym::printSummary(std::ostream& os, index_t start, index_t stop)
{
    stop = std::min(stop, size_);
    os << "SMS size " << size_ << ":";
    for ( index_t jj = start; jj < stop; ++jj )
    {
        os << "\n   " << jj << "   " << colsiz_[jj];
    }
    std::endl(os);
}


void SparMatSym::printColumn(std::ostream& os, const index_t jj)
{
    Element const* col = column_[jj];
    os << "SMS col " << jj << ":";
    for ( index_t n = 0; n < colsiz_[jj]; ++n )
    {
        os << "\n" << col[n].inx << " :";
        os << " " << col[n].val;
    }
    std::endl(os);
}


//------------------------------------------------------------------------------
#pragma mark -

bool SparMatSym::prepareForMultiply(int)
{
    return true;
}


void SparMatSym::vecMulAdd(const real* X, real* Y) const
{
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        for ( index_t kk = 0 ; kk < colsiz_[jj] ; ++kk )
        {
            const index_t ii = column_[jj][kk].inx;
            const real a = column_[jj][kk].val;
            Y[ii] += a * X[jj];
            if ( ii != jj )
                Y[jj] += a * X[ii];
        }
    }
}


void SparMatSym::vecMulAddIso2D(const real* X, real* Y) const
{
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        const index_t Djj = 2 * jj;
        for ( index_t kk = 0 ; kk < colsiz_[jj] ; ++kk )
        {
            const index_t Dii = 2 * column_[jj][kk].inx;
            const real  a = column_[jj][kk].val;
            Y[Dii  ] += a * X[Djj  ];
            Y[Dii+1] += a * X[Djj+1];
            if ( Dii != Djj )
            {
                Y[Djj  ] += a * X[Dii  ];
                Y[Djj+1] += a * X[Dii+1];
            }
        }
    }
}


void SparMatSym::vecMulAddIso3D(const real* X, real* Y) const
{
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        const index_t Djj = 3 * jj;
        for ( index_t kk = 0 ; kk < colsiz_[jj] ; ++kk )
        {
            const index_t Dii = 3 * column_[jj][kk].inx;
            const real  a = column_[jj][kk].val;
            Y[Dii  ] += a * X[Djj  ];
            Y[Dii+1] += a * X[Djj+1];
            Y[Dii+2] += a * X[Djj+2];
            if ( Dii != Djj )
            {
                Y[Djj  ] += a * X[Dii  ];
                Y[Djj+1] += a * X[Dii+1];
                Y[Djj+2] += a * X[Dii+2];
            }
        }
    }
}

