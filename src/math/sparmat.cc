// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.


#include "real.h"
#include "sparmat.h"
#include "assert_macro.h"
#include "blas.h"
#include <iomanip>
#include <sstream>

#define SPARMAT_OPTIMIZED_MULTIPLY 1

const size_t AVAILABLE_CELL = ~0U >> 1;
const size_t LAST_IN_COLUMN = ~0U;


SparMat::SparMat()
{
    size_  = 0;
    alloc_ = 0;
    mxCol  = nullptr;
    mxRow  = nullptr;
}


void SparMat::allocate(size_t sz)
{
    size_ = sz;
    if ( size_ > alloc_ )
    {
        real   ** mxCol_new = new real*[size_];
        size_t ** mxRow_new = new size_t*[size_];
        
        size_t ii = 0;
        if ( mxCol )
        {
            for ( ; ii < alloc_; ++ii )
            {
                mxCol_new[ii] = mxCol[ii];
                mxRow_new[ii] = mxRow[ii];
            }
            delete[] mxCol;
            delete[] mxRow;
        }
        
        for ( ; ii < size_; ++ii )
        {
            mxCol_new[ii] = nullptr;
            mxRow_new[ii] = nullptr;
        }
        
        mxCol = mxCol_new;
        mxRow = mxRow_new;
        alloc_ = size_;
    }
}


void SparMat::deallocate()
{
    if ( mxCol )
    {
        for ( size_t ii = 0; ii < alloc_; ++ii )
            if ( mxCol[ii] )
            {
                delete[] mxCol[ii];
                delete[] mxRow[ii];
            };
        delete[] mxCol;
        delete[] mxRow;
        mxCol = nullptr;
        mxRow = nullptr;
    }
    alloc_ = 0;
}


void SparMat::allocateColumn( const size_t jj, size_t sz )
{
    assert_true( jj < size_ );
    assert_true( sz > 0 );
    //printf("new S-COL %i %i\n", jj, sz );
    
    constexpr size_t chunk = 16;
    sz = ( sz + chunk - 1 ) & ~( chunk -1 );

    real  * mxCol_new = new real[sz];
    size_t* mxRow_new = new size_t[sz];
    
    size_t ii = 0;
    if ( mxCol[jj] )
    {
        for ( ; mxRow[jj][ii] != LAST_IN_COLUMN ; ++ii )
        {
            mxCol_new[ii] =  mxCol[jj][ii];
            mxRow_new[ii] =  mxRow[jj][ii];
        }
        
        delete[] mxCol[jj];
        delete[] mxRow[jj];
    }
    for ( ; ii < sz-1 ; ++ii )
        mxRow_new[ii] = AVAILABLE_CELL;
    mxRow_new[sz-1] = LAST_IN_COLUMN;
    
    mxCol[jj]  = mxCol_new;
    mxRow[jj]  = mxRow_new;
}


//allocate the position if necessary:
real& SparMat::operator()( size_t x, size_t y )
{
    assert_true( x < size_ );
    assert_true( y < size_ );
    
    if ( mxRow[y] )
    {
        size_t ii = 0;
        for ( ; mxRow[y][ii] != LAST_IN_COLUMN; ++ii )
            if ( mxRow[y][ii] == x )
                return mxCol[y][ii];
        
        if ( mxRow[y][ii] == LAST_IN_COLUMN )
            allocateColumn( y, ii + 1 );
        assert_true( mxRow[y][ii] == AVAILABLE_CELL );
        mxRow[y][ii] = x;
        mxCol[y][ii] = 0;
        //printf("allo. %3i %3i\n", x, y );
        return mxCol[y][ii];
    }
    
    allocateColumn( y, 1 );
    //printf("allo. %3i %3i\n", nx, ny );
    assert_true( mxRow[y][0] == AVAILABLE_CELL );
    
    //put the diagonal term first:
    mxRow[y][0] = y;
    mxCol[y][0] = 0;
    if ( x == y )
        return mxCol[y][0];
    
    mxRow[y][1] = x;
    mxCol[y][1] = 0;
    return mxCol[y][1];
}


//does not allocate the position:
real* SparMat::address( size_t x, size_t y) const
{
    size_t * row = mxRow[y];
    if ( row )
    {
        for ( ; *row != LAST_IN_COLUMN; ++row )
            if ( *row == x )
                return & mxCol[y][ row - mxRow[y] ];
    }
    return nullptr;
}


void SparMat::reset()
{
    for ( size_t ii = 0; ii < size_; ++ii )
        if ( mxRow[ii] )
            for ( int jj = 0; mxRow[ii][jj] != LAST_IN_COLUMN; ++jj )
                mxRow[ii][jj] = AVAILABLE_CELL;
}


void SparMat::scale( real a )
{
    for ( size_t ii = 0; ii < size_; ++ii )
        if ( mxRow[ii] )
            for ( int jj = 0; mxRow[ii][jj] != LAST_IN_COLUMN; ++jj )
                mxCol[ii][jj] *= a;
}


void SparMat::addDiagonalBlock(real* mat, size_t ldd, size_t start, size_t cnt,
                               const size_t mul, const size_t amp) const
{
    start *= mul;
    cnt *= mul;
    assert_true( start + cnt <= size_ );
    
    for ( size_t jj = 0; jj < cnt; ++jj )
    {
        size_t* row = mxRow[jj+start];
        if ( row != nullptr )
        {
            real* col = mxCol[jj+start];
            for ( ; *row != LAST_IN_COLUMN; ++row, ++col )
            {
                if ( *row > start )
                {
                    size_t ii = *row - start;
                    if ( ii < cnt )
                    {
                        assert_true( ii <= jj );
                        mat[amp*(ii+ldd*jj)] += *col;
                        if ( ii != jj )
                            mat[amp*(jj+ldd*ii)] += *col;
                        //printf("Sp %4i %4i % .4f\n", ii, jj, a );
                    }
                }
            }
        }
    }
}


int SparMat::bad() const
{
    if ( size_ <= 0 ) return 1;
    for ( size_t jj = 0; jj < size_; ++jj )
    {
        if ( mxRow[jj] )
            for ( int ii = 0; mxRow[jj][ii] != LAST_IN_COLUMN; ++ii )
            {
                if ( mxRow[jj][ii] <  0     ) return 2;
                if ( mxRow[jj][ii] >= size_ ) return 3;
            }
    }
    return 0;
}


void SparMat::printSparse(std::ostream& os, real, size_t start, size_t stop) const
{
    stop = std::min(stop, size_);
    std::streamsize p = os.precision(8);
    os << "% SparMat size " << size_ << ":\n";
    for ( size_t jj = start; jj < stop; ++jj )
    {
        if ( mxRow[jj] )
            for ( int ii = 0; mxRow[jj][ii] != LAST_IN_COLUMN; ++ii )
            {
                os << mxRow[jj][ii] << " " << jj << " ";
                os << std::setw(16) << mxCol[jj][ii] << '\n';
            }
    }
    os.flush();
    os.precision(p);
}


bool SparMat::notZero() const
{
    for ( size_t jj = 0; jj < size_; ++jj )
        if ( mxRow[jj] )
            for ( int ii = 0; mxRow[jj][ii] != LAST_IN_COLUMN; ++ii )
                if ( mxCol[jj][ii] != 0 )
                    return true;
    return false;
}


size_t SparMat::nbElements(size_t start, size_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_);
    //all allocated elements are counted, even if the value is zero
    size_t cnt = 0;
    for ( size_t jj = start; jj < stop; ++jj )
        if ( mxRow[jj] )
            for ( size_t ii = 0; mxRow[jj][ii] != LAST_IN_COLUMN; ++ii )
                ++cnt;
    return cnt;
}


std::string SparMat::what() const
{
    std::ostringstream msg;
#if SPARMAT_OPTIMIZED_MULTIPLY
    msg << "mS+ " << nbElements();
#else
    msg << "mS " << nbElements();
#endif
    return msg.str();
}


void SparMat::vecMulAdd( const real* X, real* Y ) const
{
    size_t kk;
    
    for ( size_t jj = 0; jj < size_; ++jj )
    {
        if ( mxRow[jj] )
            for ( size_t ii = 0; ( kk = mxRow[jj][ii] ) != LAST_IN_COLUMN; ++ii )
            {
                Y[kk] += mxCol[jj][ii] * X[jj];
            }
    }
}

//------------------------------------------------------------------------------
#if ( SPARMAT_OPTIMIZED_MULTIPLY == 0 )


void SparMat::vecMulAddIso2D( const real* X, real* Y ) const
{
    for ( size_t jj = 0; jj < size_; ++jj )
    {
        if ( mxRow[jj] )
            for ( size_t ii = 0; mxRow[jj][ii] != LAST_IN_COLUMN; ++ii )
            {
                const size_t kk = 2 * mxRow[jj][ii];
                const real a = mxCol[jj][ii];
                Y[kk  ] += a * X[kk  ];
                Y[kk+1] += a * X[kk+1];
            }
    }
}


void SparMat::vecMulAddIso3D( const real* X, real* Y ) const
{
    for ( size_t jj = 0; jj < size_; ++jj )
    {
        if ( mxRow[jj] )
            for ( size_t ii = 0; mxRow[jj][ii] != LAST_IN_COLUMN; ++ii )
            {
                const size_t kk = 3 * mxRow[jj][ii];
                const real a = mxCol[jj][ii];
                Y[kk  ] += a * X[kk  ];
                Y[kk+1] += a * X[kk+1];
                Y[kk+2] += a * X[kk+2];
            }
    }
}


#else  // SPARMAT_OPTIMIZED_MULTIPLY


void SparMat::vecMulAddIso2D( const real* X, real* Y ) const
{
    for ( size_t jj = 0; jj < size_; ++jj )
    {
        size_t* row = mxRow[jj];
        if ( row != nullptr )
        {
            real* col = mxCol[jj];
            size_t ll = 2 * jj;
            
            real X1 = X[ll  ];
            real X2 = X[ll+1];
            
            while ( *row != LAST_IN_COLUMN )
            {
                size_t kk = 2 * ( *row );
                Y[kk  ] += (*col) * X1;
                Y[kk+1] += (*col) * X2;
                
                ++row;
                ++col;
            }
        }
    }
}


void SparMat::vecMulAddIso3D( const real* X, real* Y ) const
{
    for ( size_t jj = 0; jj < size_; ++jj )
    {
        size_t* row = mxRow[jj];
        if ( row != nullptr )
        {
            real* col = mxCol[jj];
            size_t ll = 3 * jj;
            
            real X1 = X[ll  ];
            real X2 = X[ll+1];
            real X3 = X[ll+2];
            
            while ( *row != LAST_IN_COLUMN )
            {
                size_t kk = 3 * ( *row );
                Y[kk  ] += (*col) * X1;
                Y[kk+1] += (*col) * X2;
                Y[kk+2] += (*col) * X3;
                
                ++row;
                ++col;
            }
        }
    }
}

#endif

