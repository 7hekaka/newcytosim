// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.


#include "real.h"
#include "sparmat.h"
#include "assert_macro.h"
#include "blas.h"
#include <iomanip>
#include <sstream>

#define SPARMAT_OPTIMIZED_MULTIPLY 1

constexpr index_t FREE_CELL = ~0U >> 1;
constexpr index_t COLUMN_END = ~0U;


SparMat::SparMat()
: size_(0), alloc_(0)
{
    mxCol  = nullptr;
    mxRow  = nullptr;
}


void SparMat::allocate(index_t sz)
{
    size_ = sz;
    if ( size_ > alloc_ )
    {
        real    ** mxCol_new = new real*[size_];
        index_t ** mxRow_new = new index_t*[size_];
        
        index_t ii = 0;
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
        for ( index_t ii = 0; ii < alloc_; ++ii )
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


void SparMat::allocateColumn( const index_t jj, index_t sz )
{
    assert_true( jj < size_ );
    assert_true( sz > 0 );
    //printf("new S-COL %i %i\n", jj, sz );
    
    constexpr index_t chunk = 16;
    sz = ( sz + chunk - 1 ) & ~( chunk -1 );

    real   * mxCol_new = new real[sz];
    index_t* mxRow_new = new index_t[sz];
    
    index_t ii = 0;
    if ( mxCol[jj] )
    {
        for ( ; mxRow[jj][ii] != COLUMN_END ; ++ii )
        {
            mxCol_new[ii] =  mxCol[jj][ii];
            mxRow_new[ii] =  mxRow[jj][ii];
        }
        
        delete[] mxCol[jj];
        delete[] mxRow[jj];
    }
    for ( ; ii < sz-1 ; ++ii )
        mxRow_new[ii] = FREE_CELL;
    mxRow_new[sz-1] = COLUMN_END;
    
    mxCol[jj]  = mxCol_new;
    mxRow[jj]  = mxRow_new;
}


//allocate the position if necessary:
real& SparMat::operator()(index_t x, index_t y)
{
    assert_true( x < size_ );
    assert_true( y < size_ );
    
    if ( mxRow[y] )
    {
        index_t ii = 0;
        for ( ; mxRow[y][ii] != COLUMN_END; ++ii )
            if ( mxRow[y][ii] == x )
                return mxCol[y][ii];
        
        if ( mxRow[y][ii] == COLUMN_END )
            allocateColumn( y, ii + 1 );
        assert_true( mxRow[y][ii] == FREE_CELL );
        mxRow[y][ii] = x;
        mxCol[y][ii] = 0;
        //printf("allo. %3i %3i\n", x, y );
        return mxCol[y][ii];
    }
    
    allocateColumn( y, 1 );
    //printf("allo. %3i %3i\n", nx, ny );
    assert_true( mxRow[y][0] == FREE_CELL );
    
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
real* SparMat::address(index_t x, index_t y) const
{
    index_t * row = mxRow[y];
    if ( row )
    {
        for ( ; *row != COLUMN_END; ++row )
            if ( *row == x )
                return & mxCol[y][ row - mxRow[y] ];
    }
    return nullptr;
}


void SparMat::reset()
{
    for ( index_t ii = 0; ii < size_; ++ii )
        if ( mxRow[ii] )
            for ( index_t jj = 0; mxRow[ii][jj] != COLUMN_END; ++jj )
                mxRow[ii][jj] = FREE_CELL;
}


void SparMat::scale( real a )
{
    for ( index_t ii = 0; ii < size_; ++ii )
        if ( mxRow[ii] )
            for ( index_t jj = 0; mxRow[ii][jj] != COLUMN_END; ++jj )
                mxCol[ii][jj] *= a;
}


void SparMat::addDiagonalBlock(real* mat, index_t ldd, index_t start, index_t cnt,
                               const index_t mul) const
{
    assert_true( start + cnt <= size_ );
    index_t end = start + cnt;

    for ( index_t jj = start; jj < end; ++jj )
    {
        index_t* row = mxRow[jj];
        if ( row != nullptr )
        {
            real* col = mxCol[jj];
            for ( ; *row != COLUMN_END; ++row, ++col )
            {
                if ( *row > start )
                {
                    index_t ii = *row;
                    if ( start <= ii && ii < end )
                    {
                        mat[mul*(ii+ldd*jj)] += *col;
                        if ( ii != jj )
                            mat[mul*(jj+ldd*ii)] += *col;
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
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        if ( mxRow[jj] )
            for ( index_t ii = 0; mxRow[jj][ii] != COLUMN_END; ++ii )
            {
                if ( mxRow[jj][ii] <  0     ) return 2;
                if ( mxRow[jj][ii] >= size_ ) return 3;
            }
    }
    return 0;
}


void SparMat::printSparse(std::ostream& os, real, index_t start, index_t stop) const
{
    stop = std::min(stop, size_);
    std::streamsize p = os.precision(8);
    os << "% SparMat size " << size_ << ":\n";
    for ( index_t jj = start; jj < stop; ++jj )
    {
        if ( mxRow[jj] )
            for ( index_t ii = 0; mxRow[jj][ii] != COLUMN_END; ++ii )
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
    for ( index_t jj = 0; jj < size_; ++jj )
        if ( mxRow[jj] )
            for ( index_t ii = 0; mxRow[jj][ii] != COLUMN_END; ++ii )
                if ( mxCol[jj][ii] != 0 )
                    return true;
    return false;
}


size_t SparMat::nbElements(index_t start, index_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_);
    //all allocated elements are counted, even if the value is zero
    size_t cnt = 0;
    for ( index_t jj = start; jj < stop; ++jj )
        if ( mxRow[jj] )
            for ( index_t ii = 0; mxRow[jj][ii] != COLUMN_END; ++ii )
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
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        if ( mxRow[jj] )
        {
            index_t kk;
            for ( index_t ii = 0; ( kk = mxRow[jj][ii] ) != COLUMN_END; ++ii )
            {
                Y[kk] += mxCol[jj][ii] * X[jj];
            }
        }
    }
}

//------------------------------------------------------------------------------
#if ( SPARMAT_OPTIMIZED_MULTIPLY == 0 )


void SparMat::vecMulAddIso2D( const real* X, real* Y ) const
{
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        if ( mxRow[jj] )
            for ( index_t ii = 0; mxRow[jj][ii] != COLUMN_END; ++ii )
            {
                const index_t kk = 2 * mxRow[jj][ii];
                const real a = mxCol[jj][ii];
                Y[kk  ] += a * X[kk  ];
                Y[kk+1] += a * X[kk+1];
            }
    }
}


void SparMat::vecMulAddIso3D( const real* X, real* Y ) const
{
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        if ( mxRow[jj] )
            for ( index_t ii = 0; mxRow[jj][ii] != COLUMN_END; ++ii )
            {
                const index_t kk = 3 * mxRow[jj][ii];
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
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        index_t* row = mxRow[jj];
        if ( row != nullptr )
        {
            real* col = mxCol[jj];
            index_t ll = 2 * jj;
            
            real X1 = X[ll  ];
            real X2 = X[ll+1];
            
            while ( *row != COLUMN_END )
            {
                index_t kk = 2 * ( *row );
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
    for ( index_t jj = 0; jj < size_; ++jj )
    {
        index_t* row = mxRow[jj];
        if ( row != nullptr )
        {
            real* col = mxCol[jj];
            index_t ll = 3 * jj;
            
            real X1 = X[ll  ];
            real X2 = X[ll+1];
            real X3 = X[ll+2];
            
            while ( *row != COLUMN_END )
            {
                index_t kk = 3 * ( *row );
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

