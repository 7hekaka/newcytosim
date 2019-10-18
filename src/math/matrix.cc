// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "matrix.h"
#include "assert_macro.h"
#include "cblas.h"
#include <iomanip>

//------------------------------------------------------------------------------
real Matrix::value(const size_t x, const size_t y) const
{
    real* v = addr( x, y );
    if ( v == nullptr )
        return 0;
    else
        return *v;
}

real Matrix::norm_inf() const
{
    const size_t Z = size();
    real result = 0;
    for ( size_t ii = 0; ii < Z; ++ii )
    {
        for ( size_t jj = 0; jj < Z; ++jj )
        {
            real* v = addr( ii, jj );
            if ( v  &&  ( *v > result ) )
                result = *v;
        }
    }
    return result;
}

bool Matrix::nonZero() const
{
    const size_t Z = size();
    for ( size_t ii = 0; ii < Z; ++ii )
        for ( size_t jj = 0; jj < Z; ++jj )
            if ( 0 != value( ii, jj ) )
                return true;
    return false;
}

size_t Matrix::nbElements(size_t start, size_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= size_ );
    
    size_t result = 0;
    for ( size_t jj = start; jj < stop; ++jj )
        for ( size_t ii = 0; ii < size_; ++ii )
            result += ( 0 != value( ii, jj ) );
    return result;
}

//------------------------------------------------------------------------------
void Matrix::copyBlock(real* mat, unsigned ldd, size_t sx, size_t nx, size_t sy, size_t ny) const
{
    assert_true( sx + nx < size() );
    assert_true( sy + ny < size() );
    
    for ( size_t ii = 0; ii < nx; ++ii )
    for ( size_t jj = 0; jj < ny; ++jj )
        mat[ii + ldd * jj] = value( sx + ii, sy + jj );
}


void Matrix::addDiagonalBlock(real* mat, const size_t ldd, const size_t si, const size_t nb) const
{
    assert_true( si + nb < size() );

    for ( size_t jj = 0; jj < nb; ++jj )
    for ( size_t ii = 0; ii < nb; ++ii )
        mat[ ii + ldd * jj ] += value( si + ii, si + jj );
}


void Matrix::addTriangularBlock(real* mat, const size_t ldd, const size_t si, const size_t nb, const size_t dim) const
{
    assert_true( si + nb < size() );

    for ( size_t ii = 0; ii < nb; ++ii )
    for ( size_t jj = ii; jj < nb; ++jj )
        mat[ dim*ii + ldd * dim*jj ] += value( si + ii, si + jj );
}

//------------------------------------------------------------------------------
void Matrix::vecMul( const real* X, real* Y ) const
{
    zero_real(size(), Y);
    vecMulAdd( X, Y );
}


//------------------------------------------------------------------------------
void Matrix::printFull(std::ostream& os) const
{
    char str[32];
    const size_t Z = size();
    //printf("%i %i\n", size, size);
    for ( size_t ii = 0; ii < Z; ++ii )
    {
        for ( size_t jj = 0; jj < Z; ++jj )
        {
            real * a = addr(ii,jj);
            if ( a )
            {
                snprintf(str, sizeof(str), " %9.3f", *a);
                os << str;
            }
            else
                os << "       .  ";
        }
        std::endl(os);
    }
}

void Matrix::printSparse(std::ostream& os) const
{
    char str[256];
    const size_t Z = size();
    for ( size_t ii = 0; ii < Z; ++ii )
        for ( size_t jj = 0; jj < Z; ++jj )
        {
            real * v = addr(ii, jj);
            if ( v && *v != 0 )
            {
                snprintf(str, sizeof(str), "%6lu %6lu %16.6f\n", ii, jj, *v);
                os << str;
            }
        }
}

