// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#include "matsym.h"
#include "assert_macro.h"
#include "blas.h"


MatrixSymmetric::MatrixSymmetric()
{
    allocated_ = 0;
    val        = nullptr;
    in_charge  = true;
}


void MatrixSymmetric::allocate(size_t alc)
{
    if ( alc > allocated_ )
    {
        allocated_ = alc;
        free_real(val);
        val = new_real(alc*alc);
    }
}


void MatrixSymmetric::deallocate()
{
    if ( in_charge )
        free_real(val);
    allocated_ = 0;
    val = nullptr;
}


void MatrixSymmetric::reset()
{
    for ( size_t i = 0; i < size_ * size_; ++i )
        val[i] = 0;
}


void MatrixSymmetric::scale( real alpha )
{
    for ( size_t i = 0; i < size_ * size_; ++i )
        val[i] *= alpha;
}

//------------------------------------------------------------------------------
real& MatrixSymmetric::operator()( size_t x, size_t y)
{
    assert_true( x < size_ );
    assert_true( y < size_ );
    return val[ std::max(x,y) + dimension_ * std::min(x,y) ];
}


real* MatrixSymmetric::addr( size_t x, size_t y) const
{
    assert_true( x < size_ );
    assert_true( y < size_ );
    return val + ( std::max(x,y) + dimension_ * std::min(x,y) );
}


bool MatrixSymmetric::notZero() const
{
    return true;
}


size_t MatrixSymmetric::nbElements(size_t start, size_t stop) const
{
    assert_true( start <= stop );
    stop = std::min(stop, size_);
    start = std::min(start, size_);

    return size_ * ( stop - start );
}


std::string MatrixSymmetric::what() const
{
    return "full-symmetric";
}

//------------------------------------------------------------------------------
void MatrixSymmetric::vecMulAdd( const real* X, real* Y ) const
{
    blas::xsymv('L', size_, 1.0, val, size_, X, 1, 1.0, Y, 1);
}


void MatrixSymmetric::vecMulAddIso2D( const real* X, real* Y ) const
{
    blas::xsymv('L', size_, 1.0, val, size_, X+0, 2, 1.0, Y+0, 2);
    blas::xsymv('L', size_, 1.0, val, size_, X+1, 2, 1.0, Y+1, 2);
}


void MatrixSymmetric::vecMulAddIso3D( const real* X, real* Y ) const
{
    blas::xsymv('L', size_, 1.0, val, size_, X+0, 3, 1.0, Y+0, 3);
    blas::xsymv('L', size_, 1.0, val, size_, X+1, 3, 1.0, Y+1, 3);
    blas::xsymv('L', size_, 1.0, val, size_, X+2, 3, 1.0, Y+2, 3);
}


