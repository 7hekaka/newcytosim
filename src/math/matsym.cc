// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#include "matsym.h"
#include "assert_macro.h"
#include "blas.h"


MatrixSymmetric::MatrixSymmetric()
{
    allocated_ = 0;
    sym_       = nullptr;
    in_charge_ = true;
}


void MatrixSymmetric::allocate(size_t alc)
{
    if ( alc > allocated_ )
    {
        constexpr size_t chunk = 4;
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );
        allocated_ = alc;
        free_real(sym_);
        sym_ = new_real(alc*alc);
    }
}


void MatrixSymmetric::deallocate()
{
    if ( in_charge_ )
        free_real(sym_);
    allocated_ = 0;
    sym_ = nullptr;
}


void MatrixSymmetric::reset()
{
    for ( size_t i = 0; i < size_ * size_; ++i )
        sym_[i] = 0;
}


void MatrixSymmetric::scale(real alpha)
{
    for ( size_t i = 0; i < size_ * size_; ++i )
        sym_[i] *= alpha;
}

//------------------------------------------------------------------------------
real& MatrixSymmetric::operator()(size_t x, size_t y)
{
    assert_true( x < size_ );
    assert_true( y < size_ );
    return sym_[ std::max(x,y) + ldim_ * std::min(x,y) ];
}


real* MatrixSymmetric::address(size_t x, size_t y) const
{
    assert_true( x < size_ );
    assert_true( y < size_ );
    return sym_ + ( std::max(x,y) + ldim_ * std::min(x,y) );
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
    blas::xsymv('L', size_, 1.0, sym_, ldim_, X, 1, 1.0, Y, 1);
}


void MatrixSymmetric::vecMulAddIso2D( const real* X, real* Y ) const
{
    blas::xsymv('L', size_, 1.0, sym_, ldim_, X+0, 2, 1.0, Y+0, 2);
    blas::xsymv('L', size_, 1.0, sym_, ldim_, X+1, 2, 1.0, Y+1, 2);
}


void MatrixSymmetric::vecMulAddIso3D( const real* X, real* Y ) const
{
    blas::xsymv('L', size_, 1.0, sym_, ldim_, X+0, 3, 1.0, Y+0, 3);
    blas::xsymv('L', size_, 1.0, sym_, ldim_, X+1, 3, 1.0, Y+1, 3);
    blas::xsymv('L', size_, 1.0, sym_, ldim_, X+2, 3, 1.0, Y+2, 3);
}


