// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "matfull.h"
#include "blas.h"

#include <iomanip>
#include <sstream>


MatrixFull::MatrixFull()
{
    allo_ = 0;
    size_ = 0;
    mat_  = nullptr;
}


void MatrixFull::deallocate()
{
    free_real(mat_);
    mat_  = nullptr;
    allo_ = 0;
    size_ = 0;
}


void MatrixFull::allocate(size_t alc)
{
    assert_true( alc > 0 );
    if ( alc > allo_ )
    {
        constexpr size_t chunk = 4 * sizeof(double);
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );
        
        //printf("new block-matrix sz %i\n", nbblock );
        real* ptr = new_real(alc*alc);
        
        size_t i = 0;
        if ( mat_ )
        {
            for ( i = 0; i < size_*allo_; ++i )
                ptr[i] = mat_[i];
            free_real(mat_);
        }
        
        for ( ; i < alc*alc; ++i )
            ptr[i] = 0;
        
        mat_  = ptr;
        allo_ = alc;
    }
}

//------------------------------------------------------------------------------
#pragma mark - Access

real* MatrixFull::addr(size_t i, size_t j) const
{
    assert_true( i < size_ );
    assert_true( j < size_ );
    
    return & mat_[ i*allo_ + j ];
}


void MatrixFull::reset()
{
    for ( size_t i = 0; i < size_*allo_ ; ++i )
        mat_[i] = 0;
}


void MatrixFull::scale(const real a)
{
    for ( size_t i = 0; i < size_*allo_ ; ++i )
        mat_[i] *= a;
}

//------------------------------------------------------------------------------
#pragma mark - Vector Multiplication

void MatrixFull::vecMulAdd(const real* X, real* Y)  const
{
    for ( size_t i = 0; i < size_; ++i )
    {
        real val = 0;
        for ( size_t j = 0; j < size_; ++j )
            val += value(i,j) * X[j];
        Y[i] += val;
    }
}


void MatrixFull::vecMul(const real* X, real* Y)  const
{
    for ( size_t i = 0; i < size_; ++i )
    {
        real val = 0;
        for ( size_t j = 0; j < size_; ++j )
            val += value(i,j) * X[j];
        Y[i] = val;
    }
}


//------------------------------------------------------------------------------
#pragma mark - I/O

real MatrixFull::norm_inf() const
{
    real res = 0;
    for ( size_t i = 0; i < size_*allo_; ++i )
        res = std::max(res, abs_real(mat_[i]));
    return res;
}


void MatrixFull::print(std::ostream& os) const
{
    const int w = (int)os.width();
    os << "MatrixFull " << size_ << " [";
    for ( size_t i = 0; i < size_; ++i )
    {
        os << "\nline " << i << ":";
        for ( size_t j = 0; j < size_; ++j )
            os << " " << std::setw(w) << std::showpos << value(i, j);
    }
    os << " ]\n";
}


std::string MatrixFull::what() const
{
    std::ostringstream msg;
    msg << "mBF " << size_;
    return msg.str();
}

