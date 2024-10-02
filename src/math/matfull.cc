// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#include "matfull.h"
#include "blas.h"
#include "simd.h"

#include <iomanip>
#include <sstream>


#ifdef __AVX__
#  define MATFULL_USES_AVX REAL_IS_DOUBLE
#else
#  define MATFULL_USES_AVX 0
#endif



MatrixFull::MatrixFull()
{
    allo_ = 0;
    size_ = 0;
    nblk_ = 0;
    mat_  = nullptr;
}


void MatrixFull::deallocate()
{
    free_real(mat_);
    mat_  = nullptr;
    allo_ = 0;
    size_ = 0;
    nblk_ = 0;
}


void MatrixFull::allocate(index_t alc)
{
    assert_true( alc > 0 );
    if ( alc > allo_ )
    {
        constexpr index_t chunk = 4;
        alc = ( alc + chunk - 1 ) & ~( chunk -1 );
        
        //printf("new block-matrix sz %i\n", nbblock );
        real* ptr = new_real(alc*alc);
        
        size_t i = 0;
        if ( mat_ )
        {
            for ( i = 0; i < allo_*allo_; ++i )
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


real* MatrixFull::address(index_t i, index_t j) const
{
    assert_true( i < size_ );
    assert_true( j < size_ );
    size_t b = ( j >> 2 ) + nblk_ * ( i >> 2 );
    size_t x = b * 16 + ( i & 3UL ) * 4 + ( j & 3UL );
    return & mat_[ x ];
    //return & mat_[ i*allo_ + j ];
}


void MatrixFull::reset(real dia, real off)
{
    for ( index_t i = 0; i < allo_*allo_ ; ++i )
        mat_[i] = off;
    for ( index_t i = 0; i < size_ ; ++i )
        *address(i,i) = dia;
}


void MatrixFull::truncate(index_t kl, index_t ku)
{
    for ( index_t j = 0; j < size_; ++j )
    {
        //zero out terms above the diagonal:
        for ( index_t i = 0; i+ku < j; ++i )
            *address(i,j) = 0;
        
        //zero out terms below the diagonal:
        for ( index_t i = j+kl+1; i < size_; ++i )
            *address(i,j) = 0;
    }
}


void MatrixFull::importMatrix(index_t size, real const* ptr, index_t lld)
{
    resize(size);
    for ( index_t i = 0; i < size; ++i )
    for ( index_t j = 0; j < size; ++j )
        *address(i,j) = ptr[i+lld*j];
}


void MatrixFull::scale(const real a)
{
    for ( index_t i = 0; i < allo_*allo_ ; ++i )
        mat_[i] *= a;
}


void MatrixFull::transpose()
{
    for ( index_t i = 0; i < size_; ++i )
    for ( index_t j = i; j < size_; ++j )
    {
        real tmp = value(i,j);
        *address(i,j) = *address(j,i);
        *address(j,i) = tmp;
    }
}

//------------------------------------------------------------------------------
#pragma mark - Vector Multiplication


void MatrixFull::vecMulAdd(const real* X, real* Y)  const
{
    for ( index_t i = 0; i < size_; ++i )
    {
        real val = 0;
        for ( index_t j = 0; j < size_; ++j )
            val += value(i,j) * X[j];
        Y[i] += val;
    }
}

void MatrixFull::vecMul0(const real* X, real* Y)  const
{
    for ( index_t i = 0; i < size_; ++i )
    {
        real val = 0;
        for ( index_t j = 0; j < size_; ++j )
            val += value(i,j) * X[j];
        Y[i] = val;
    }
}

#if MATFULL_USES_AVX

void MatrixFull::vecMul(const real* X, real* Y)  const
{
    for ( index_t i = 0; i < size_; i += 4 )
    {
        index_t end = ~3 & size_; // last multiple of 4 <= size_
        __m256i msk = makemask(size_-end);
        vec4 y0 = setzero4();
        vec4 y1 = setzero4();
        vec4 y2 = setzero4();
        vec4 y3 = setzero4();
        real const* ptr = mat_ + SB * block(i, 0);
        const index_t last = ~7 & size_; // last multiple of 8 <= size_
        {
            vec4 s0 = setzero4();
            vec4 s1 = setzero4();
            vec4 s2 = setzero4();
            vec4 s3 = setzero4();
            for ( index_t j = 0; j < last; j += 8 )
            {
                vec4 xx = loadu4(X+j);
                vec4 yy = loadu4(X+j+4);
                y0 = fmadd4(streamload4(ptr   ), xx, y0);
                y1 = fmadd4(streamload4(ptr+ 4), xx, y1);
                y2 = fmadd4(streamload4(ptr+ 8), xx, y2);
                y3 = fmadd4(streamload4(ptr+12), xx, y3);
                s0 = fmadd4(streamload4(ptr+ SB), yy, s0);
                s1 = fmadd4(streamload4(ptr+(SB+4)), yy, s1);
                s2 = fmadd4(streamload4(ptr+(SB+8)), yy, s2);
                s3 = fmadd4(streamload4(ptr+(SB+12)), yy, s3);
                ptr += 2*SB;
            }
            y0 = add4(y0, s0);
            y1 = add4(y1, s1);
            y2 = add4(y2, s2);
            y3 = add4(y3, s3);
        }
        for ( index_t j = last; j < end; j += 4 )
        {
            vec4 xx = loadu4(X+j);
            y0 = fmadd4(streamload4(ptr   ), xx, y0);
            y1 = fmadd4(streamload4(ptr+ 4), xx, y1);
            y2 = fmadd4(streamload4(ptr+ 8), xx, y2);
            y3 = fmadd4(streamload4(ptr+12), xx, y3);
            ptr += SB;
        }
        if ( end < size_ )
        {
            vec4 xx = maskload4(X+end, msk);
            y0 = fmadd4(streamload4(ptr   ), xx, y0);
            y1 = fmadd4(streamload4(ptr+ 4), xx, y1);
            y2 = fmadd4(streamload4(ptr+ 8), xx, y2);
            y3 = fmadd4(streamload4(ptr+12), xx, y3);
        }
        // sum y0 = { Y0 Y0 Y0 Y0 }, y1 = { Y1 Y1 Y1 Y1 }, y2 = { Y2 Y2 Y2 Y2 }
        y0 = add4(unpacklo4(y0, y1), unpackhi4(y0, y1));
        y2 = add4(unpacklo4(y2, y3), unpackhi4(y2, y3));
        y0 = add4(catshift2(y0, y2), blend22(y0, y2));
        maskstore4(Y+i, makemask(size_-i), y0);
    }
}


void MatrixFull::transVecMulAdd(const real* X, real* Y)  const
{
    // as we accumulate in Y[], the columns cannot be process independently
    for ( index_t i = 0; i < size_; i += 4 )
    {
        index_t end = ~3 & size_; // last multiple of 4 <= size_
        vec4 x0, x1, x2, x3;
        {
            x0 = maskload4(X+i, makemask(size_-i));
            x1 = duplo2f128(x0);
            x3 = duphi2f128(x0);
            x0 = duplo4(x1); // = broadcast1(xxx  );
            x1 = duphi4(x1); // = broadcast1(xxx+1);
            x2 = duplo4(x3); // = broadcast1(xxx+2);
            x3 = duphi4(x3); // = broadcast1(xxx+3);
        }
        real const* ptr = mat_ + SB * block(i, 0);
        for ( index_t j = 0; j < end; j += 4 )
        {
            vec4 s = fmadd4(streamload4(ptr), x0, load4(Y+j));
            s = fmadd4(streamload4(ptr+ 4), x1, s);
            s = fmadd4(streamload4(ptr+ 8), x2, s);
            s = fmadd4(streamload4(ptr+12), x3, s);
            store4(Y+j, s);
            ptr += SB;
        }
        if ( end < size_ )
        {
            __m256i msk = makemask(size_-end);
            vec4 s = maskload4(Y+end, msk);
            s = fmadd4(streamload4(ptr  ), x0, s);
            s = fmadd4(streamload4(ptr+4), x1, s);
            s = fmadd4(streamload4(ptr+8), x2, s);
            maskstore4(Y+end, msk, s);
        }
    }
}

#else

void MatrixFull::vecMul(const real* X, real* Y)  const
{
    for ( index_t i = 0; i < size_; ++i )
    {
        real val = 0;
        for ( index_t j = 0; j < size_; ++j )
            val += value(i,j) * X[j];
        Y[i] = val;
    }
}

void MatrixFull::transVecMulAdd(const real* X, real* Y)  const
{
    for ( index_t i = 0; i < size_; ++i )
    {
        real val = 0;
        for ( index_t j = 0; j < size_; ++j )
            val += value(j,i) * X[j];
        Y[i] += val;
    }
}

#endif


//------------------------------------------------------------------------------
#pragma mark - I/O

real MatrixFull::norm_inf() const
{
    real res = 0;
    for ( index_t i = 0; i < size_*allo_; ++i )
        res = max_real(res, abs_real(mat_[i]));
    return res;
}


void MatrixFull::print(std::ostream& os, index_t imin, index_t imax, index_t jmin, index_t jmax) const
{
    imax = std::min(imax, size_);
    jmax = std::min(jmax, size_);
    
    const int digits = 3;
    char str[32] = { 0 }, fmt[32] = " %4.0f";
    snprintf(fmt, sizeof(fmt), " %%%i.%if", digits+5, digits);

    os << "MatrixFull " << size_ << " (" << nblk_ << ") [";
    for ( index_t i = imin; i < imax; ++i )
    {
        os << "\n" << std::setw(2) << i;
        for ( index_t j = jmin; j < jmax; ++j )
        {
            snprintf(str, sizeof(str), fmt, value(i, j));
            os << str;
        }
    }
    os << std::noshowpos << " ]";
}


std::string MatrixFull::what() const
{
    std::ostringstream msg;
#if MATFULL_USES_AVX
    msg << "mFx " << size_ << "x" << size_;
#else
    msg << "mF " << size_ << "x" << size_;
#endif
    return msg.str();
}

