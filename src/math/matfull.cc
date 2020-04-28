// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "matfull.h"
#include "blas.h"
#include "simd.h"

#include <iomanip>
#include <sstream>


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


void MatrixFull::allocate(size_t alc)
{
    assert_true( alc > 0 );
    if ( alc > allo_ )
    {
        constexpr size_t chunk = 4;
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


real* MatrixFull::addr(size_t i, size_t j) const
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
    for ( size_t i = 0; i < allo_*allo_ ; ++i )
        mat_[i] = off;
    for ( size_t i = 0; i < size_ ; ++i )
        *addr(i,i) = dia;
}


void MatrixFull::importMatrix(size_t size, real const* ptr, size_t lld)
{
    resize(size);
    for ( size_t i = 0; i < size; ++i )
    for ( size_t j = 0; j < size; ++j )
        *addr(i,j) = ptr[i+lld*j];
}


void MatrixFull::scale(const real a)
{
    for ( size_t i = 0; i < allo_*allo_ ; ++i )
        mat_[i] *= a;
}


void MatrixFull::transpose()
{
    for ( size_t i = 0; i < size_; ++i )
    for ( size_t j = i; j < size_; ++j )
    {
        real tmp = value(i,j);
        *addr(i,j) = *addr(j,i);
        *addr(j,i) = tmp;
    }
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

void MatrixFull::vecMul0(const real* X, real* Y)  const
{
    for ( size_t i = 0; i < size_; ++i )
    {
        real val = 0;
        for ( size_t j = 0; j < size_; ++j )
            val += value(i,j) * X[j];
        Y[i] = val;
    }
}

#ifdef __AVX__

void MatrixFull::vecMul(const real* X, real* Y)  const
{
    for ( size_t i = 0; i < size_; i += 4 )
    {
        size_t end = ~3 & size_; // last multiple of 4 <= size_
        __m256i msk = makemask(size_-end);
        vec4 y0 = setzero4();
        vec4 y1 = setzero4();
        vec4 y2 = setzero4();
        vec4 y3 = setzero4();
        real const* ptr = mat_ + 16 * block(i, 0);
        const size_t last = ~7 & size_; // last multiple of 8 <= size_
        {
            vec4 s0 = setzero4();
            vec4 s1 = setzero4();
            vec4 s2 = setzero4();
            vec4 s3 = setzero4();
            for ( size_t j = 0; j < last; j += 8 )
            {
                vec4 xx = loadu4(X+j);
                vec4 yy = loadu4(X+j+4);
                y0 = fmadd4(streamload4(ptr   ), xx, y0);
                y1 = fmadd4(streamload4(ptr+ 4), xx, y1);
                y2 = fmadd4(streamload4(ptr+ 8), xx, y2);
                y3 = fmadd4(streamload4(ptr+12), xx, y3);
                s0 = fmadd4(streamload4(ptr+16), yy, s0);
                s1 = fmadd4(streamload4(ptr+20), yy, s1);
                s2 = fmadd4(streamload4(ptr+24), yy, s2);
                s3 = fmadd4(streamload4(ptr+28), yy, s3);
                ptr += 32;
            }
            y0 = add4(y0, s0);
            y1 = add4(y1, s1);
            y2 = add4(y2, s2);
            y3 = add4(y3, s3);
        }
        for ( size_t j = last; j < end; j += 4 )
        {
            vec4 xx = loadu4(X+j);
            y0 = fmadd4(streamload4(ptr   ), xx, y0);
            y1 = fmadd4(streamload4(ptr+ 4), xx, y1);
            y2 = fmadd4(streamload4(ptr+ 8), xx, y2);
            y3 = fmadd4(streamload4(ptr+12), xx, y3);
            ptr += 16;
        }
        if ( end < size_ )
        {
            vec4 xx = maskload(X+end, msk);
            y0 = fmadd4(streamload4(ptr   ), xx, y0);
            y1 = fmadd4(streamload4(ptr+ 4), xx, y1);
            y2 = fmadd4(streamload4(ptr+ 8), xx, y2);
            y3 = fmadd4(streamload4(ptr+12), xx, y3);
        }
        // sum y0 = { Y0 Y0 Y0 Y0 }, y1 = { Y1 Y1 Y1 Y1 }, y2 = { Y2 Y2 Y2 Y2 }
        y0 = add4(unpacklo4(y0, y1), unpackhi4(y0, y1));
        y2 = add4(unpacklo4(y2, y3), unpackhi4(y2, y3));
        y0 = add4(permute2f128(y0, y2, 0x21), blend4(y0, y2, 0b1100));
        maskstore(Y+i, makemask(size_-i), y0);
    }
}


void MatrixFull::transVecMulAdd(const real* X, real* Y)  const
{
    // as we accumulate in Y[], the columns cannot be process independently
    for ( size_t i = 0; i < size_; i += 4 )
    {
        size_t end = ~3 & size_; // last multiple of 4 <= size_
        vec4 x0, x1, x2, x3;
        {
            x0 = maskload(X+i, makemask(size_-i));
            //x1 = permute2f128(x0, x0, 0x00);
            //x3 = permute2f128(x0, x0, 0x11);
            x2 = permute2f128(x0, x0, 0x21);
            x1 = blend4(x0, x2, 0b1100);
            x3 = blend4(x0, x2, 0b0011);
            x0 = duplo4(x1); // = broadcast1(xxx  );
            x1 = duphi4(x1); // = broadcast1(xxx+1);
            x2 = duplo4(x3); // = broadcast1(xxx+2);
            x3 = duphi4(x3); // = broadcast1(xxx+3);
        }
        real const* ptr = mat_ + 16 * block(i, 0);
        for ( size_t j = 0; j < end; j += 4 )
        {
            vec4 s = fmadd4(streamload4(ptr), x0, load4(Y+j));
            s = fmadd4(streamload4(ptr+ 4), x1, s);
            s = fmadd4(streamload4(ptr+ 8), x2, s);
            s = fmadd4(streamload4(ptr+12), x3, s);
            store4(Y+j, s);
            ptr += 16;
        }
        if ( end < size_ )
        {
            __m256i msk = makemask(size_-end);
            vec4 s = maskload(Y+end, msk);
            s = fmadd4(streamload4(ptr  ), x0, s);
            s = fmadd4(streamload4(ptr+4), x1, s);
            s = fmadd4(streamload4(ptr+8), x2, s);
            maskstore(Y+end, msk, s);
        }
    }
}

#else

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

void MatrixFull::transVecMulAdd(const real* X, real* Y)  const
{
    for ( size_t i = 0; i < size_; ++i )
    {
        real val = 0;
        for ( size_t j = 0; j < size_; ++j )
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
    for ( size_t i = 0; i < size_*allo_; ++i )
        res = std::max(res, abs_real(mat_[i]));
    return res;
}


void MatrixFull::print(std::ostream& os) const
{
    const int w = (int)os.width();
    os << "MatrixFull " << size_ << " (" << nblk_ << ") [";
    for ( size_t i = 0; i < size_; ++i )
    {
        os << "\n   line " << std::setw(2) << i << ":";
        for ( size_t j = 0; j < size_; ++j )
            os << " " << std::setw(w) << std::showpos << value(i, j);
    }
    os << std::noshowpos << " ]";
}


std::string MatrixFull::what() const
{
    std::ostringstream msg;
#ifdef __AVX__
    msg << "mFX " << size_ << "x" << size_;
#else
    msg << "mF " << size_ << "x" << size_;
#endif
    return msg.str();
}

