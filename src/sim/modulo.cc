// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "modulo.h"
#include "dim.h"
#include "exceptions.h"

constexpr int PERIODIC_XYZ = 7;
constexpr int PERIODIC_XY  = 3;
constexpr int PERIODIC_X   = 1;


/// enable periodicity in dimension 'd'
void Modulo::enable(size_t d, real size)
{
    if ( size > 0 )
    {
        mMode |= 1<<d;
        period_[d] = size;
    }
    else
        ;//throw InvalidParameter("periodic:length[",d,"] must be > 0");
}


Vector Modulo::period(size_t d) const
{
    Vector vec(0,0,0);
    if ( d < DIM && ( mMode & 1<<d ))
        vec[d] = period_[d];
    return vec;
}


void Modulo::fold(Vector& vec) const
{
#if ( DIM > 2 )
    if ( mMode == PERIODIC_XYZ )
    {
        vec.XX = fold_real(vec.XX, period_[0]);
        vec.YY = fold_real(vec.YY, period_[1]);
        vec.ZZ = fold_real(vec.ZZ, period_[2]);
        return;
    }
#endif
#if ( DIM > 1 )
    if ( mMode == PERIODIC_XY )
    {
        vec.XX = fold_real(vec.XX, period_[0]);
        vec.YY = fold_real(vec.YY, period_[1]);
        return;
    }
#endif
    if ( mMode == PERIODIC_X )
    {
        vec.XX = fold_real(vec.XX, period_[0]);
        return;
    }

    if ( mMode & 1 ) vec.XX = fold_real(vec.XX, period_[0]);
#if ( DIM > 1 )
    if ( mMode & 2 ) vec.YY = fold_real(vec.YY, period_[1]);
#endif
#if ( DIM > 2 )
    if ( mMode & 4 ) vec.ZZ = fold_real(vec.ZZ, period_[2]);
#endif
}


//this makes modulo around the center 'ref'
void Modulo::fold(Vector& pos, Vector const& ref) const
{
    Vector img = pos - ref;
    fold(img);
    pos = img + ref;
}


//calculate the offset from the canonical image to actual 'pos'
Vector Modulo::offset(Vector const& pos) const
{
    Vector img = pos;
    fold(img);
    return pos - img;
}


//calculate the canonical image of 'pos' and return the associated shift
void Modulo::foldOffset(Vector& pos, Vector& off) const
{
    Vector vec = pos;
    fold(pos);
    off = vec - pos;
}


void Modulo::fold_float(float* vec) const
{
#if ( DIM > 2 )
    if ( mMode == PERIODIC_XYZ )
    {
        vec[0] = std::remainderf(vec[0], period_[0]);
        vec[1] = std::remainderf(vec[1], period_[1]);
        vec[2] = std::remainderf(vec[2], period_[2]);
        return;
    }
#endif
#if ( DIM > 1 )
    if ( mMode == PERIODIC_XY )
    {
        //printf("fold(%6.2f %6.2f / %6.2f", vec[0], ref[0], period_[0]);
        vec[0] = std::remainderf(vec[0], period_[0]);
        vec[1] = std::remainderf(vec[1], period_[1]);
        //printf(" : %6.2f)\n", vec[0]);
        return;
    }
#endif
    if ( mMode == PERIODIC_X )
    {
        vec[0] = std::remainderf(vec[0], period_[0]);
        return;
    }

    if ( mMode & 1 ) vec[0] = std::remainderf(vec[0], period_[0]);
#if ( DIM > 1 )
    if ( mMode & 2 ) vec[1] = std::remainderf(vec[1], period_[1]);
#endif
#if ( DIM > 2 )
    if ( mMode & 4 ) vec[2] = std::remainderf(vec[2], period_[2]);
#endif
}


void Modulo::fold_float(float* vec, float const* ref) const
{
#if ( DIM > 2 )
    if ( mMode == PERIODIC_XYZ )
    {
        vec[0] = ref[0] + std::remainderf(vec[0]-ref[0], period_[0]);
        vec[1] = ref[1] + std::remainderf(vec[1]-ref[1], period_[1]);
        vec[2] = ref[2] + std::remainderf(vec[2]-ref[2], period_[2]);
        return;
    }
#endif
#if ( DIM > 1 )
    if ( mMode == PERIODIC_XY )
    {
        //printf("fold(%6.2f %6.2f / %6.2f", vec[0], ref[0], period_[0]);
        vec[0] = ref[0] + std::remainderf(vec[0]-ref[0], period_[0]);
        vec[1] = ref[1] + std::remainderf(vec[1]-ref[1], period_[1]);
        //printf(" : %6.2f)\n", vec[0]);
        return;
    }
#endif
    if ( mMode == PERIODIC_X )
    {
        vec[0] = ref[0] + std::remainderf(vec[0]-ref[0], period_[0]);
        return;
    }

    if ( mMode & 1 ) vec[0] = ref[0] + std::remainderf(vec[0]-ref[0], period_[0]);
#if ( DIM > 1 )
    if ( mMode & 2 ) vec[1] = ref[1] + std::remainderf(vec[1]-ref[1], period_[1]);
#endif
#if ( DIM > 2 )
    if ( mMode & 4 ) vec[2] = ref[2] + std::remainderf(vec[2]-ref[2], period_[2]);
#endif
}
