// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "modulo.h"
#include "dim.h"
#include "exceptions.h"

#if ENABLE_PERIODIC_BOUNDARIES
/// global Modulo object
Modulo const* modulo = nullptr;
#endif


constexpr int PERIODIC_XYZ = 7;
constexpr int PERIODIC_XY  = 3;
constexpr int PERIODIC_X   = 1;


/// enable periodicity in dimension 'd'
void Modulo::enablePeriodic(int d, real len)
{
    if ( d < 3 && len > REAL_EPSILON )
    {
        mMode |= 1 << d;
        period_[d] = len;
        inv_period_[d] = 1 / len;
    }
    else
        ;//throw InvalidParameter("periodic:length[",d,"] must be > 0");
}


Vector3 Modulo::period(int d) const
{
    Vector3 vec(0,0,0);
    if ( d < 3 && ( mMode & 1<<d ))
        vec[d] = period_[d];
    return vec;
}


void Modulo::fold(Vector& vec) const
{
    switch ( mMode )
    {
#if ( DIM > 2 )
        case PERIODIC_XYZ:
            vec.XX = fold_(vec.XX, 0);
            vec.YY = fold_(vec.YY, 1);
            vec.ZZ = fold_(vec.ZZ, 2);
            return;
#endif
#if ( DIM > 1 )
        case PERIODIC_XY:
            vec.XX = fold_(vec.XX, 0);
            vec.YY = fold_(vec.YY, 1);
            return;
#endif
        case PERIODIC_X:
            vec.XX = fold_(vec.XX, 0);
            return;
        default:
            throw InvalidParameter("periodic boundary conditions");
    }
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


//set 'pos' to its canonical image and return the associated shift
Vector Modulo::foldOffset(Vector& pos) const
{
    Vector vec = pos;
    fold(pos);
    return vec - pos;
}


void Modulo::fold_float(float* vec) const
{
    switch ( mMode )
    {
#if ( DIM > 2 )
        case PERIODIC_XYZ:
            vec[0] = foldf(vec[0], 0);
            vec[1] = foldf(vec[1], 1);
            vec[2] = foldf(vec[2], 2);
            return;
#endif
#if ( DIM > 1 )
        case PERIODIC_XY:
            vec[0] = foldf(vec[0], 0);
            vec[1] = foldf(vec[1], 1);
            return;
#endif
        case PERIODIC_X:
            vec[0] = foldf(vec[0], 0);
            return;
        default:
            throw InvalidParameter("periodic boundary conditions");
    }
}


void Modulo::fold_float(float* vec, float const* ref) const
{
    switch ( mMode )
    {
#if ( DIM > 2 )
        case PERIODIC_XYZ:
            vec[0] = ref[0] + foldf(vec[0]-ref[0], 0);
            vec[1] = ref[1] + foldf(vec[1]-ref[1], 1);
            vec[2] = ref[2] + foldf(vec[2]-ref[2], 2);
            return;
#endif
#if ( DIM > 1 )
        case PERIODIC_XY:
            vec[0] = ref[0] + foldf(vec[0]-ref[0], 0);
            vec[1] = ref[1] + foldf(vec[1]-ref[1], 1);
            return;
#endif
        case PERIODIC_X:
            vec[0] = ref[0] + foldf(vec[0]-ref[0], 0);
            return;
        default:
            throw InvalidParameter("periodic boundary conditions");
    }
}
