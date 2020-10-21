// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "interpolation4.h"
#include "interpolation.h"
#include "mecapoint.h"
#include "simul.h"
#include "meca.h"


void Interpolation4::clear()
{
    mec_ = nullptr;
    prime_ = 0;
    
    rank_ = 0;
    coef_[0] = 1.0;
    coef_[1] = 0.0;
    coef_[2] = 0.0;
    coef_[3] = 0.0;
}


/**
Set a point of index 'P' on Mecable
*/
void Interpolation4::set(Mecable const* m, size_t p)
{
    assert_true( !m || p < m->nbPoints() );
    
    mec_ = m;
    prime_ = p;
    
    rank_ = 1;
    coef_[0] = 1.0;
    coef_[1] = 0.0;
    coef_[2] = 0.0;
    coef_[3] = 0.0;
}


size_t Interpolation4::calcRank() const
{
    size_t res = 4;
    while ((res > 0) & (abs_real(coef_[res-1]) < REAL_EPSILON))
        --res;
    return res;
}

/**
Coefficient 'c' determines the position between 'P' and 'Q':
- with ( c == 0 ) the point is equal to Q
- with ( c == 1 ) the point is equal to P
.
*/
void Interpolation4::set(Mecable const* m, size_t p, size_t q, real c)
{
    assert_true(m);
    assert_true(q == p+1);
    assert_true(q < m->nbPoints());
    
    mec_ = m;
    prime_ = p;
    
    rank_ = 2;
    coef_[0] = c;
    coef_[1] = 1.0 - c;
    coef_[2] = 0.0;
    coef_[3] = 0.0;
}

/**
The Vector 'vec' determines the interpolation coefficients:
    - if ( vec == [0, 0, 0] ) this interpolates P exactly
    - if ( vec == [1, 0, 0] ) this interpolates P+1
    - if ( vec == [0, 1, 0] ) this interpolates P+2
    - if ( vec == [0, 0, 1] ) this interpolates P+3
    .
This is used when the four vertices [P ... P+3] define a orthonormal reference,
as the components of 'vec' are then simply the coordinates of the position of the
interpolation in this reference frame.
*/
void Interpolation4::set(Mecable const* m, size_t p, Vector const& vec)
{
    assert_true(m);
    
    mec_ = m;
    prime_ = p;

    coef_[1] = vec.XX;
#if ( DIM == 1 )
    coef_[2] = 0.0;
    coef_[3] = 0.0;
    coef_[0] = 1.0 - vec.XX;
#elif ( DIM == 2 )
    coef_[2] = vec.YY;
    coef_[3] = 0.0;
    coef_[0] = 1.0 - vec.XX - vec.YY;
#else
    coef_[2] = vec.YY;
    coef_[3] = vec.ZZ;
    coef_[0] = 1.0 - vec.XX - vec.YY - vec.ZZ;
#endif

    rank_ = calcRank();

    // the last point to be interpolated is ( prime_ + rank_ -1 )
    assert_true(prime_+rank_ <= m->nbPoints());
}


Vector Interpolation4::position() const
{
    assert_true(mec_);
    size_t top = std::min(rank_, mec_->nbPoints()-prime_);
    Vector res = coef_[0] * mec_->posPoint(prime_);
    for ( size_t i = 1; i < top; ++i )
        res += coef_[i] * mec_->posPoint(prime_+i);
    return res;
}


void Interpolation4::addLink(Meca& meca, Interpolation const& arg, const real weight) const
{
    assert_true(mec_);
    size_t off = mec_->matIndex() + prime_;
    size_t pts[] = { off, off+1, off+2, off+3 };
    
    switch ( rank_ )
    {
        case 0:
            break;
        case 1:
            meca.addLink1(arg, off, weight);
            break;
        case 2:
            meca.addLink2(arg, pts, coef_, weight);
            break;
        case 3:
            meca.addLink3(arg, pts, coef_, weight);
            break;
        case 4:
            meca.addLink4(arg, pts, coef_, weight);
        break;
    }
}


void Interpolation4::addLink(Meca& meca, Mecapoint const& arg, const real weight) const
{
    size_t off = mec_->matIndex() + prime_;
    size_t pts[] = { off, off+1, off+2, off+3 };
    
    switch ( rank_ )
    {
        case 0:
            break;
        case 1:
            meca.addLink(arg, Mecapoint(mec_, prime_), weight);
            break;
        case 2:
            meca.addLink2(arg, pts, coef_, weight);
            break;
        case 3:
            meca.addLink3(arg, pts, coef_, weight);
            break;
        case 4:
            meca.addLink4(arg, pts, coef_, weight);
        break;
    }
}


int Interpolation4::bad() const
{
    if ( !mec_ )
        return 1;

    if ( prime_ >= mec_->nbPoints() )
        return 2;

    if ( prime_+rank_ > mec_->nbPoints() )
        return 3;

    // the sum of the coefficients should equal 1:
    real s = -1;
    for ( int d = 0; d < 4; ++d )
        s += coef_[d];
    
    if ( abs_real(s) > 0.1 )
        return 4;

    return 0;
}


void Interpolation4::write(Outputter& out) const
{
    out.writeSoftSpace();
    Object::writeReference(out, mec_);
    if ( mec_ )
    {
        out.writeUInt16(prime_);
        for ( int d = 1; d < 4; ++d )
            out.writeFloat(coef_[d]);
    }
}


void Interpolation4::read(Inputter& in, Simul& sim)
{
    ObjectTag g;
    mec_ = Simul::toMecable(sim.readReference(in, g));
    
#ifdef BACKWARD_COMPATIBILITY
    if ( mec_ || in.formatID() < 55 )
#else
    if ( mec_ )
#endif
    {
        prime_ = in.readUInt16();
        for ( int d = 1; d < 4; ++d )
            coef_[d] = in.readFloat();
        // set derived variables:
        coef_[0] = 1.0 - coef_[1] - coef_[2] - coef_[3];
        rank_ = calcRank();
    }
    else
        clear();
}


void Interpolation4::print(std::ostream& out) const
{
    const int w = 9;
    if ( mec_ )
    {
        out << "(" << mec_->reference() << "  " << std::setw(w) << coef_[0];
        for ( int d = 1; d < 4; ++d )
            out << " " << std::setw(w) << coef_[d];
        out << ")";
    }
    else
        out << "(void)";
}

