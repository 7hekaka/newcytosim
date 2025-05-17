// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University

#include "interpolation4.h"
#include "interpolation.h"
#include "mecapoint.h"
#include "simul.h"
#include "meca.h"


void Interpolation4::clear()
{
    mec_ = nullptr;
    prime_ = 0;
    set_coef(0, 0, 0);
}


void Interpolation4::polish()
{
    // ensures that sum of coefficients is 1
    coef_[0] = ( real(1) - coef_[1] ) - ( coef_[2] + coef_[3] );

    unsigned R = 4;
    while ((R > 0) & (abs_real(coef_[R-1]) < REAL_EPSILON))
        --R;
    rank_ = R;
    
    // the last point to be interpolated is ( prime_ + rank_ -1 )
    if ( mec_ && prime_+rank_ > mec_->nbPoints() )
        throw InvalidParameter("out-of-range Interpolation4");
}


/// (position of interpolation) - (position of prime point = center of sphere)
Vector Interpolation4::normal() const
{
    if ( rank_ < 2 )
        return Vector::randU();
    real C[4] = { coef_[0] - real(1), coef_[1], coef_[2], coef_[3] };
    return mec_->interpolatePoints(prime_, C, rank_);
}

/**
Set a point of index 'P' on Mecable
*/
void Interpolation4::set(Mecable const* m, index_t p)
{
    assert_true( !m || p < m->nbPoints() );
    
    mec_ = m;
    prime_ = p;
    set_coef(0, 0, 0);
}

/**
Coefficient 'c' determines the position between 'P' and 'Q':
- with ( c == 0 ) the point is equal to Q
- with ( c == 1 ) the point is equal to P
.
*/
void Interpolation4::set(Mecable const* m, index_t p, index_t q, real c)
{
    assert_true(m);
    assert_true(q < m->nbPoints());
    
    mec_ = m;
    prime_ = p;
    if ( q == p + 1 )
        set_coef(c, 0, 0);
    else if ( q == p + 2 )
        set_coef(0, c, 0);
    else if ( q == p + 3 )
        set_coef(0, 0, c);
    else
        throw InvalidSyntax("out-of-range index in Interpolation");
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
void Interpolation4::set(Mecable const* m, index_t p, Vector const& vec)
{
    assert_true(m);
    
    mec_ = m;
    prime_ = p;

#if ( DIM == 1 )
    set_coef(vec.XX, 0, 0);
#elif ( DIM == 2 )
    set_coef(vec.XX, vec.YY, 0);
#else
    set_coef(vec.XX, vec.YY, vec.ZZ);
#endif
}


void Interpolation4::addLink(Meca& meca, Mecapoint const& arg, const real weight) const
{
    index_t off = mec_->matIndex() + prime_;
    
    switch ( rank_ )
    {
        case 0:
            break;
        case 1:
            meca.addLink(arg, Mecapoint(mec_, prime_), weight);
            break;
        case 2:
            meca.addLink2(arg, off, coef_[0], coef_[1], weight);
            break;
        case 3:
            meca.addLink3(arg, off, coef_[0], coef_[1], coef_[2], weight);
            break;
        case 4:
            meca.addLink4(arg, off, coef_[0], coef_[1], coef_[2], coef_[3], weight);
        break;
    }
}


void Interpolation4::addLink(Meca& meca, Interpolation const& arg, const real weight) const
{
    assert_true(mec_);
    index_t off = mec_->matIndex() + prime_;
    
    switch ( rank_ )
    {
        case 0:
            break;
        case 1:
            meca.addLink1(arg, off, weight);
            break;
        case 2:
            meca.addLink2(arg, off, coef_[0], coef_[1], weight);
            break;
        case 3:
            meca.addLink3(arg, off, coef_[0], coef_[1], coef_[2], weight);
            break;
        case 4:
            meca.addLink4(arg, off, coef_[0], coef_[1], coef_[2], coef_[3], weight);
        break;
    }
}


int Interpolation4::invalid() const
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
    Object * obj = sim.readReference(in, g);
    mec_ = Simul::toMecable(obj);
    
#if BACKWARD_COMPATIBILITY < 55
    if ( mec_ || in.formatID() < 55 )
#else
    if ( mec_ )
#endif
    {
        prime_ = in.readUInt16();
        for ( int d = 1; d < 4; ++d )
            coef_[d] = in.readFloat();
        polish();
    }
    else
    {
        clear();
        if ( obj )
            throw InvalidIO("invalid pointer while reading Interpolation4");
    }
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

