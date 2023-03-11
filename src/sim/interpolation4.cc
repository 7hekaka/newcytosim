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

    size_t R = 4;
    while ((R > 0) & (abs_real(coef_[R-1]) < REAL_EPSILON))
        --R;
    rank_ = R;
    
    // the last point to be interpolated is ( prime_ + rank_ -1 )
    if ( mec_ && prime_+rank_ > mec_->nbPoints() )
        throw InvalidParameter("out-of-range Interpolation4");
}

/// (position of interpolation) - (position of prime point = center of sphere)
Vector Interpolation4::dir() const
{
    if ( rank_ < 2 )
        return Vector::randU();
    real C[4] = { coef_[0] - real(1), coef_[1], coef_[2], coef_[3] };
    return mec_->interpolatePoints(prime_, C, rank_);
}

/**
Set a point of index 'P' on Mecable
*/
void Interpolation4::set(Mecable const* m, size_t p)
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
void Interpolation4::set(Mecable const* m, size_t p, size_t q, real c)
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
void Interpolation4::set(Mecable const* m, size_t p, Vector const& vec)
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
    size_t off = mec_->matIndex() + prime_;
    
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
    size_t off = mec_->matIndex() + prime_;
    
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

/**
 This creates a second link, which is offset with respect to the point defined by `this`,
 in the radial direction from the origin of the local reference system of the Solid.
 Thus together with addLink() this will align a fiber in a radial direction.
 */

void Interpolation4::addOffsetLink(Meca& meca, real len, Mecapoint const& arg, const real weight) const
{
    assert_true(mec_);
    size_t off = mec_->matIndex() + prime_;
    // coefficients to extract ( surface_point - center )
    real alp[4] = { coef_[0] - real(1), coef_[1], coef_[2], coef_[3] };
    //printf("offsetLink %9.3f %9.3f %9.3f\n", alp[1], alp[2], alp[3]);
    // extract distance in current configuration:
    real rad = mec_->interpolatePoints(prime_, alp, rank_).norm();
    // distance for distal point, offset by `len` from the bead surface
    real alpha = 1.0 + len / rad;
    // calculate coefficients to interpolate the offset point:
    real c1 = alpha * coef_[1];
    real c2 = alpha * coef_[2];
    real c3 = alpha * coef_[3];
    real c0 = ( real(1) - alp[1] ) - ( alp[2] + alp[3] );
    switch ( rank_ )
    {
        case 0:
            break;
        case 1:
            break;
        case 2:
            meca.addLink2(arg, off, c0, c1, weight);
            break;
        case 3:
            meca.addLink3(arg, off, c0, c1, c2, weight);
            break;
        case 4:
            meca.addLink4(arg, off, c0, c1, c2, c3, weight);
        break;
    }
}


/**
 This creates a second link, which is offset with respect to the point defined by `this`,
 in the direction (-1,-1,-1) relative to the local reference system of the Solid.
 Thus together with addLink() this will align a fiber with the diagonal axis of the Solid.
 */
void Interpolation4::addAlignedOffsetLink(Meca& meca, real len, Mecapoint const& arg, const real weight) const
{
    assert_true(mec_);
    size_t off = mec_->matIndex() + prime_;
    // coefficients to extract ( surface_point - center )
    real alp[4] = { coef_[0] - real(1), coef_[1], coef_[2], coef_[3] };
    // extract distance in current configuration:
    real rad = mec_->interpolatePoints(prime_, alp, rank_).norm();
    // coefficients to offset the point by 'len' along the diagonal:
    switch ( rank_ )
    {
        case 0:
            break;
        case 1:
            break;
        case 2: {
            real inc = len / rad;
            meca.addLink2(arg, off, coef_[0]-inc, coef_[1]+inc, weight);
        } break;
        case 3: {
            real inc = len / ( -M_SQRT2 * rad );
            meca.addLink3(arg, off, coef_[0]-inc*2, coef_[1]+inc, coef_[2]+inc, weight);
        } break;
        case 4: {
            real inc = len / ( -M_SQRT3 * rad );
            meca.addLink4(arg, off, coef_[0]-inc*3, coef_[1]+inc, coef_[2]+inc, coef_[3]+inc, weight);
        } break;
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

