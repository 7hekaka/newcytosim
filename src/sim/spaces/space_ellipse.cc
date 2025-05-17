// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.
#include "space_ellipse.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include "project_ellipse.h"


SpaceEllipse::SpaceEllipse(SpaceProp const* p)
: Space(p)
{
#if ELLIPSE_HAS_SPHEROID
    spheroid_ = -1;
#endif
    for ( int d = 0; d < 3; ++d )
        radius_[d] = 0;
    thickness_ = 0;
}


void SpaceEllipse::update()
{
    for ( unsigned d = 0; d < DIM; ++d )
        radiusSqr_[d] = square(radius_[d]);
    
#if ( DIM > 2 ) && ELLIPSE_HAS_SPHEROID
    spheroid_ = -1;
    
    // if any two dimensions are similar, then the ellipsoid is a spheroid
    for ( int zz = 0; zz < DIM; ++zz )
    {
        int xx = ( zz + 1 ) % DIM;
        int yy = ( zz + 2 ) % DIM;
        if ( abs_real(radius_[xx]-radius_[yy]) < REAL_EPSILON*(radius_[xx]+radius_[yy]) )
            spheroid_ = zz;
    }
#endif
}


void SpaceEllipse::resize(Glossary& opt)
{
    for ( unsigned d = 0; d < DIM; ++d )
    {
        real len = radius_[d];
        if ( opt.set(len, "diameter", d) || opt.set(len, "length", d) )
            len *= 0.5;
        else opt.set(len, "radius", d);
        if ( len < REAL_EPSILON )
            throw InvalidParameter("ellipse:radius[] must be > 0");
        radius_[d] = len;
        opt.set(thickness_, "thickness");
    }
    update();
}


void SpaceEllipse::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-radius_[0],-radius_[1],-radius_[2]);
    sup.set( radius_[0], radius_[1], radius_[2]);
}

/**
 A vector orthogonal to the ellipse at position ( X, Y, Z ) is
 
    ( X / lenX^2, Y / lenY^2, Z / lenZ^2 )
 
And we need to normalize this vector
*/

Vector SpaceEllipse::normalToEdge(Vector const& pos) const
{
#if ( DIM == 1 )
    return Vector(std::copysign(real(1.0), pos.XX));
#elif ( DIM == 2 )
    return normalize(Vector(pos.XX/radiusSqr_[0], pos.YY/radiusSqr_[1]));
#else
    return normalize(Vector(pos.XX/radiusSqr_[0], pos.YY/radiusSqr_[1], pos.ZZ/radiusSqr_[2]));
#endif
}


real SpaceEllipse::volume() const
{
#if ( DIM == 1 )
    return 2 * radius_[0];
#elif ( DIM == 2 )
    return M_PI * radius_[0] * radius_[1];
#else
    constexpr real C = 4 * M_PI / 3.0;
    return (C * radius_[0]) * (radius_[1] * radius_[2]);
#endif
}


real SpaceEllipse::surface() const
{
#if ( DIM == 1 )
    return 2;
#elif ( DIM == 2 )
    // approximate formula
    real h = square(radius_[0]-radius_[1]) / square(radius_[0]+radius_[1]);
    real S = M_PI * ( radius_[0] + radius_[1] );
    return S * ( 1.0 + 0.25 * h * ( 1.0 + 0.0625 * h * ( 1.0 + 0.25 * h )));
#else
    // approximate formula
    constexpr real POW = 1.6075;
    real AB = radius_[0]*radius_[1];
    real AC = radius_[0]*radius_[2];
    real BC = radius_[1]*radius_[2];
    real S = std::pow(AB,POW) + std::pow(AC,POW) + std::pow(BC,POW);
    return (4.0*M_PI) * std::pow(S/3.0, 1.0/POW);
#endif
}


bool SpaceEllipse::inside(Vector const& W) const
{
#if ( DIM == 1 )
    return abs_real(W.XX) <= radius_[0];
#elif ( DIM == 2 )
    return square(W.XX/radius_[0]) + square(W.YY/radius_[1]) <= 1;
#else
    return square(W.XX/radius_[0]) + square(W.YY/radius_[1]) + square(W.ZZ/radius_[2]) <= 1;
#endif
}


Vector1 SpaceEllipse::project1D(Vector1 const& W) const
{
    if ( W.XX >= 0 )
        return Vector1(radius_[0], 0, 0);
    else
        return Vector1(-radius_[0], 0, 0);
}


Vector2 SpaceEllipse::project2D(Vector2 const& W) const
{
    Vector2 P(W);
    projectEllipse(P.XX, P.YY, W.XX, W.YY, radius_[0], radius_[1]);
    // check that results are valid numbers:
    assert_true(P.valid());
    return P;
}


Vector3 SpaceEllipse::project3D(Vector3 const& W) const
{
    Vector3 P(W);
#if ( DIM > 2 ) && ELLIPSE_HAS_SPHEROID
    /*
     If the ellipsoid has two equal axes, we can reduce the problem to 2D,
     because it is symmetric by rotation around the remaining axis, which
     is here indicated by 'spheroid_'.
     */
    if ( spheroid_ >= 0 )
    {
        const int zz = spheroid_;
        const int xx = ( zz + 1 ) % DIM;
        const int yy = ( zz + 2 ) % DIM;
        
        if ( radius_[xx] != radius_[yy] )
            throw InvalidParameter("Inconsistent spheroid_ dimensions");
        
        //rotate point around the xx axis to bring it into the yy-zz plane:
        real pR, rr = std::sqrt( W[xx]*W[xx] + W[yy]*W[yy] );
        projectEllipse(pR, P[zz], rr, W[zz], radius_[xx], radius_[zz]);
        // back-rotate to get the projection in 3D:
        if ( rr > 0 ) {
            real s = pR / rr;
            P[xx] = W[xx] * s;
            P[yy] = W[yy] * s;
        }
        else {
            P[xx] = 0;
            P[yy] = 0;
        }
        return P;
    }
#endif
    
    projectEllipsoid(P.data(), W.data(), radius_);
    return P;
}


//------------------------------------------------------------------------------

void SpaceEllipse::write(Outputter& out) const
{
    writeMarker(out, TAG);
    writeShape(out, "LLLE");
    out.writeUInt16(4);
    out.writeFloat(radius_[0]);
    out.writeFloat(radius_[1]);
    out.writeFloat(radius_[2]);
    out.writeFloat(thickness_);
}

void SpaceEllipse::setLengths(const real len[8])
{
    radius_[0] = len[0];
    radius_[1] = len[1];
    radius_[2] = len[2];
    thickness_ = len[3];
    update();
}

void SpaceEllipse::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "LLLE");
    setLengths(len);
}

//------------------------------------------------------------------------------
#pragma mark - OpenGL display

#ifdef DISPLAY

#include "gle.h"
#include "point_disp.h"
#include "gym_view.h"


void SpaceEllipse::draw2D(float width) const
{
    const float X(radius_[0]);
    const float Y((DIM>1)?radius_[1]:1);
    const float Z((DIM>2)?radius_[2]:1);

    gym::scale(X, Y, Z);
    gle::circle1(width);
    
    if ( prop->disp && prop->disp->visible & 2 )
    {
        gym::color(prop->disp->color2);
        gle::disc1();
    }
}


void SpaceEllipse::draw3D() const
{
    const float T = thickness_;
    const float X(T+radius_[0]);
    const float Y(T+((DIM>1)?radius_[1]:1));
    const float Z(T+((DIM>2)?radius_[2]:1));

    gym::scale(X, Y, Z);
    gle::sphere();
}

#else

void SpaceEllipse::draw2D(float) const {}
void SpaceEllipse::draw3D() const {}

#endif
