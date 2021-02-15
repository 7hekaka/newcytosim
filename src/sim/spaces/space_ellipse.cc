// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_ellipse.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include "project_ellipse.h"


SpaceEllipse::SpaceEllipse(SpaceProp const* p)
: Space(p)
{
#ifdef ELLIPSE_HAS_SPHEROID
    mSpheroid = -1;
#endif
    for ( int d = 0; d < 3; ++d )
        length_[d] = 0;
}


void SpaceEllipse::update()
{
    for ( unsigned d = 0; d < DIM; ++d )
        lengthSqr_[d] = square(length_[d]);
    
#if ( DIM > 2 ) && defined ELLIPSE_HAS_SPHEROID
    mSpheroid = -1;
    
    // if any two dimensions are similar, then the ellipsoid is a spheroid
    for ( int zz = 0; zz < DIM; ++zz )
    {
        int xx = ( zz + 1 ) % DIM;
        int yy = ( zz + 2 ) % DIM;
        if ( abs_real( (length(xx)-length(yy)) / (length(xx)+length(yy)) ) < REAL_EPSILON )
            mSpheroid = zz;
    }
#endif
}


void SpaceEllipse::resize(Glossary& opt)
{
    for ( unsigned d = 0; d < DIM; ++d )
    {
        real len = length_[d];
        if ( opt.set(len, "diameter", d) || opt.set(len, "length", d) )
            len *= 0.5;
        else opt.set(len, "radius", d);
        if ( len < REAL_EPSILON )
            throw InvalidParameter("ellipse:radius[] must be > 0");
        length_[d] = len;
    }
    update();
}


void SpaceEllipse::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-length_[0],-length_[1],-length_[2]);
    sup.set( length_[0], length_[1], length_[2]);
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
    return normalize(Vector(pos.XX/lengthSqr_[0], pos.YY/lengthSqr_[1]));
#else
    return normalize(Vector(pos.XX/lengthSqr_[0], pos.YY/lengthSqr_[1], pos.ZZ/lengthSqr_[2]));
#endif
}


real SpaceEllipse::volume() const
{
#if ( DIM == 1 )
    return 2 * length_[0];
#elif ( DIM == 2 )
    return M_PI * length_[0] * length_[1];
#else
    constexpr real C = 4 * M_PI / 3.0;
    return (C * length_[0]) * (length_[1] * length_[2]);
#endif
}


real SpaceEllipse::surface() const
{
#if ( DIM == 1 )
    return 2;
#elif ( DIM == 2 )
    // approximate formula
    real h = square(length_[0]-length_[1]) / square(length_[0]+length_[1]);
    real S = M_PI * ( length_[0] + length_[1] );
    return S * ( 1.0 + 0.25 * h * ( 1.0 + 0.0625 * h * ( 1.0 + 0.25 * h )));
#else
    // approximate formula
    constexpr real POW = 1.6075;
    real AB = length_[0]*length_[1];
    real AC = length_[0]*length_[2];
    real BC = length_[1]*length_[2];
    real S = std::pow(AB,POW) + std::pow(AC,POW) + std::pow(BC,POW);
    return (4.0*M_PI) * std::pow(S/3.0, 1.0/POW);
#endif
}


bool SpaceEllipse::inside(Vector const& W) const
{
#if ( DIM == 1 )
    return abs_real(W.XX) < length_[0];
#elif ( DIM == 2 )
    return square(W.XX/length_[0]) + square(W.YY/length_[1]) <= 1;
#else
    return square(W.XX/length_[0]) + square(W.YY/length_[1]) + square(W.ZZ/length_[2]) <= 1;
#endif
}


Vector1 SpaceEllipse::project1D(Vector1 const& W) const
{
    if ( W.XX >= 0 )
        return Vector1(length_[0], 0, 0);
    else
        return Vector1(-length_[0], 0, 0);
}


Vector2 SpaceEllipse::project2D(Vector2 const& W) const
{
    Vector2 P(W);
    projectEllipse(P.XX, P.YY, W.XX, W.YY, length_[0], length_[1]);
    // check that results are valid numbers:
    assert_true(P.valid());
    return P;
}


Vector3 SpaceEllipse::project3D(Vector3 const& W) const
{
    Vector3 P(W);
#if ( DIM > 2 ) && defined ELLIPSE_HAS_SPHEROID
    /*
     If the ellipsoid has two equal axes, we can reduce the problem to 2D,
     because it is symmetric by rotation around the remaining axis, which
     is here indicated by 'mSpheroid'.
     */
    if ( mSpheroid >= 0 )
    {
        const int zz = mSpheroid;
        const int xx = ( zz + 1 ) % DIM;
        const int yy = ( zz + 2 ) % DIM;
        
        if ( length(xx) != length(yy) )
            throw InvalidParameter("Inconsistent mSpheroid dimensions");
        
        //rotate point around the xx axis to bring it into the yy-zz plane:
        real pR, rr = std::sqrt( W[xx]*W[xx] + W[yy]*W[yy] );
        projectEllipse(pR, P[zz], rr, W[zz], length(xx), length(zz), 8*REAL_EPSILON);
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
        return;
    }
#endif
    
    projectEllipsoid(P.data(), W.data(), length_);
    return P;
}


//------------------------------------------------------------------------------

void SpaceEllipse::write(Outputter& out) const
{
    writeShape(out, "ellipse");
    out.writeUInt16(4);
    out.writeFloat(length_[0]);
    out.writeFloat(length_[1]);
    out.writeFloat(length_[2]);
    out.writeFloat(0.f);
}

void SpaceEllipse::setLengths(const real len[])
{
    length_[0] = len[0];
    length_[1] = len[1];
    length_[2] = len[2];
    update();
}

void SpaceEllipse::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "ellipse");
    setLengths(len);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY

#include "gle.h"

void SpaceEllipse::draw2D() const
{
    GLfloat X(length_[0]);
    GLfloat Y((DIM>1)?length_[1]:1);
    GLfloat Z((DIM>2)?length_[2]:1);

    glPushMatrix();
    glScalef(X, Y, Z);
    gle::circle();
    glPopMatrix();
}


void SpaceEllipse::draw3D() const
{
    GLfloat X(length_[0]);
    GLfloat Y((DIM>1)?length_[1]:1);
    GLfloat Z((DIM>2)?length_[2]:1);
#if 1
    glPushMatrix();
    glScalef(X, Y, Z);
    gle::sphere8();
    glPopMatrix();
#else
    gle::ellipse(X, Y, Z);
    /* Add decoration */
    glLineWidth(0.5);
    for ( GLfloat u = -0.9f; u < 1.0f; u += 0.2f )
        gle::ellipse_circle(X, Y, Z, u);
#endif
}

#else

void SpaceEllipse::draw2D() const {}
void SpaceEllipse::draw3D() const {}

#endif
