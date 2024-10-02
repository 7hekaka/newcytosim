// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#include "space_square.h"
#include "exceptions.h"
#include "mecapoint.h"
#include "iowrapper.h"
#include "glossary.h"
#include "meca.h"


SpaceSquare::SpaceSquare(SpaceProp const* p)
: Space(p)
{
    for ( size_t d = 0; d < 4; ++d )
        half_[d] = 0;
}


void SpaceSquare::resize(Glossary& opt)
{
    for ( size_t d = 0; d < DIM; ++d )
    {
        real len = half_[d];
        if ( opt.set(len, "length", d) )
            len *= 0.5;
        if ( len < 0 )
            throw InvalidParameter("square:length[] must be >= 0");
        half_[d] = len;
    }
#if ( DIM == 2 )
    // that is for impersonating a 'cylinder' in 2D:
    if ( half_[1] <= 0 )
    {
        real rad = 0;
        if ( opt.set(rad, "radius") )
            half_[1] = rad;
        else if ( opt.set(rad, "diameter") )
            half_[1] = rad * 0.5;
    }
#endif
}


void SpaceSquare::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-half_[0],-half_[1],-half_[2]);
    sup.set( half_[0], half_[1], half_[2]);
}

//------------------------------------------------------------------------------
#pragma mark - DIM=1

real SpaceSquare::volume() const
{
#if ( DIM == 1 )
    return 2 * half_[0];
#elif ( DIM == 2 )
    return 4 * half_[0] * half_[1];
#else
    return 8 * half_[0] * half_[1] * half_[2];
#endif
}

real SpaceSquare::surface() const
{
#if ( DIM == 1 )
    return 2 * half_[0];
#elif ( DIM == 2 )
    return 4 * ( half_[0] + half_[1] );
#else
    return 8 * ( half_[0] * ( half_[1] + half_[2] ) + half_[1] * half_[2] );
#endif
}

bool SpaceSquare::inside(Vector const& W) const
{
#if ( DIM == 1 )
    return abs_real(W.XX) <= half_[0];
#elif ( DIM == 2 )
    return (abs_real(W.XX) <= half_[0]) &&
           (abs_real(W.YY) <= half_[1]);
#else
    return (abs_real(W.XX) <= half_[0]) &&
           (abs_real(W.YY) <= half_[1]) &&
           (abs_real(W.ZZ) <= half_[2]);
#endif
}

bool SpaceSquare::allInside(Vector const& W, const real rad) const
{
    assert_true( rad >= 0 );
#if ( DIM == 1 )
    return std::max(rad-W.XX, W.XX+rad) <= half_[0];
#elif ( DIM == 2 )
    return (std::max(rad-W.XX, W.XX+rad) <= half_[0]) &&
           (std::max(rad-W.YY, W.YY+rad) <= half_[1]);
#else
    return (std::max(rad-W.XX, W.XX+rad) <= half_[0]) &&
           (std::max(rad-W.YY, W.YY+rad) <= half_[1]) &&
           (std::max(rad-W.ZZ, W.ZZ+rad) <= half_[2]);
#endif

}


#if ( DIM == 1 )

Vector SpaceSquare::project(Vector const& W) const
{
    return Vector(std::copysign(half_[0], W.XX), 0, 0);
}

#else

Vector SpaceSquare::project(Vector const& W) const
{
    Vector P(W);
    bool in = true;
    
    if ( abs_real(P.XX) > half_[0] )
    {
        P.XX = std::copysign(half_[0], P.XX);
        in = false;
    }
    if ( abs_real(P.YY) > half_[1] )
    {
        P.YY = std::copysign(half_[1], P.YY);
        in = false;
    }
#if ( DIM > 2 )
    if ( abs_real(P.ZZ) > half_[2] )
    {
        P.ZZ = std::copysign(half_[2], P.ZZ);
        in = false;
    }
#endif

    if ( in )
    {
        // find the dimensionality corresponding to the closest face
        real d0 = half_[0] - abs_real(W.XX);
        real d1 = half_[1] - abs_real(W.YY);
#if ( DIM > 2 )
        real d2 = half_[2] - abs_real(W.ZZ);
        if ( d2 < d1 )
        {
            if ( d0 < d2 )
                P.XX = std::copysign(half_[0], W.XX);
            else
                P.ZZ = std::copysign(half_[2], W.ZZ);
        }
        else
#endif
        {
            if ( d0 < d1 )
                P.XX = std::copysign(half_[0], W.XX);
            else
                P.YY = std::copysign(half_[1], W.YY);
        }
    }
    return P;
}
#endif

//------------------------------------------------------------------------------
#pragma mark - Interaction

/// apply a force directed towards the edge of the box
/**
 When the point is in the center of the box.
 
 When a point is along the edge of the cube, the interaction
 is flat in one direction, and curved in the two others.

 */

void SpaceSquare::setConfinement(const real pos[], Mecapoint const& mp,
                                 Meca& meca, real stiff, const real dim[])
{
    bool in = true;
    
    for ( int d = 0; d < DIM; ++d )
    {
        if ( abs_real(pos[d]) > dim[d] )
        {
            meca.addPlaneClampXYZ(mp, d, std::copysign(dim[d], pos[d]), stiff);
            in = false;
        }
    }

    if ( in ) 
    {
        // find the dimensionality 'dip' corresponding to the closest face
        int dip = 0;
        
        real l = dim[0] - abs_real(pos[0]);
#if ( DIM > 1 )
        real u = dim[1] - abs_real(pos[1]);
        if ( u < l ) { dip = 1; l = u; };
#endif
#if ( DIM > 2 )
        u = dim[2] - abs_real(pos[2]);
        if ( u < l )  dip = 2;
#endif
        meca.addPlaneClampXYZ(mp, dip, std::copysign(dim[dip], pos[dip]), stiff);
    }
}


void SpaceSquare::setConfinement(Vector const& pos, Mecapoint const& mp,
                                 Meca& meca, real stiff) const
{
    setConfinement(pos, mp, meca, stiff, half_);
}


void SpaceSquare::setConfinement(Vector const& pos, Mecapoint const& mp,
                                 real rad, Meca& meca, real stiff) const
{
    real dim[DIM];
    for ( size_t d = 0; d < DIM; ++d )
        dim[d] = max_real(0, half_[d] - rad);

    setConfinement(pos, mp, meca, stiff, dim);
}

//------------------------------------------------------------------------------

void SpaceSquare::write(Outputter& out) const
{
    writeMarker(out, TAG);
    writeShape(out, "LLL");
    out.writeUInt16(4);
    out.writeFloat(half_[0]);
    out.writeFloat(half_[1]);
    out.writeFloat(half_[2]);
    out.writeFloat(0.f);
}


void SpaceSquare::setLengths(const real len[8])
{
    half_[0] = len[0];
    half_[1] = len[1];
    half_[2] = len[2];
}

void SpaceSquare::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "LLL");
    setLengths(len);
}

//------------------------------------------------------------------------------
#pragma mark - OpenGL display

#ifdef DISPLAY

#include "gle.h"
#include "gym_view.h"
#include "gym_cap.h"

void SpaceSquare::drawEdges(float width) const
{
#if ( DIM > 2 )
    gym::disableLighting();
    gle::cubeEdges(width);
    gym::restoreLighting();
#else
    const float X(half_[0]);
    const float Y(half_[1]);
    const float Z((DIM>2) ? half_[2] : 0);
    gym::scale(X, Y, Z);
    gle::square1(width);
#endif
}

void SpaceSquare::drawFaces() const
{
#if ( DIM > 2 )
    const float X(half_[0]);
    const float Y(half_[1]);
    const float Z((DIM>2) ? half_[2] : 0);
    gym::scale(X, Y, Z);
    gle::cube();
#endif
}

#else

void SpaceSquare::drawEdges(float) const {}
void SpaceSquare::drawFaces() const {}

#endif
