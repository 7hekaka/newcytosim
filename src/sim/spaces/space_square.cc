// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_square.h"
#include "exceptions.h"
#include "mecapoint.h"
#include "iowrapper.h"
#include "glossary.h"
#include "meca.h"


SpaceSquare::SpaceSquare(SpaceProp const* p)
: Space(p)
{
    for ( size_t d = 0; d < 3; ++d )
        length_[d] = 0;
}


void SpaceSquare::resize(Glossary& opt)
{
    for ( size_t d = 0; d < DIM; ++d )
    {
        real len = length_[d];
        if ( opt.set(len, "length", d) )
            len *= 0.5;
        if ( len < 0 )
            throw InvalidParameter("square:length[] must be >= 0");
        length_[d] = len;
    }
#if ( DIM == 2 )
    // that is for impersonating a 'cylinder' in 2D:
    if ( length_[1] <= 0 )
    {
        real rad = 0;
        if ( opt.set(rad, "radius") )
            length_[1] = rad;
    }
#endif
}


void SpaceSquare::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-length_[0],-length_[1],-length_[2]);
    sup.set( length_[0], length_[1], length_[2]);
}

//------------------------------------------------------------------------------
#pragma mark - DIM=1

real SpaceSquare::volume() const
{
#if ( DIM == 1 )
    return 2 * length_[0];
#elif ( DIM == 2 )
    return 4 * length_[0] * length_[1];
#else
    return 8 * length_[0] * length_[1] * length_[2];
#endif
}

bool SpaceSquare::inside(Vector const& W) const
{
#if ( DIM == 1 )
    return abs_real(W.XX) <= length_[0];
#elif ( DIM == 2 )
    return (abs_real(W.XX) <= length_[0]) &
           (abs_real(W.YY) <= length_[1]);
#else
    return (abs_real(W.XX) <= length_[0]) &
           (abs_real(W.YY) <= length_[1]) &
           (abs_real(W.ZZ) <= length_[2]);
#endif
}

bool SpaceSquare::allInside(Vector const& W, const real rad) const
{
    assert_true( rad >= 0 );
#if ( DIM == 1 )
    return std::max(rad-W.XX, W.XX+rad) <= length_[0];
#elif ( DIM == 2 )
    return (std::max(rad-W.XX, W.XX+rad) <= length_[0]) &
           (std::max(rad-W.YY, W.YY+rad) <= length_[1]);
#else
    return (std::max(rad-W.XX, W.XX+rad) <= length_[0]) &
           (std::max(rad-W.YY, W.YY+rad) <= length_[1]) &
           (std::max(rad-W.ZZ, W.ZZ+rad) <= length_[2]);
#endif

}


#if ( DIM == 1 )

Vector SpaceSquare::project(Vector const& W) const
{
    return Vector(std::copysign(length_[0], W.XX), 0, 0);
}

#else

Vector SpaceSquare::project(Vector const& W) const
{
    Vector P(W);
    bool in = true;
    
    if ( abs_real(P.XX) > length_[0] )
    {
        P.XX = std::copysign(length_[0], P.XX);
        in = false;
    }
    if ( abs_real(P.YY) > length_[1] )
    {
        P.YY = std::copysign(length_[1], P.YY);
        in = false;
    }
#if ( DIM > 2 )
    if ( abs_real(P.ZZ) > length_[2] )
    {
        P.ZZ = std::copysign(length_[2], P.ZZ);
        in = false;
    }
#endif

    if ( in )
    {
        // find the dimensionality corresponding to the closest face
        real d0 = length_[0] - abs_real(W.XX);
        real d1 = length_[1] - abs_real(W.YY);
#if ( DIM > 2 )
        real d2 = length_[2] - abs_real(W.ZZ);
        if ( d2 < d1 )
        {
            if ( d0 < d2 )
                P.XX = std::copysign(length_[0], W.XX);
            else
                P.ZZ = std::copysign(length_[2], W.ZZ);
        }
        else
#endif
        {
            if ( d0 < d1 )
                P.XX = std::copysign(length_[0], W.XX);
            else
                P.YY = std::copysign(length_[1], W.YY);
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

void SpaceSquare::setInteraction(const real pos[], Mecapoint const& pe, Meca& meca, real stiff, const real dim[])
{
    bool in = true;
    
    for ( size_t d = 0; d < DIM; ++d )
    {
        if ( abs_real(pos[d]) > dim[d] )
        {
            meca.addPlaneClamp(DIM*pe.matIndex()+d, std::copysign(dim[d], pos[d]), stiff);
            in = false;
        }
    }

    if ( in ) 
    {
        // find the dimensionality 'dip' corresponding to the closest face
        size_t dip = 0;
        
        real l = dim[0] - abs_real(pos[0]);
#if ( DIM > 1 )
        real u = dim[1] - abs_real(pos[1]);
        if ( u < l ) { dip = 1; l = u; };
#endif
#if ( DIM > 2 )
        u = dim[2] - abs_real(pos[2]);
        if ( u < l )  dip = 2;
#endif
        meca.addPlaneClamp(DIM*pe.matIndex()+dip, std::copysign(dim[dip], pos[dip]), stiff);
    }
}


void SpaceSquare::setInteraction(Vector const& pos, Mecapoint const& pe, Meca& meca, real stiff) const
{
    setInteraction(pos, pe, meca, stiff, length_);
}


void SpaceSquare::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca& meca, real stiff) const
{
    real dim[DIM];
    for ( size_t d = 0; d < DIM; ++d )
        dim[d] = max_real(0, length_[d] - rad);

    setInteraction(pos, pe, meca, stiff, dim);
}

//------------------------------------------------------------------------------

void SpaceSquare::write(Outputter& out) const
{
    out.put_characters("square", 16);
    out.writeUInt16(4);
    out.writeFloat(length_[0]);
    out.writeFloat(length_[1]);
    out.writeFloat(length_[2]);
    out.writeFloat(0.f);
}


void SpaceSquare::setLengths(const real len[])
{
    length_[0] = len[0];
    length_[1] = len[1];
    length_[2] = len[2];
}

void SpaceSquare::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    read_data(in, 8, len, "square");
    setLengths(len);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"
using namespace gle;

bool SpaceSquare::draw() const
{
    const real X = length_[0];
    const real Y = length_[1];
    const real Z = ( DIM > 2 ) ? length_[2] : 0;

#if ( DIM > 2 )

    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  X,  Y, -Z );
    gleVertex(  X,  Y,  Z );
    gleVertex(  X, -Y, -Z );
    gleVertex(  X, -Y,  Z );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex( -X, -Y, -Z );
    gleVertex( -X, -Y,  Z );
    gleVertex( -X,  Y, -Z );
    gleVertex( -X,  Y,  Z );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  X,  Y, -Z );
    gleVertex( -X,  Y, -Z );
    gleVertex(  X,  Y,  Z );
    gleVertex( -X,  Y,  Z );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  X, -Y,  Z );
    gleVertex( -X, -Y,  Z );
    gleVertex(  X, -Y, -Z );
    gleVertex( -X, -Y, -Z );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  X,  Y,  Z );
    gleVertex( -X,  Y,  Z );
    gleVertex(  X, -Y,  Z );
    gleVertex( -X, -Y,  Z );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  X,  Y, -Z );
    gleVertex(  X, -Y, -Z );
    gleVertex( -X,  Y, -Z );
    gleVertex( -X, -Y, -Z );
    glEnd();

    glDisable(GL_LIGHTING);
    glLineWidth(0.5);
    
    glBegin(GL_LINE_LOOP);
    gleVertex(  X,  Y, Z );
    gleVertex(  X, -Y, Z );
    gleVertex( -X, -Y, Z );
    gleVertex( -X,  Y, Z );
    glEnd();
    
    glBegin(GL_LINES);
    gleVertex(  X,  Y, -Z );
    gleVertex(  X,  Y,  Z );
    gleVertex(  X, -Y, -Z );
    gleVertex(  X, -Y,  Z );
    gleVertex( -X, -Y, -Z );
    gleVertex( -X, -Y,  Z );
    gleVertex( -X,  Y, -Z );
    gleVertex( -X,  Y,  Z );
    glEnd();

#endif
  
    glDisable(GL_LIGHTING);
    glBegin(GL_LINE_LOOP);
    gleVertex(  X,  Y, -Z );
    gleVertex(  X, -Y, -Z );
    gleVertex( -X, -Y, -Z );
    gleVertex( -X,  Y, -Z );
    glEnd();

    return true;
}

#else

bool SpaceSquare::draw() const
{
    return false;
}

#endif
