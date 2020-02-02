// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_dice.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include "meca.h"


SpaceDice::SpaceDice(SpaceProp const* p)
: Space(p)
{
    if ( DIM == 1 )
        throw InvalidParameter("dice is not usable in 1D");

    for ( int d = 0; d < 3; ++d )
        length_[d] = 0;
    edge_ = 0;
}


void SpaceDice::resize(Glossary& opt)
{
    real edg = edge_;
    
    opt.set(edg, "edge");
    if ( edg < 0 )
        throw InvalidParameter("dice:radius must be >= 0");

    for ( unsigned d = 0; d < DIM; ++d )
    {
        real len = length_[d];
        if ( opt.set(len, "length", d) )
            len *= 0.5;
        if ( len < edg )
            throw InvalidParameter("dice:length[] must be >= 2 * radius");
        length_[d] = len;
    }
    
    edge_ = edg;
    update();
}


/**
 The `dice` is included in the rectangle
 */
void SpaceDice::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-length_[0],-length_[1],-length_[2]);
    sup.set( length_[0], length_[1], length_[2]);
}


/**
 If `radius==0`, the volume should be the volume of a rectangle
 */
real SpaceDice::volume() const
{
#if ( DIM == 1 )
    return 2 * length_[0];
#elif ( DIM == 2 )
    return 4 * length_[0]*length_[1] + (M_PI-4) * square(edge_);
#else
    return 8 * length_[0]*length_[1]*length_[2]
    + 2 * (M_PI-4) * ( length_[0] + length_[1] + length_[2] - 3 * edge_ ) * square(edge_)
    + (4*M_PI/3.0 - 8) * cube(edge_);
#endif
}


//------------------------------------------------------------------------------

bool SpaceDice::inside(Vector const& w) const
{
    real dis = 0;
    for ( unsigned d = 0; d < DIM; ++d )
        dis += square(std::max((real)0, fabs(w[d]) - length_[d] + edge_));
    return ( dis <= edgeSqr_ );
}


bool SpaceDice::allInside(Vector const& w, real rad) const
{
    assert_true( rad >= 0 );
    real dis = 0;
    for ( unsigned d = 0; d < DIM; ++d )
        dis += square(std::max((real)0, fabs(w[d]) - length_[d] + edge_ + rad));
    return ( dis <= edgeSqr_ );
}

//------------------------------------------------------------------------------

#if ( DIM == 1 )

Vector SpaceDice::project(Vector const& w) const
{
    return Vector(std::copysign(length_[0], w.XX), 0, 0);
}

#else

Vector SpaceDice::project(Vector const& w) const
{
    Vector p = w;
    bool in = true;

    real X = length_[0] - fabs(w.XX);
    if ( X < edge_ ) { p.XX = std::copysign(length_[0]-edge_, w.XX); in=false; }
    
    real Y = length_[1] - fabs(w.YY);
    if ( Y < edge_ ) { p.YY = std::copysign(length_[1]-edge_, w.YY); in=false; }

#if ( DIM > 2 )
    real Z = length_[2] - fabs(w.ZZ);
    if ( Z < edge_ ) { p.ZZ = std::copysign(length_[2]-edge_, w.ZZ); in=false; }
#endif
    
    if ( in )
    {
        // find the dimensionality corresponding to the closest face
#if ( DIM > 2 )
        if ( Z < Y )
        {
            if ( Z < X )
                p.ZZ = std::copysign(length_[2], w.ZZ);
            else
                p.XX = std::copysign(length_[0], w.XX);
        }
        else
#endif
        {
            if ( Y < X )
                p.YY = std::copysign(length_[1], w.YY);
            else
                p.XX = std::copysign(length_[0], w.XX);
        }
        return p;
    }

    //normalize to radius(), and add to p to get the real projection
    real dis = edge_ / norm(w-p);
    return dis * ( w - p ) + p;
}
#endif

//------------------------------------------------------------------------------

#if ADVANCED_DICE
 
void SpaceDice::setInteraction(Vector const& w, Mecapoint const& pe, Meca& meca, real stiff, const real dim[], real E) const
{
#if ( DIM == 1 )
    meca.addPlaneClampX(pe, std::copysign(length_[0], w.XX), stiff);
#else
    real dX = length_[0] - fabs(w.XX);
    real dY = length_[1] - fabs(w.YY);
#if ( DIM > 2 )
    real dZ = length_[2] - fabs(w.ZZ);
#endif
    
#if ( DIM > 2 )
    if ( dX > edge_ && dY > edge_ && dZ > edge_ )
#else
    if ( dX > edge_ && dY > edge_ )
#endif
    {
        // find the dimensionality corresponding to the closest face
#if ( DIM > 2 )
        if ( dZ < dY )
        {
            if ( dZ < dX )
                meca.addPlaneClampZ(pe, std::copysign(length_[2], w.ZZ), stiff);
            else
                meca.addPlaneClampX(pe, std::copysign(length_[0], w.XX), stiff);
        }
        else
#endif
        {
            if ( dY < dX )
                meca.addPlaneClampY(pe, std::copysign(length_[1], w.YY), stiff);
            else
                meca.addPlaneClampX(pe, std::copysign(length_[0], w.XX), stiff);
        }
    }
#if ( DIM > 2 )
    else if ( dY > edge_ && dZ > edge_ )
    {
        meca.addPlaneClampX(pe, std::copysign(length_[0], w.XX), stiff);
    }
    else if ( dX > edge_ && dZ > edge_ )
    {
        meca.addPlaneClampY(pe, std::copysign(length_[1], w.YY), stiff);
    }
    else if ( dX > edge_ && dY > edge_ )
    {
        meca.addPlaneClampZ(pe, std::copysign(length_[2], w.ZZ), stiff);
    }
#endif
    else if ( dX > edge_ )
    {
#if ( DIM > 2 )
        real cY = std::copysign(length_[1]-edge_, w.YY);
        real cZ = std::copysign(length_[2]-edge_, w.ZZ);
        meca.addCylinderClamp(pe, Vector(1, 0, 0), Vector(0, cY, cZ), edge_, stiff);
#else
        meca.addPlaneClampY(pe, std::copysign(length_[1], w.YY), stiff);
#endif
    }
    else if ( dY > edge_ )
    {
#if ( DIM > 2 )
        real cX = std::copysign(length_[0]-edge_, w.XX);
        real cZ = std::copysign(length_[2]-edge_, w.ZZ);
        meca.addCylinderClamp(pe, Vector(0, 1, 0), Vector(cX, 0, cZ), edge_, stiff);
#else
        meca.addPlaneClampX(pe, std::copysign(length_[0], w.XX), stiff);
#endif
    }
#if ( DIM > 2 )
    else if ( dZ > edge_ )
    {
        real cX = std::copysign(length_[0]-edge_, w.XX);
        real cY = std::copysign(length_[1]-edge_, w.YY);
        meca.addCylinderClamp(pe, Vector(0, 0, 1), Vector(cX, cY, 0), edge_, stiff);
    }
#endif
    else
    {
        real cX = std::copysign(length_[0]-edge_, w.XX);
        real cY = std::copysign(length_[1]-edge_, w.YY);
#if ( DIM > 2 )
        real cZ = std::copysign(length_[2]-edge_, w.ZZ);
#else
        real cZ = 0;
#endif
        meca.addSphereClamp(pe, Vector(cX, cY, cZ), edge_, stiff);
    }
#endif
}


void SpaceDice::setInteraction(Vector const& pos, Mecapoint const& pe, Meca& meca, real stiff) const
{
    setInteraction(pos, pe, meca, stiff, length_, edge_);
}


void SpaceDice::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca& meca, real stiff) const
{
    real E = std::max((real)0, edge_-rad);  // remaining edge
    real R = std::min((real)0, edge_-rad);  // size reduction

    real dim[DIM];
    for ( unsigned d = 0; d < DIM; ++d )
        dim[d] = std::max((real)0, length_[d]+R);

    setInteraction(pos, pe, meca, stiff, dim, E);
}
#endif

//------------------------------------------------------------------------------

void SpaceDice::write(Outputter& out) const
{
    out.put_characters("dice", 16);
    out.writeUInt16(4);
    out.writeFloat(length_[0]);
    out.writeFloat(length_[1]);
    out.writeFloat(length_[2]);
    out.writeFloat(edge_);
}


void SpaceDice::setLengths(const real len[])
{
    length_[0] = len[0];
    length_[1] = len[1];
    length_[2] = len[2];
    edge_ = len[4];
    update();
}

void SpaceDice::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    read_data(in, len, "dice");
    setLengths(len);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"
using namespace gle;

bool SpaceDice::draw() const
{
#if ( DIM > 2 )
    
    const real X = length_[0] - edge_;
    const real Y = length_[1] - edge_;
    const real Z = length_[2] - edge_;
 
    const real XR = length_[0];
    const real YR = length_[1];
    const real ZR = length_[2];

    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  XR,  Y, -Z );
    gleVertex(  XR,  Y,  Z );
    gleVertex(  XR, -Y, -Z );
    gleVertex(  XR, -Y,  Z );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex( -XR, -Y, -Z );
    gleVertex( -XR, -Y,  Z );
    gleVertex( -XR,  Y, -Z );
    gleVertex( -XR,  Y,  Z );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  X,  YR, -Z );
    gleVertex( -X,  YR, -Z );
    gleVertex(  X,  YR,  Z );
    gleVertex( -X,  YR,  Z );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  X, -YR,  Z );
    gleVertex( -X, -YR,  Z );
    gleVertex(  X, -YR, -Z );
    gleVertex( -X, -YR, -Z );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  X,  Y,  ZR );
    gleVertex( -X,  Y,  ZR );
    gleVertex(  X, -Y,  ZR );
    gleVertex( -X, -Y,  ZR );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  X,  Y, -ZR );
    gleVertex(  X, -Y, -ZR );
    gleVertex( -X,  Y, -ZR );
    gleVertex( -X, -Y, -ZR );
    glEnd();
    
    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    
    glLineStipple(1, 0x000F);
    glEnable(GL_LINE_STIPPLE);
    drawSection( 0, -X, 0.01 );
    drawSection( 0,  X, 0.01 );
    drawSection( 1, -Y, 0.01 );
    drawSection( 1,  Y, 0.01 );
    drawSection( 2, -Z, 0.01 );
    drawSection( 2,  Z, 0.01 );
    glDisable(GL_LINE_STIPPLE);
    glPopAttrib();
    
#else

    drawSection( 2, 0, 0.01 );

#endif
    
    return true;
}

#else

bool SpaceDice::draw() const
{
    return false;
}

#endif


