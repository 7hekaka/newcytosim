// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#include "dim.h"
#include "space_strip.h"
#include "exceptions.h"
#include "mecapoint.h"
#include "iowrapper.h"
#include "glossary.h"
#include "meca.h"


SpaceStrip::SpaceStrip(SpaceProp const* p)
: Space(p)
{
    if ( DIM == 1 )
        throw InvalidParameter("strip is not usable in 1D");
    half_[0] = 0;
    half_[1] = 0;
    bot_ = 0;
    top_ = 0;
}


void SpaceStrip::resize(Glossary& opt)
{
    for ( unsigned d = 0; d < DIM-1; ++d )
    {
        real len = half_[d];
        if ( opt.set(len, "length", d) )
            len *= 0.5;
        if ( len < 0 )
            throw InvalidParameter("strip:length[] must be >= 0");
        half_[d] = len;
    }
    
    real bot = bot_, top = top_;
    if ( opt.set(top, "length", DIM-1) )
    {
        bot = -0.5 * top;
        top =  0.5 * top;
    }
    else
    {
        opt.set(bot, "bottom");
        opt.set(top, "top");
    }

#if ( DIM == 2 )
    // that is for impersonating a 'cylinderP' in 2D:
    if ( top <= bot )
    {
        real rad = 0;
        if ( opt.set(rad, "radius") )
        {
            top =  rad;
            bot = -rad;
        }
    }
#endif

    if ( top < bot )
        throw InvalidParameter("strip:top must be >= strip:bottom");
    
    bot_ = bot;
    top_ = top;
    
    update();
}


void SpaceStrip::update()
{
    modulo_.reset();
    for ( unsigned d = 0; d < DIM-1; ++d )
        modulo_.enable(d, 2*half_[d]);
}


void SpaceStrip::boundaries(Vector& inf, Vector& sup) const
{
#if ( DIM >= 3 )
    inf.set(-half_[0],-half_[1], bot_);
    sup.set( half_[0], half_[1], top_);
#else
    inf.set(-half_[0], bot_, 0);
    sup.set( half_[0], top_, 0);
#endif
}


/** bounce within [bot_, top_] in the last dimension, and periodic in the others */
void SpaceStrip::bounce(Vector& pos) const
{
    real W = top_ - bot_;
#if ( DIM >= 3 )
    pos.XX = fold_real(pos.XX, modulo_.period_[0]);
    pos.YY = fold_real(pos.YY, modulo_.period_[1]);
    real Z = ( pos.ZZ - bot_ ) / W;
    int i = (int)floor(Z);
    W = std::copysign(W, (i&1)?-1:1);
    pos.ZZ = bot_ + W * ( Z - ((i+1)&~1) );
#elif ( DIM > 1 )
    pos.XX = fold_real(pos.XX, modulo_.period_[0]);
    real Z = ( pos.YY - bot_ ) / W;
    int i = (int)floor(Z);
    W = std::copysign(W, (i&1)?-1:1);
    pos.YY = bot_ + W * ( Z - ((i+1)&~1) );
#endif
}


/**
 place only at top/bottom surface. This overrides the function in Space
 */
Vector SpaceStrip::placeOnEdge(real) const
{
    real Z = RNG.choice(bot_, top_);
#if ( DIM >= 3 )
    return Vector(RNG.sfloat()*half_[0], RNG.sfloat()*half_[0], Z);
#elif ( DIM > 1 )
    return Vector(RNG.sfloat()*half_[0], Z, 0);
#else
    return Vector(Z, 0, 0);
#endif
}


//------------------------------------------------------------------------------
#pragma mark -


real SpaceStrip::volume() const
{
#if ( DIM == 1 )
    return ( top_ - bot_ );
#elif ( DIM == 2 )
    return 2 * half_[0] * ( top_ - bot_ );
#else
    return 4 * half_[0] * half_[1] * ( top_ - bot_ );
#endif
}


bool SpaceStrip::inside(Vector const& pos) const
{
#if ( DIM == 1 )
    return (( bot_ <= pos.XX ) & ( pos.XX <= top_ ));
#elif ( DIM == 2 )
    return (( bot_ <= pos.YY ) & ( pos.YY <= top_ ));
#else
    return (( bot_ <= pos.ZZ ) & ( pos.ZZ <= top_ ));
#endif
}


bool SpaceStrip::allInside(Vector const& pos, const real rad) const
{
    assert_true( rad >= 0 );
#if ( DIM == 1 )
    return (( bot_+rad <= pos.XX ) & ( pos.XX+rad <= top_ ));
#elif ( DIM == 2 )
    return (( bot_+rad <= pos.YY ) & ( pos.YY+rad <= top_ ));
#else
    return (( bot_+rad <= pos.ZZ ) & ( pos.ZZ+rad <= top_ ));
#endif
}


bool SpaceStrip::allOutside(Vector const& pos, const real rad) const
{
    assert_true( rad >= 0 );
#if ( DIM == 1 )
    return (( bot_ > pos.XX+rad ) | ( pos.XX > top_+rad ));
#elif ( DIM == 2 )
    return (( bot_ > pos.YY+rad ) | ( pos.YY > top_+rad ));
#else
    return (( bot_ > pos.ZZ+rad ) | ( pos.ZZ > top_+rad ));
#endif
}


Vector SpaceStrip::project(Vector const& pos) const
{
#if ( DIM == 1 )
    real X = sign_select(2 * pos.XX - (bot_+top_), bot_, top_);
    return Vector(X);
#elif ( DIM == 2 )
    real Y = sign_select(2 * pos.YY - (bot_+top_), bot_, top_);
    return Vector(pos.XX, Y);
#else
    real Z = sign_select(2 * pos.ZZ - (bot_+top_), bot_, top_);
    return Vector(pos.XX, pos.YY, Z);
#endif
}


//------------------------------------------------------------------------------
#pragma mark - setConfinement


void SpaceStrip::setConfinement(Vector const& pos, Mecapoint const& pe, Meca& meca, real stiff) const
{
#if ( DIM == 2 )
    real Y = sign_select(2 * pos.YY - (bot_+top_), bot_, top_);
    meca.addPlaneClampY(pe, Y, stiff);
#elif ( DIM > 2 )
    real Z = sign_select(2 * pos.ZZ - (bot_+top_), bot_, top_);
    meca.addPlaneClampZ(pe, Z, stiff);
#endif
}


void SpaceStrip::setConfinement(Vector const& pos, Mecapoint const& pe, real rad, Meca& meca, real stiff) const
{
#if ( DIM == 2 )
    real Y = sign_select(2 * pos.YY - (bot_+top_), bot_+rad, top_-rad);
    meca.addPlaneClampY(pe, Y, stiff);
#elif ( DIM > 2 )
    real Z = sign_select(2 * pos.ZZ - (bot_+top_), bot_+rad, top_-rad);
    meca.addPlaneClampZ(pe, Z, stiff);
#endif
}

//------------------------------------------------------------------------------

void SpaceStrip::write(Outputter& out) const
{
    writeShape(out, "strip");
    out.writeUInt16(4);
    out.writeFloat(half_[0]);
    out.writeFloat(half_[1]);
    out.writeFloat(bot_);
    out.writeFloat(top_);
}


void SpaceStrip::setLengths(const real len[])
{
    half_[0] = len[0];
    half_[1] = len[1];
    bot_ = len[2];
    top_ = len[3];
#if BACKWARD_COMPATIBILITY < 50
    // changed from 'length[2]' to 'bot_' & 'top_' on 12.06.2020
    if ( bot_ > top_ && top_ == 0 )
    {
        top_ =  0.5 * bot_;
        bot_ = -0.5 * bot_;
    }
#endif
    update();
}


void SpaceStrip::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "strip");
    setLengths(len);
}


//------------------------------------------------------------------------------
#pragma mark - Display

#ifdef DISPLAY
#include "opengl.h"

void SpaceStrip::draw2D() const
{
    const GLfloat X(half_[0]);
    const GLfloat T(top_);
    const GLfloat B(bot_);
    
    GLfloat pts[16] = {
        -X, T, X, T, X, B,-X, B,
        +X, T, X, B,-X, T,-X, B };
    glVertexPointer(2, GL_FLOAT, 0, pts);
    glDrawArrays(GL_LINES, 0, 4);
    glLineStipple(1, 0x000F);
    glEnable(GL_LINE_STIPPLE);
    glDrawArrays(GL_LINES, 4, 4);
    glDisable(GL_LINE_STIPPLE);
}

void SpaceStrip::draw3D() const
{
    const GLfloat X(half_[0]);
    const GLfloat Y(half_[1]);
    const GLfloat T(top_);
    const GLfloat B(bot_);

    // draw faces:
    GLfloat pts[24] = {
        -X, Y, B, X, Y, B,-X,-Y, B, X,-Y, B,
        -X, Y, T,-X,-Y, T, X, Y, T, X,-Y, T };
    glVertexPointer(3, GL_FLOAT, 0, pts);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    glDrawArrays(GL_TRIANGLE_STRIP, 4, 4);
    
    // draw outline:
    GLfloat lin[30] = {
        -X, Y, B, X, Y, B, X,-Y, B,-X,-Y, B,-X, Y, B,
        -X, Y, T,-X,-Y, T, X,-Y, T, X, Y, T,-X, Y, T };
    glVertexPointer(3, GL_FLOAT, 0, lin);
    glDrawArrays(GL_LINE_STRIP, 0, 5);
    glDrawArrays(GL_LINE_STRIP, 5, 5);
    
    // draw edges on periodic boundaries:
    GLfloat edg[24] = {
        +X, Y, T, X, Y, B, X,-Y, T, X,-Y, B,
        -X, Y, T,-X, Y, B,-X,-Y, T,-X,-Y, B };
    glVertexPointer(3, GL_FLOAT, 0, edg);
    glLineStipple(1, 0x000F);
    glEnable(GL_LINE_STIPPLE);
    glDrawArrays(GL_LINES, 0, 8);
    glDisable(GL_LINE_STIPPLE);
}

#else

void SpaceStrip::draw2D() const {}
void SpaceStrip::draw3D() const {}

#endif

