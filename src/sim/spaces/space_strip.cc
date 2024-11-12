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
    pot_ = 0;
    no_top_ = 0;
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

    opt.set(no_top_, "no_top");

    update();
}


void SpaceStrip::update()
{
#if ENABLE_PERIODIC_BOUNDARIES
    modulo_.reset();
    for ( int d = 0; d < DIM-1; ++d )
        modulo_.enablePeriodic(d, 2*half_[d]);
#endif
    mid_ = ( top_ + bot_ ) * 0.5;
    // option to limit to bottom edge:
    if ( no_top_ )
        pot_ = bot_;
    else
        pot_ = top_;
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
Vector SpaceStrip::bounce(Vector const& pos) const
{
#if ENABLE_PERIODIC_BOUNDARIES
#if ( DIM >= 3 )
    real X = modulo_.fold_(pos.XX, 0);
    real Y = modulo_.fold_(pos.YY, 1);
    real Z = bounce1(pos.ZZ, bot_, top_-bot_);
    return Vector(X, Y, Z);
#elif ( DIM > 1 )
    real X = modulo_.fold_(pos.XX, 0);
    real Y = bounce1(pos.YY, bot_, top_-bot_);
    return Vector(X, Y);
#endif
#endif
    return pos;
}


/**
 place only at top/bottom surface. This overrides the function in Space
 */
Vector SpaceStrip::placeOnEdge(real) const
{
    real Z = RNG.choice(bot_, pot_);
#if ( DIM > 2 )
    return Vector(RNG.sreal()*half_[0], RNG.sreal()*half_[0], Z);
#elif ( DIM > 1 )
    return Vector(RNG.sreal()*half_[0], Z, 0);
#else
    return Vector(Z, 0, 0);
#endif
}


/**
 directed away at the edge: up at the top surface, down at the bottom surface
 */
Vector SpaceStrip::normalToEdge(Vector const& pos) const
{
#if ( DIM > 2 )
    real Z = sign_real(pos.ZZ - mid_);
    return Vector(0, 0, Z);
#elif ( DIM > 1 )
    real Y = sign_real(pos.YY - mid_);
    return Vector(0, Y, 0);
#else
    return Vector(1, 0, 0);
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


real SpaceStrip::surface() const
{
#if ( DIM == 1 )
    return -1;
#elif ( DIM == 2 )
    return (no_top_?2:4) * half_[0];
#else
    return (no_top_?4:8) * half_[0] * half_[1];
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
    real X = sign_select(pos.XX - mid_, bot_, top_);
    return Vector(X);
#elif ( DIM == 2 )
    real Y = sign_select(pos.YY - mid_, bot_, pot_);
    return Vector(pos.XX, Y);
#else
    real Z = sign_select(pos.ZZ - mid_, bot_, pot_);
    return Vector(pos.XX, pos.YY, Z);
#endif
}


//------------------------------------------------------------------------------
#pragma mark - setConfinement


void SpaceStrip::setConfinement(Vector const& pos, Mecapoint const& mp,
                                Meca& meca, real stiff) const
{
#if ( DIM == 2 )
    real Y = sign_select(pos.YY - mid_, bot_, pot_);
    meca.addPlaneClampY(mp, Y, stiff);
#elif ( DIM > 2 )
    real Z = sign_select(pos.ZZ - mid_, bot_, pot_);
    meca.addPlaneClampZ(mp, Z, stiff);
#endif
}


void SpaceStrip::setConfinement(Vector const& pos, Mecapoint const& mp,
                                real rad, Meca& meca, real stiff) const
{
#if ( DIM == 2 )
    real Y = sign_select(pos.YY - mid_, bot_+rad, pot_-rad);
    meca.addPlaneClampY(mp, Y, stiff);
#elif ( DIM > 2 )
    real Z = sign_select(pos.ZZ - mid_, bot_+rad, pot_-rad);
    meca.addPlaneClampZ(mp, Z, stiff);
#endif
}

//------------------------------------------------------------------------------

void SpaceStrip::write(Outputter& out) const
{
    writeMarker(out, TAG);
    writeShape(out, "LLBT");
    out.writeUInt16(4);
    out.writeFloat(half_[0]);
    out.writeFloat(half_[1]);
    out.writeFloat(bot_);
    out.writeFloat(top_);
}


void SpaceStrip::setLengths(const real len[8])
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
    readShape(in, 8, len, "LLBT");
    setLengths(len);
}


//------------------------------------------------------------------------------
#pragma mark - Display

#ifdef DISPLAY

#include "gle.h"
#include "gym_flute.h"
#include "gym_draw.h"
#include "gym_view.h"
#include "gym_cap.h"

void SpaceStrip::draw2D(float width) const
{
    const float X(half_[0]);
    const float T(top_);
    const float B(bot_);

    flute2 * flu = gym::mapBufferV2(5);
    flu[0] = { X, T };
    flu[1] = { X, B };
    flu[2] = {-X, B };
    flu[3] = {-X, T };
    flu[4] = { X, T };
    gym::unmapBufferV2();
    gym::drawLines(width, 1, 4);
    gym::enableLineStipple(0x000F);
    gym::drawLines(width, 0, 4);
    gym::disableLineStipple();
}

void SpaceStrip::draw3D() const
{
    const float WIDTH = 2;
    const float X(half_[0]);
    const float Y((DIM>1) ? half_[1] : 1);
    const float T(top_);
    const float B(bot_);
    gym::shift(0, 0, 0.5*(B+T));
    gym::scale(X, Y, 0.5*(T-B));
    // draw top and bottom faces:
    gym::enableLighting();
    gle::cubeFaces();
    // draw vertical edges:
    gym::disableLighting();
    gym::enableLineStipple(0x000F);
    gle::cubeVerticalEdges(WIDTH);
    gym::disableLineStipple();
}

#else

void SpaceStrip::draw2D(float) const {}
void SpaceStrip::draw3D() const {}

#endif

