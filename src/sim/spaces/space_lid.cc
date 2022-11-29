// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#include "dim.h"
#include "space_lid.h"
#include "exceptions.h"
#include "mecapoint.h"
#include "iowrapper.h"
#include "glossary.h"
#include "random.h"
#include "meca.h"


SpaceLid::SpaceLid(SpaceDynamicProp const* p)
: Space(p)
{
    if ( DIM == 1 )
        throw InvalidParameter("lid  is not usable in 1D");
    half_[0] = 0;
    half_[1] = 0;
    bot_ = 0;
    top_ = 0;
    force_ = 0;
}


void SpaceLid::resize(Glossary& opt)
{
    for ( unsigned d = 0; d < DIM; ++d )
    {
        real len = half_[d];
        if ( opt.set(len, "length", d) )
            len *= 0.5;
        if ( len <= 0 )
            throw InvalidParameter("lid:length must be > 0");
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

    if ( top < bot )
        throw InvalidParameter("lid:top must be >= lid:bottom");
    
    bot_ = bot;
    top_ = top;
    
    update();
}


void SpaceLid::update()
{
    modulo_.reset();
    for ( unsigned d = 0; d < DIM-1; ++d )
        modulo_.enable(d, 2*half_[d]);
}


void SpaceLid::boundaries(Vector& inf, Vector& sup) const
{
#if ( DIM >= 3 )
    inf.set(-half_[0],-half_[1], bot_);
    sup.set( half_[0], half_[1], top_);
#else
    inf.set(-half_[0], bot_, 0);
    sup.set( half_[0], top_, 0);
#endif
}


Vector SpaceLid::bounce(Vector const& pos) const
{
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
    return pos;
}


/**
 place only at upper boundary. This overrides the function in Space
 */
Vector SpaceLid::placeOnEdge(real) const
{
    return Vector(RNG.sreal()*half_[0], top_, 0);
}


//------------------------------------------------------------------------------
#pragma mark -


real SpaceLid::volume() const
{
#if ( DIM == 1 )
    return ( top_ - bot_ );
#elif ( DIM == 2 )
    return 2 * half_[0] * ( top_ - bot_ );
#else
    return 4 * half_[0] * half_[1] * ( top_ - bot_ );
#endif
}


bool SpaceLid::inside(Vector const& pos) const
{
#if ( DIM == 1 )
    return (( bot_ <= pos.XX ) & ( pos.XX <= top_ ));
#elif ( DIM == 2 )
    return (( bot_ <= pos.YY ) & ( pos.YY <= top_ ));
#else
    return (( bot_ <= pos.ZZ ) & ( pos.ZZ <= top_ ));
#endif
}


bool SpaceLid::allInside(Vector const& pos, const real rad) const
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


bool SpaceLid::allOutside(Vector const& pos, const real rad) const
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


Vector SpaceLid::project(Vector const& pos) const
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


void SpaceLid::setConfinement(Vector const& pos, Mecapoint const& mp,
                              Meca& meca, real stiff) const
{
#if ( DIM == 2 )
    real S = 2 * pos.YY - (bot_+top_);
#elif ( DIM > 2 )
    real S = 2 * pos.ZZ - (bot_+top_);
#endif

    // record force only on top edge:
#if ( DIM == 2 )
    meca.addPlaneClampY(mp, sign_select(S, bot_, top_), stiff);
    force_ += ( S > 0 ) * stiff * ( pos.YY - top_ );
#elif ( DIM > 2 )
    meca.addPlaneClampZ(mp, sign_select(S, bot_, top_), stiff);
    force_ += ( S > 0 ) * stiff * ( pos.ZZ - top_ );
#endif
}


void SpaceLid::setConfinement(Vector const& pos, Mecapoint const& mp, real rad,
                              Meca& meca, real stiff) const
{
#if ( DIM == 2 )
    real S = 2 * pos.YY - (bot_+top_);
#elif ( DIM > 2 )
    real S = 2 * pos.ZZ - (bot_+top_);
#endif

    // record force only on top edge:
#if ( DIM == 2 )
    real Y = sign_select(S, bot_+rad, top_-rad);
    meca.addPlaneClampY(mp, Y, stiff);
    force_ += ( S > 0 ) * stiff * ( pos.YY - Y );
#elif ( DIM > 2 )
    real Z = sign_select(S, bot_+rad, top_-rad);
    meca.addPlaneClampZ(mp, Z, stiff);
    force_ += ( S > 0 ) * stiff * ( pos.ZZ - Z );
#endif
}

void SpaceLid::setInteractions(Meca& meca, Simul const&) const
{
    force_ = 0;
}

void SpaceLid::step()
{
    real dc = prop()->mobility_dt * force_;
    
    if ( abs_real(dc) < 1 )
        top_ += dc;
    else
        std::cerr << "Error: lid displacement is too fast: " << dc << '\n';
    std::cerr << "SpaceLid: top = " << top_ << " force = " << force_ << '\n';
}


//------------------------------------------------------------------------------

void SpaceLid::write(Outputter& out) const
{
    writeMarker(out, TAG);
    writeShape(out, "LLBT");
    out.writeUInt16(6);
    out.writeFloat(half_[0]);
    out.writeFloat(half_[1]);
    out.writeFloat(bot_);
    out.writeFloat(top_);
    out.writeFloat(0.f);
    out.writeFloat(force_);
}


void SpaceLid::setLengths(const real len[])
{
    half_[0] = len[0];
    half_[1] = len[1];
    bot_   = len[2];
    top_   = len[3];
    force_ = len[5];
    update();
}


void SpaceLid::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "LLBT");
    setLengths(len);
}


//------------------------------------------------------------------------------
#pragma mark - Display

#ifdef DISPLAY

#include "gym_flute.h"
#include "gym_draw.h"
#include "gym_view.h"
#include "gym_cap.h"


void SpaceLid::draw2D(float width) const
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

void SpaceLid::draw3D() const
{
    const float WIDTH = 2;
    const float X(half_[0]);
    const float Y((DIM>1) ? half_[1] : 1);
    const float T(top_);
    const float B(bot_);

    flute3 * flu = gym::mapBufferV3(8);
    flu[0] = {-X, Y, B};
    flu[1] = {-X, Y, T};
    flu[2] = {-X,-Y, B};
    flu[3] = {-X,-Y, T};
    flu[4] = { X, Y, B};
    flu[5] = { X, Y, T};
    flu[6] = { X,-Y, B};
    flu[7] = { X,-Y, T};
    gym::unmapBufferV3();
    // draw top and bottom faces:
    gym::enableLighting();
    gym::rebindBufferV3(2, 0);
    gym::drawTriangleStrip(0, 4);
    gym::rebindBufferV3(2, 1);
    gym::drawTriangleStrip(0, 4);
    // draw vertical edges:
    gym::disableLighting();
    gym::enableLineStipple(0x000F);
    gym::rebindBufferV3(1, 0);
    gym::drawLines(WIDTH, 0, 8);
    gym::disableLineStipple();
    gym::cleanup();
}

#else

void SpaceLid::draw2D(float) const {}
void SpaceLid::draw3D() const {}

#endif

