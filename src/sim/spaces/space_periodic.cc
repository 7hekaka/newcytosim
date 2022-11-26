// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#include "dim.h"
#include "space_periodic.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"


SpacePeriodic::SpacePeriodic(SpaceProp const* p)
: Space(p)
{
    for ( int d = 0; d < 4; ++d )
        half_[d] = 0;
}


void SpacePeriodic::resize(Glossary& opt)
{
    for ( unsigned d = 0; d < DIM; ++d )
    {
        real len = half_[d];
        if ( opt.set(len, "length", d) )
            len *= 0.5;
        if ( len <= 0 )
            throw InvalidParameter("periodic:length[",d,"] must be > 0");
        half_[d] = len;
    }
    update();
}


void SpacePeriodic::update()
{
    modulo_.reset();
    for ( unsigned d = 0; d < DIM; ++d )
        modulo_.enable(d, 2*half_[d]);
}


void SpacePeriodic::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-half_[0],-half_[1],-half_[2]);
    sup.set( half_[0], half_[1], half_[2]);
}


Vector SpacePeriodic::bounce(Vector const& pos) const
{
    real X = modulo_.fold_(pos.XX, 0);
#if ( DIM == 2 )
    real Y = modulo_.fold_(pos.YY, 1);
    return Vector(X, Y);
#endif
#if ( DIM > 2 )
    real Y = modulo_.fold_(pos.YY, 1);
    real Z = modulo_.fold_(pos.ZZ, 2);
    return Vector(X, Y, Z);
#endif
    return Vector(X, 0, 0);
}


real SpacePeriodic::volume() const
{
#if ( DIM == 1 )
    return 2 * half_[0];
#elif ( DIM == 2 )
    return 4 * half_[0] * half_[1];
#else
    return 8 * half_[0] * half_[1] * half_[2];
#endif
}


bool SpacePeriodic::inside(Vector const&) const
{
    return true;
}


Vector SpacePeriodic::project(Vector const&) const
{
    throw InvalidParameter("A periodic space has no edge!");
    return Vector(0, 0, 0);
}

//------------------------------------------------------------------------------
#pragma mark - I/O

void SpacePeriodic::write(Outputter& out) const
{
    writeMarker(out, TAG);
    writeShape(out, "LLL");
    out.writeUInt16(4);
    out.writeFloat(half_[0]);
    out.writeFloat(half_[1]);
    out.writeFloat(half_[2]);
    out.writeFloat(0.f);
}


void SpacePeriodic::setLengths(const real len[])
{
    half_[0] = len[0];
    half_[1] = len[1];
    half_[2] = len[2];
    update();
}


void SpacePeriodic::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "LLL");
    setLengths(len);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------
#pragma mark -

#ifdef DISPLAY

#include "gym_flute.h"
#include "gym_draw.h"
#include "gym_view.h"
#include "gym_cap.h"

// draw periodic edges of the box
void SpacePeriodic::draw2D(float width) const
{
    const float X(half_[0]);
    const float T((DIM>1) ? half_[1] : 1);
    const float B((DIM>1) ?-half_[1] :-1);

    flute2 * flu = gym::mapBufferV2(5);
    flu[0] = { X, T };
    flu[1] = { X, B };
    flu[2] = {-X, B };
    flu[3] = {-X, T };
    flu[4] = { X, T };
    gym::unmapBufferV2();
    gym::enableLineStipple(0x000F);
    gym::drawLineStrip(width, 0, 5);
    gym::disableLineStipple();
}


// draw periodic edges of the box
void SpacePeriodic::draw3D() const
{
    const float WIDTH = 2;
    const float X(half_[0]);
    const float Y((DIM>1) ? half_[1] : 1);
    const float T((DIM>2) ? half_[2] : 0);
    const float B(-T);

    flute3 * flu = gym::mapBufferV3(10);
    flu[0] = { X, Y, B};
    flu[1] = { X, Y, T};
    flu[2] = { X,-Y, B};
    flu[3] = { X,-Y, T};
    flu[4] = {-X,-Y, B};
    flu[5] = {-X,-Y, T};
    flu[6] = {-X, Y, B};
    flu[7] = {-X, Y, T};
    flu[8] = { X, Y, B};
    flu[9] = { X, Y, T};
    gym::unmapBufferV3();
    gym::disableLighting();
    gym::enableLineStipple(0x000F);
    gym::drawLines(WIDTH, 0, 8);
    gym::rebindBufferV3(2, 0);
    gym::drawLineStrip(WIDTH, 0, 5);
    gym::rebindBufferV3(2, 1);
    gym::drawLineStrip(WIDTH, 0, 5);
    gym::disableLineStipple();
}

#else

void SpacePeriodic::draw2D(float) const {}
void SpacePeriodic::draw3D() const {}

#endif

