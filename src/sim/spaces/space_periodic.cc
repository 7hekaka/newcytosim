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


void SpacePeriodic::bounce(Vector& pos) const
{
    pos.XX = fold_real(pos.XX, modulo_.period_[0]);
#if ( DIM > 1 )
    pos.YY = fold_real(pos.YY, modulo_.period_[1]);
#endif
#if ( DIM > 2 )
    pos.ZZ = fold_real(pos.ZZ, modulo_.period_[2]);
#endif
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
    writeShape(out, "periodic");
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
    readShape(in, 8, len, "periodic");
    setLengths(len);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------
#pragma mark -

#ifdef DISPLAY
#include "opengl.h"

void SpacePeriodic::draw3D() const
{
    const GLfloat X(half_[0]);
    const GLfloat Y(( DIM > 1 ) ? half_[1] : 1);
    const GLfloat T(( DIM > 2 ) ? half_[2] : 0);
    const GLfloat B(-T);

    glLineStipple(1, 0x000F);
    glEnable(GL_LINE_STIPPLE);

#if ( DIM == 1 )
    GLfloat lin[8] = { X,-Y, X, Y,-X, Y,-X,-Y };
    glVertexPointer(2, GL_FLOAT, 0, lin);
    glDrawArrays(GL_LINES, 0, 4);
#elif ( DIM > 1 )
    // draw edges of the box
    GLfloat lin[48] = {
        +X, Y, B, X,-Y, B,-X,-Y, B,-X, Y, B,
        +X, Y, T, X,-Y, T,-X,-Y, T,-X, Y, T,
        +X, Y, B, X, Y, T, X,-Y, B, X,-Y, T,
        -X,-Y, B,-X,-Y, T,-X, Y, B,-X, Y, T
    };
    glVertexPointer(3, GL_FLOAT, 0, lin);
    glDisable(GL_LIGHTING);
    glDrawArrays(GL_LINE_LOOP, 0, 4);
    if ( DIM > 2 )
    {
        glDrawArrays(GL_LINE_LOOP, 4, 4);
        glDrawArrays(GL_LINES, 8, 8);
    }
#endif
    glDisable(GL_LINE_STIPPLE);
}

#else

void SpacePeriodic::draw3D() const {}

#endif

