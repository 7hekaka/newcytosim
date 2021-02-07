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
    halflength_[0] = 0;
    halflength_[1] = 0;
    bot_ = 0;
    top_ = 0;
}


void SpaceStrip::resize(Glossary& opt)
{
    for ( unsigned d = 0; d < DIM-1; ++d )
    {
        real len = halflength_[d];
        if ( opt.set(len, "length", d) )
            len *= 0.5;
        if ( len < 0 )
            throw InvalidParameter("strip:length[] must be >= 0");
        halflength_[d] = len;
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
        modulo_.enable(d, 2*halflength_[d]);
}


void SpaceStrip::boundaries(Vector& inf, Vector& sup) const
{
#if ( DIM >= 3 )
    inf.set(-halflength_[0],-halflength_[1], bot_);
    sup.set( halflength_[0], halflength_[1], top_);
#else
    inf.set(-halflength_[0], bot_, 0);
    sup.set( halflength_[0], top_, 0);
#endif
}


void SpaceStrip::bounce(Vector& pos) const
{
    if ( !SpaceStrip::inside(pos) )
        bounceOnEdges(pos);
    
    // periodic in all except the last dimension:
#if ( DIM > 1 )
    pos.XX = fold_real(pos.XX, modulo_.period_[0]);
#endif
#if ( DIM > 2 )
    pos.YY = fold_real(pos.YY, modulo_.period_[1]);
#endif
}


//------------------------------------------------------------------------------
#pragma mark -


real SpaceStrip::volume() const
{
#if ( DIM == 1 )
    return ( top_ - bot_ );
#elif ( DIM == 2 )
    return 2 * halflength_[0] * ( top_ - bot_ );
#else
    return 4 * halflength_[0] * halflength_[1] * ( top_ - bot_ );
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
    real X = sign_select(2 * pos.XX - bot_ - top_, bot_, top_);
    return Vector(X);
#elif ( DIM == 2 )
    real Y = sign_select(2 * pos.YY - bot_ - top_, bot_, top_);
    return Vector(pos.XX, Y);
#else
    real Z = sign_select(2 * pos.ZZ - bot_ - top_, bot_, top_);
    return Vector(pos.XX, pos.YY, Z);
#endif
}


//------------------------------------------------------------------------------
#pragma mark - setInteraction


void SpaceStrip::setInteraction(Vector const& pos, Mecapoint const& pe, Meca& meca, real stiff) const
{
#if ( DIM == 2 )
    real Y = sign_select(2 * pos.YY - bot_ - top_, bot_, top_);
    meca.addPlaneClampY(pe, Y, stiff);
#elif ( DIM > 2 )
    real Z = sign_select(2 * pos.ZZ - bot_ - top_, bot_, top_);
    meca.addPlaneClampZ(pe, Z, stiff);
#endif
}


void SpaceStrip::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca& meca, real stiff) const
{
#if ( DIM == 2 )
    real Y = sign_select(2 * pos.YY - bot_ - top_, bot_+rad, top_-rad);
    meca.addPlaneClampY(pe, Y, stiff);
#elif ( DIM > 2 )
    real Z = sign_select(2 * pos.ZZ - bot_ - top_, bot_+rad, top_-rad);
    meca.addPlaneClampZ(pe, Z, stiff);
#endif
}

//------------------------------------------------------------------------------

void SpaceStrip::write(Outputter& out) const
{
    writeShape(out, "strip");
    out.writeUInt16(4);
    out.writeFloat(halflength_[0]);
    out.writeFloat(halflength_[1]);
    out.writeFloat(bot_);
    out.writeFloat(top_);
}


void SpaceStrip::setLengths(const real len[])
{
    halflength_[0] = len[0];
    halflength_[1] = len[1];
    bot_ = len[2];
    top_ = len[3];
#ifdef BACKWARD_COMPATIBILITY
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
    const GLfloat X(halflength_[0]);
    const GLfloat T(top_);
    const GLfloat B(bot_);
    
    glBegin(GL_LINES);
    glVertex3f(-X, T, 0);
    glVertex3f( X, T, 0);
    glVertex3f( X, B, 0);
    glVertex3f(-X, B, 0);
    glEnd();
    
    // draw periodic boundaries:
    glLineStipple(1, 0x000F);
    glEnable(GL_LINE_STIPPLE);
    glBegin(GL_LINES);
    glVertex3f( X, T, 0);
    glVertex3f( X, B, 0);
    glVertex3f(-X, T, 0);
    glVertex3f(-X, B, 0);
    glEnd();
    glDisable(GL_LINE_STIPPLE);
}

void SpaceStrip::draw3D() const
{
    const GLfloat X(halflength_[0]);
    const GLfloat T(top_);
    const GLfloat B(bot_);

    const GLfloat Y(halflength_[1]);
    // draw faces:
    glBegin(GL_TRIANGLE_STRIP);
    glNormal3f(0, 0, 1);
    glVertex3f(-X,  Y, B);
    glVertex3f( X,  Y, B);
    glVertex3f(-X, -Y, B);
    glVertex3f( X, -Y, B);
    glEnd();
    glBegin(GL_TRIANGLE_STRIP);
    glNormal3f(0, 0, -1);
    glVertex3f(-X,  Y, T);
    glVertex3f(-X, -Y, T);
    glVertex3f( X,  Y, T);
    glVertex3f( X, -Y, T);
    glEnd();
    // draw outline:
    glBegin(GL_LINE_STRIP);
    glVertex3f(-X,  Y, B);
    glVertex3f( X,  Y, B);
    glVertex3f( X, -Y, B);
    glVertex3f(-X, -Y, B);
    glVertex3f(-X,  Y, B);
    glEnd();
    glBegin(GL_LINE_STRIP);
    glVertex3f(-X,  Y, T);
    glVertex3f(-X, -Y, T);
    glVertex3f( X, -Y, T);
    glVertex3f( X,  Y, T);
    glVertex3f(-X,  Y, T);
    glEnd();
    
    // draw periodic boundaries:
    glLineStipple(1, 0x000F);
    glEnable(GL_LINE_STIPPLE);
    glBegin(GL_LINES);
    glVertex3f( X,  Y, T);
    glVertex3f( X,  Y, B);
    glVertex3f( X, -Y, T);
    glVertex3f( X, -Y, B);
    glVertex3f(-X,  Y, T);
    glVertex3f(-X,  Y, B);
    glVertex3f(-X, -Y, T);
    glVertex3f(-X, -Y, B);
    glEnd();
    glDisable(GL_LINE_STIPPLE);
}

#else

void SpaceStrip::draw2D() const {}
void SpaceStrip::draw3D() const {}

#endif

