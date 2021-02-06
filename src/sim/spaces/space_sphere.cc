// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_sphere.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include "random.h"
#include "meca.h"

SpaceSphere::SpaceSphere(SpaceProp const* p)
: Space(p), radius_(0)
{
}


void SpaceSphere::resize(Glossary& opt)
{
    real rad = radius_;
    
    if ( opt.set(rad, "diameter") )
        rad *= 0.5;
    else opt.set(rad, "radius");
    
    if ( rad < 0 )
        throw InvalidParameter(prop->name()+":radius must be >= 0");
    
    radius_ = rad;
}

void SpaceSphere::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-radius_,-radius_,-radius_);
    sup.set( radius_, radius_, radius_);
}


real SpaceSphere::volume() const
{
#if ( DIM == 1 )
    return 2 * radius_;
#elif ( DIM == 2 )
    return M_PI * square(radius_);
#else
    return ( 4 * M_PI / 3.0 ) * cube(radius_);
#endif
}


real SpaceSphere::surface() const
{
#if ( DIM == 1 )
    return 2;
#elif ( DIM == 2 )
    return ( 2 * M_PI ) * radius_;
#else
    return ( 4 * M_PI ) * square(radius_);
#endif
}


bool SpaceSphere::inside(Vector const& pos) const
{
    return pos.normSqr() <= square(radius_);
}

Vector SpaceSphere::project(Vector const& pos) const
{
    real n = pos.normSqr();
    
    if ( n > 0 ) {
        return pos * ( radius_ / std::sqrt(n) );
    }
    else {
        //select a random point on the surface
        return radius_ * Vector::randU();
    }
}

//------------------------------------------------------------------------------

void SpaceSphere::setInteraction(Vector const& pos, Mecapoint const& pe, Meca& meca, real stiff) const
{
    meca.addSphereClamp(pos, pe, Vector(0,0,0), radius_, stiff);
}


void SpaceSphere::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca& meca, real stiff) const
{
    if ( radius_ > rad )
        meca.addSphereClamp(pos, pe, Vector(0,0,0), radius_-rad, stiff);
    else {
        meca.addPointClamp(pe, Vector(0,0,0), stiff);
        std::cerr << "object is too big to fit in SpaceSphere\n";
    }
}


//------------------------------------------------------------------------------

void SpaceSphere::write(Outputter& out) const
{
    writeShape(out, "sphere");
    out.writeUInt16(2);
    out.writeFloat(radius_);
    out.writeFloat(0.f);
}


void SpaceSphere::setLengths(const real len[])
{
    radius_ = len[0];
}


void SpaceSphere::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "sphere");
    setLengths(len);
}


//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY

#include "gle.h"
#include "point_disp.h"

void SpaceSphere::draw2D() const
{
    //number of sections in the quarter-circle
    constexpr size_t fin = ((DIM==2) ? 32 : 8) * gle::finesse;
    
    GLfloat cir[2*fin+6];
    cir[2] = 0;
    cir[3] = 0;
    gle::circle(fin, cir+4, (GLfloat)radius_);
    
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(2, GL_FLOAT, 0, cir+4);
    glDrawArrays(GL_LINE_STRIP, 0, fin+1);
    glDisableClientState(GL_VERTEX_ARRAY);

    if ( prop->disp->visible & 2 )
    {
        prop->disp->color2.load_load();
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(2, GL_FLOAT, 0, cir+2);
        glDrawArrays(GL_TRIANGLE_FAN, 0, fin+2);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
}

void SpaceSphere::draw3D() const
{
    GLfloat R = (GLfloat)radius_;
    glPushMatrix();
    glScalef(R, R, R);
    gle::sphere8();
    gle::drawThreeBands(128);
    glPopMatrix();
}

#else

void SpaceSphere::draw2D() const {}
void SpaceSphere::draw3D() const {}

#endif
