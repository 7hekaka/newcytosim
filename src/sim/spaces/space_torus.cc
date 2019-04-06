// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_torus.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"

SpaceTorus::SpaceTorus(const SpaceProp* p)
: Space(p)
{
    if ( DIM == 1 )
        throw InvalidParameter("torus is not usable in 1D");
}


void SpaceTorus::resize(Glossary& opt)
{
    opt.set(bCurvature, "radius");
    opt.set(bRadius,   "width");

    if ( bCurvature <= 0 )
        throw InvalidParameter("torus:width must be < radius");
    if ( bRadius <= 0 )
        throw InvalidParameter("torus:width must be < radius");
    if ( bRadius > bCurvature )
        throw InvalidParameter("torus:width must be < radius");
}


real SpaceTorus::volume() const
{
#if ( DIM == 2 )
    return 4 * M_PI * bCurvature * bRadius;
#else
    return 2 * M_PI * M_PI * bCurvature * bRadiusSqr;
#endif
}


void SpaceTorus::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-bCurvature-bRadius,-bCurvature-bRadius,-bRadius);
    sup.set( bCurvature+bRadius, bCurvature+bRadius, bRadius);
}


///project on the backbone circle in the XY plane:
Vector SpaceTorus::project0(Vector const& pos) const
{
#if ( DIM > 1 )
    real n = bCurvature / pos.normXY();
    return Vector(n * pos.XX, n * pos.YY, 0);
#else
    return Vector(0, 0, 0);
#endif
}


bool SpaceTorus::inside(Vector const& pos) const
{
    Vector prj = project0(pos);
    return ( distanceSqr(prj, pos) <= bRadiusSqr );
}


Vector SpaceTorus::project(Vector const& pos) const
{
    Vector cen = project0(pos);
    Vector ax = pos - cen;
    real n = ax.normSqr();
    n = bRadius / sqrt(n);
    return cen + n * ax;
}


//------------------------------------------------------------------------------

void SpaceTorus::write(Outputter& out) const
{
    out.put_line(" "+prop->shape+" ");
    out.writeUInt16(2);
    out.writeFloat(bCurvature);
    out.writeFloat(bRadius);
}


void SpaceTorus::setLengths(const real len[])
{
    bCurvature = len[0];
    bRadius = len[2];
    update();
}

void SpaceTorus::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    read_data(in, len);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY

#include "gle.h"
using namespace gle;

bool SpaceTorus::draw() const
{
#if ( DIM == 2 )
    
    constexpr size_t fin = 16 * gle::finesse;
    GLfloat cir[2*fin+2];

    glEnableClientState(GL_VERTEX_ARRAY);

    gle::circle(fin, cir, GLfloat(bCurvature-bRadius));
    glVertexPointer(2, GL_FLOAT, 0, cir);
    glDrawArrays(GL_LINE_STRIP, 0, fin+1);

    gle::circle(fin, cir, GLfloat(bCurvature+bRadius));
    glVertexPointer(2, GL_FLOAT, 0, cir);
    glDrawArrays(GL_LINE_STRIP, 0, fin+1);

    glDisableClientState(GL_VERTEX_ARRAY);
    return true;
    
#elif ( DIM > 2 )
    
    gleTorus(bCurvature, bRadius);
    return true;
    
#else
    return false;
#endif
}

#else

bool SpaceTorus::draw() const
{
    return false;
}


#endif
