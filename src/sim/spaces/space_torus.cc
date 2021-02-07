// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
#include "dim.h"
#include "space_torus.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"

SpaceTorus::SpaceTorus(SpaceProp const* p)
: Space(p)
{
    if ( DIM == 1 )
        throw InvalidParameter("torus is not usable in 1D");
    bRadius = 2;
    bCurve = INFINITY;
}


void SpaceTorus::resize(Glossary& opt)
{
    real rad = bRadius, cur = bCurve;
    
    if ( opt.set(rad, "width") )
        rad /= 2;
    else opt.set(rad, "radius");
    
    opt.set(cur, "curvature");

    if ( cur <= 0 )
        throw InvalidParameter("torus:curve must be >= 0");
    if ( rad < 0 )
        throw InvalidParameter("torus:radius must be > 0");
    if ( rad > bCurve )
        throw InvalidParameter("torus:radius must be <= curve");
    
    bRadius = rad;
    bCurve = cur;
    update();
}


real SpaceTorus::volume() const
{
#if ( DIM == 2 )
    return 4 * M_PI * bCurve * bRadius;
#else
    return 2 * M_PI * M_PI * bCurve * bRadiusSqr;
#endif
}


void SpaceTorus::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-bCurve-bRadius,-bCurve-bRadius,-bRadius);
    sup.set( bCurve+bRadius, bCurve+bRadius, bRadius);
}


///project on the backbone circle in the XY plane:
Vector SpaceTorus::backbone(Vector const& pos) const
{
#if ( DIM > 1 )
    real n = bCurve / pos.normXY();
    return Vector(n * pos.XX, n * pos.YY, 0);
#else
    return Vector(0, 0, 0);
#endif
}


bool SpaceTorus::inside(Vector const& pos) const
{
    Vector prj = backbone(pos);
    return ( distanceSqr(prj, pos) <= bRadiusSqr );
}


Vector SpaceTorus::project(Vector const& pos) const
{
    Vector cen = backbone(pos);
    Vector ax = pos - cen;
    real n = ax.normSqr();
    n = bRadius / std::sqrt(n);
    return cen + n * ax;
}


//------------------------------------------------------------------------------

void SpaceTorus::write(Outputter& out) const
{
    writeShape(out, "torus");
    out.writeUInt16(2);
    out.writeFloat(bCurve);
    out.writeFloat(bRadius);
}


void SpaceTorus::setLengths(const real len[])
{
    bCurve = len[0];
    bRadius = len[2];
    update();
}

void SpaceTorus::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "torus");
    setLengths(len);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY

#include "gle.h"
using namespace gle;

void SpaceTorus::draw2D() const
{
    constexpr size_t fin = 16 * gle::finesse;
    GLfloat cir[2*fin+2];

    gle::circle(fin, cir, GLfloat(bCurve-bRadius));
    glVertexPointer(2, GL_FLOAT, 0, cir);
    glDrawArrays(GL_LINE_STRIP, 0, fin+1);

    gle::circle(fin, cir, GLfloat(bCurve+bRadius));
    glVertexPointer(2, GL_FLOAT, 0, cir);
    glDrawArrays(GL_LINE_STRIP, 0, fin+1);
}

void SpaceTorus::draw3D() const
{
    gle::torusZ(bCurve, bRadius);
}

#else

void SpaceTorus::draw2D() const {}
void SpaceTorus::draw3D() const {}

#endif
