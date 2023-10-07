// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
#include "dim.h"
#include "space_banana.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"


SpaceBanana::SpaceBanana(SpaceProp const* p)
: Space(p)
{
    if ( DIM == 1 )
        throw InvalidParameter("banana is not edible in 1D");
    bLength = 10;
    bRadius = 2;
    bCurve = 20;
}


void SpaceBanana::resize(Glossary& opt)
{
    real len = bLength, rad = bRadius, cur = bCurve;
    
    if ( opt.set(rad, "width") )
        rad /= 2;
    else opt.set(rad, "radius");
    
    opt.set(cur, "curvature");
    opt.set(len, "length");

    len = len - 2 * rad;

    if ( len <= 0 )
        throw InvalidParameter("banana:length must be >= 2 * width");
    if ( cur <= 0 )
        throw InvalidParameter("banana:curve must be >= 0");
    if ( rad < 0 )
        throw InvalidParameter("banana:radius must be > 0");
    if ( rad > bCurve )
        throw InvalidParameter("banana:radius must be <= curve");
    
    bLength = len;
    bRadius = rad;
    bCurve = cur;
    update();
}


void SpaceBanana::update()
{
    bAngle = 0.5 * bLength / bCurve;
    
    if ( bAngle > M_PI )
    {
        bAngle = M_PI;
        std::cerr << "banana:length should not exceed 2*PI*radius\n";
    }

    bRadiusSqr = bRadius * bRadius;
    bEnd[0] = bCurve * std::sin(bAngle);
    bEnd[1] = 0.5*bCurve*(1-std::cos(bAngle));
    
    bCenter[0] = 0;
    bCenter[1] = bCurve - bEnd[1];
    bCenter[2] = 0;
}


real SpaceBanana::volume() const
{
#if ( DIM > 2 )
    return (2*M_PI*bAngle*bCurve + 4*M_PI/3.*bRadius)*bRadiusSqr;
#else
    return 4*bAngle*bCurve*bRadius + M_PI*bRadiusSqr;
#endif
}


void SpaceBanana::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-bEnd[0]-bRadius,-bRadius,-bRadius);
    sup.set( bEnd[0]+bRadius, bEnd[1]+bRadius, bRadius);
}


/// project on the backbone circular arc in the XY plane:
Vector SpaceBanana::backbone(Vector const& pos) const
{
    Vector cp = pos - bCenter;
    
    real n = bCurve / cp.normXY();
    Vector prj;

    prj[0] = bCenter[0] + n * cp[0];
    prj[1] = bCenter[1] + n * cp[1];
    
    if ( prj[1] > bEnd[1] )
    {
        prj[0] = std::copysign(bEnd[0], pos[0]);
        prj[1] = bEnd[1];
    }
    
    if ( DIM > 2 )
        prj[2] = 0;
    return prj;
}


bool SpaceBanana::inside(Vector const& pos) const
{
    Vector prj = backbone(pos);
    return ( distanceSqr(pos, prj) <= bRadiusSqr );
}


Vector SpaceBanana::project(Vector const& pos) const
{
    Vector cen = backbone(pos);
    Vector dif = pos - cen;
    real n = dif.normSqr();
    return cen + (bRadius / std::sqrt(n)) * dif;
}


//------------------------------------------------------------------------------

void SpaceBanana::write(Outputter& out) const
{
    writeMarker(out, TAG);
    writeShape(out, "LRC");
    out.writeUInt16(4);
    out.writeFloat(bLength);
    out.writeFloat(bRadius);
    out.writeFloat(bCurve);
    out.writeFloat(0.f);
}


void SpaceBanana::setLengths(const real len[8])
{
    bLength = len[0];
    bRadius = len[1];
    bCurve = len[2];
    update();
}

void SpaceBanana::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "LRC");
    setLengths(len);
}

//------------------------------------------------------------------------------
#pragma mark - OpenGL display

#ifdef DISPLAY

#include "gle.h"
#include "gym_flute.h"
#include "gym_view.h"
#include "gym_draw.h"
#include "gym_cap.h"


void SpaceBanana::draw2D(float width) const
{
    float R(bRadius);
    float C(bCurve);
    float cX(bCenter[0]);
    float cY(bCenter[1]);
    float eX(bEnd[0]);
    float eY(bEnd[1]);
    float A(-bAngle + M_PI_2);
    float B( bAngle - M_PI_2);

    //number of sections in the quarter-circle
    constexpr size_t fin = 8 * gle::finesse;
    float* arc = (float*)gym::mapBufferV2(4*fin+4);
    // lower swing
    gle::compute_arc(fin, arc      , C+R, A-M_PI, M_PI+B-A, cX, cY);
    // right cap
    gle::compute_arc(fin, arc+2*fin, R, B, M_PI, eX, eY);
    // upper swing
    gle::compute_arc(fin, arc+4*fin, C-R, B, A-B-M_PI, cX, cY);
    // left cap
    gle::compute_arc(fin, arc+6*fin, R, A, M_PI, -eX, eY);
    arc[0+8*fin] = arc[0];
    arc[1+8*fin] = arc[1];
    gym::unmapBufferV2();
    gym::drawLineStrip(width, 0, 4*fin+1);
}

void SpaceBanana::draw3D() const
{
    float U(bCurve);
    float R(bRadius);
    
    double C = std::cos(bAngle), S = std::sin(bAngle);
    
    gym::enableClipPlane(4);
    gym::enableClipPlane(5);
    
    //center part:
    gym::transScale(bCenter[0], bCenter[1], 0, R);
    gym::setClipPlane(4, -C, -S, 0, 0);
    gym::setClipPlane(5,  C, -S, 0, 0);
    gle::torusZ(U, 1);

    gym::disableClipPlane(5);

    //right cap:
    gym::transScale(bEnd[0], bEnd[1], 0, R);
    gym::setClipPlane(4, C, S, 0, 0);
    gle::sphere();

    //left cap:
    gym::transScale(-bEnd[0], bEnd[1], 0, R);
    gym::setClipPlane(4, -C, S, 0, 0);
    gle::sphere();

    gym::disableClipPlane(4);
}

#else

void SpaceBanana::draw2D(float) const {}
void SpaceBanana::draw3D() const {}

#endif
