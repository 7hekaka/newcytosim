// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_banana.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"


SpaceBanana::SpaceBanana(const SpaceProp* p)
: Space(p)
{
    if ( DIM == 1 )
        throw InvalidParameter("banana is not edible in 1D");
}


void SpaceBanana::resize(Glossary& opt)
{
    real total = 0;
    opt.set(total,   "length");
    opt.set(bCurvature,   "curvature");
    opt.set(bRadius, "radius");

    bLength = total - 2 * bCurvature;
    bRadiusSqr = bRadius * bRadius;

    if ( bLength <= 0 )
        throw InvalidParameter("banana:length must be specified and greater than 2 * width");
    if ( bCurvature <= 0 )
        throw InvalidParameter("banana:radius must be >= 0");
    if ( bRadius < 0 )
        throw InvalidParameter("banana:width must be > 0");
    if ( bRadius > bCurvature )
        throw InvalidParameter("banana:width must be smaller than radius");
   
    bAngle = 0.5 * bLength / bCurvature;
    
    if ( bAngle > M_PI )
    {
        bAngle = M_PI;
        std::cerr << "banana:length should not exceed 2*PI*radius\n";
    }
    update();
}


void SpaceBanana::update()
{
    bRadiusSqr = bRadius * bRadius;
    bEnd[0] = bCurvature * sin(bAngle);
    bEnd[1] = 0.5*bCurvature*(1-cos(bAngle));
    
    bCenter[0] = 0;
    bCenter[1] = bCurvature - bEnd[1];
    bCenter[2] = 0;
}


real SpaceBanana::volume() const
{
#if ( DIM > 2 )
    return 2*bAngle*M_PI*bCurvature*bRadiusSqr + 4/3.*M_PI*bRadiusSqr*bRadius;
#else
    return 4*bAngle*bCurvature*bRadius + M_PI*bRadiusSqr;
#endif
}


void SpaceBanana::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-bEnd[0]-bRadius,-bRadius,-bRadius);
    sup.set( bEnd[0]+bRadius, bEnd[1]+bRadius, bRadius);
}


/// project on the backbone circular arc in the XY plane:
Vector SpaceBanana::project0(Vector const& pos) const
{
    Vector cp = pos - bCenter;
    
    real n = bCurvature / cp.normXY();
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
    Vector prj = project0(pos);
    return ( distanceSqr(pos, prj) <= bRadiusSqr );
}


Vector SpaceBanana::project(Vector const& pos) const
{
    Vector cen = project0(pos);
    Vector dif = pos - cen;
    real n = dif.normSqr();
    return cen + (bRadius / sqrt(n)) * dif;
}


//------------------------------------------------------------------------------

void SpaceBanana::write(Outputter& out) const
{
    out.put_line(" "+prop->shape+" ");
    out.writeUInt16(3);
    out.writeFloat(bLength);
    out.writeFloat(bCurvature);
    out.writeFloat(bRadius);
}


void SpaceBanana::setLengths(const real len[])
{
    bLength    = len[0];
    bCurvature = len[1];
    bRadius    = len[2];
    update();
}

void SpaceBanana::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    read_data(in, len);
    setLengths(len);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY

#include "gle.h"
using namespace gle;

bool SpaceBanana::draw() const
{
#if ( DIM == 2 )
    
    //number of sections in the quarter-circle
    constexpr size_t fin = 8 * gle::finesse;
    GLfloat c[fin+1], s[fin+1];
    
    GLfloat A = -bAngle + M_PI_2;
    GLfloat B =  bAngle - M_PI_2;
    
    glBegin(GL_LINE_LOOP);
    // lower swing
    gle::arc(fin, c, s, bCurvature+bRadius, A-M_PI, B, bCenter[0], bCenter[1]);
    for ( size_t i = 0; i < fin; ++i )
        gleVertex(c[i], s[i], 0);
    
    // right cap
    gle::arc(fin, c, s, bRadius, B, B+M_PI, bEnd[0], bEnd[1]);
    for ( size_t i = 0; i < fin; ++i )
        gleVertex(c[i], s[i], 0);
    
    // upper swing
    gle::arc(fin, c, s, bCurvature-bRadius, B, A-M_PI, bCenter[0], bCenter[1]);
    for ( size_t i = 0; i < fin; ++i )
        gleVertex(c[i], s[i], 0);
        
    // left cap
    gle::arc(fin, c, s, bRadius, A, A+M_PI, -bEnd[0], bEnd[1]);
    for ( size_t i = 0; i < fin; ++i )
        gleVertex(c[i], s[i], 0);

    glEnd();
    
#elif ( DIM > 2 )
    
    glMatrixMode(GL_MODELVIEW);

    GLdouble plane1[] = { -cos(bAngle), -sin(bAngle), 0, 0 };
    GLdouble plane2[] = {  cos(bAngle), -sin(bAngle), 0, 0 };
    GLdouble plane1i[] = {  cos(bAngle), sin(bAngle), 0, 0 };
    GLdouble plane2i[] = { -cos(bAngle), sin(bAngle), 0, 0 };
    
    const GLenum glp1 = GL_CLIP_PLANE4;
    const GLenum glp2 = GL_CLIP_PLANE5;
    
    glEnable(glp1);
    glEnable(glp2);
    
    //center part:
    glPushMatrix();
    glTranslated(bCenter[0], bCenter[1], 0);
    glClipPlane(glp1, plane1);
    glClipPlane(glp2, plane2);
    gleTorus(bCurvature, bRadius);
    glPopMatrix();

    glDisable(glp2);

    //right cap:
    glPushMatrix();
    glTranslated(bEnd[0], bEnd[1], 0);
    glClipPlane(glp1, plane1i);
    gleScale(bRadius);
    gleSphere8B();
    glPopMatrix();

    //left cap:
    glPushMatrix();
    glTranslated(-bEnd[0], bEnd[1], 0);
    glClipPlane(glp1, plane2i);
    gleScale(bRadius);
    gleSphere8B();
    glPopMatrix();
    
    glDisable(glp1);
    
#endif
    return true;
}

#else

bool SpaceBanana::draw() const
{
    return false;
}


#endif
