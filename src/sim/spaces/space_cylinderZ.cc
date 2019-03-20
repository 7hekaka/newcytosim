// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_cylinderZ.h"
#include "mecapoint.h"
#include "exceptions.h"
#include "meca.h"


SpaceCylinderZ::SpaceCylinderZ(const SpaceProp* p)
: Space(p), radius(mLength[0]), bot(mLength[1]), top(mLength[2])
{
    if ( DIM < 3 )
        throw InvalidParameter("cylinderZ is only valid in 3D: use sphere instead");
}


void SpaceCylinderZ::resize()
{
    Space::checkLengths(1, false);
    if ( top < bot )
        throw InvalidParameter("bottom must be lower than top ( size[1] <= size[2] )");
}


void SpaceCylinderZ::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-radius,-radius, bot);
    sup.set( radius, radius, top);
}



real  SpaceCylinderZ::volume() const
{
    return 2 * M_PI * ( top - bot ) * radius * radius;
}


bool  SpaceCylinderZ::inside(Vector const& w) const
{
#if ( DIM > 2 )
    if ( w.ZZ < bot ) return false;
    if ( w.ZZ > top ) return false;
#endif
#if ( DIM > 1 )
    return ( w.XX*w.XX + w.YY*w.YY <= radius * radius );
#else
    return ( fabs(w.XX) <= radius );
#endif
}

bool  SpaceCylinderZ::allInside(Vector const& w, const real rad ) const
{
    assert_true( rad >= 0 );
#if ( DIM > 2 )
    if ( w.ZZ - rad < bot ) return false;
    if ( w.ZZ + rad > top ) return false;
#endif
#if ( DIM > 1 )
    return ( w.XX*w.XX + w.YY*w.YY <= square(radius-rad) );
#else
    return ( fabs(w.XX) <= radius-rad );
#endif
}

Vector SpaceCylinderZ::randomPlace() const
{
    Vector2 sec = Vector2::randB(radius);
    return Vector(sec.XX, sec.YY, bot+RNG.preal()*(top-bot));
}

//------------------------------------------------------------------------------
Vector SpaceCylinderZ::project(Vector const& w) const
{
    Vector p = w;
#if ( DIM >= 3 )
    bool inZ = true;
    
    if ( w.ZZ > top )
    {
        p.ZZ = top;
        inZ = false;
    }
    else if ( w.ZZ < bot )
    {
        p.ZZ = bot;
        inZ = false;
    }

    real n = w.normXY();
    
    if ( n > radius )
    {
        n = radius / n;
        p.XX = n * w.XX;
        p.YY = n * w.YY;
    }
    else
    {
        if ( inZ )
        {
            if ( top - w.ZZ < radius - n )
                p.ZZ = top;
            else if ( w.ZZ - bot < radius - n )
                p.ZZ = bot;
            else
            {
                n = radius / n;
                p.XX = n * w.XX;
                p.YY = n * w.YY;
            }
        }
    }
#endif
    return p;
}

//------------------------------------------------------------------------------

/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCylinderZ::setInteraction(Vector const& pos, Mecapoint const& pe, Meca & meca, real stiff,
                                    const real rad, const real bot, const real top)
{
#if ( DIM >= 3 )
    bool cap = false;
    bool cyl = false;
    real zzz;

    // inside cylinder radius
    if ( 2 * pos.ZZ - bot > top )
    {
        zzz = top;
        cap = ( pos.ZZ > top );
    }
    else
    {
        zzz = bot;
        cap = ( pos.ZZ < bot );
    }
    
    real dis = pos.XX*pos.XX + pos.YY*pos.YY;
    
    if ( rad*rad < dis )
    {
        // outside cylinder in XY plane
        cyl = true;
    }
    else if ( ! cap )
    {
        // inside cylinder in XY plane and also inside in Z:
        if ( fabs(pos.ZZ-zzz) > rad - sqrt(dis) )
        //if ( dis > rad*rad + square(pos.ZZ-p) - 2 * rad * fabs(pos.ZZ-p) )
            cyl = true;
        else
            cap = true;
    }
    
    if ( cap )
    {
        const Matrix::index_t inx = 2 + DIM * pe.matIndex();
        meca.mC(inx, inx) -= stiff;
        meca.base(inx)    += stiff * zzz;
    }
    
    if ( cyl )
        meca.addCylinderClampZ(pe, rad, stiff);
#endif
}


/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCylinderZ::setInteraction(Vector const& pos, Mecapoint const& pe, Meca & meca, real stiff) const
{
    setInteraction(pos, pe, meca, stiff, radius, bot, top);
}

/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCylinderZ::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca & meca, real stiff) const
{
    real R = std::max((real)0, radius - rad);
    real T = top - rad;
    real B = bot + rad;
    
    if ( B > T )
    {
        B = 0.5 * ( top + bot );
        T = B;
    }
    
    setInteraction(pos, pe, meca, stiff, R, B, T);
}


//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"

bool SpaceCylinderZ::draw() const
{
#if ( DIM > 2 )
    
    GLfloat T = top;
    GLfloat B = bot;
    GLfloat R = radius;
    
    const size_t fin = 512;
    GLfloat c[fin+1], s[fin+1];
    gle::circle(fin, c, s, 1);
    
    glBegin(GL_TRIANGLE_STRIP);
    //display strips along the side of the volume:
    for ( size_t n = 0; n <= fin; ++n )
    {
        glNormal3f(c[n], s[n], 0);
        glVertex3f(R*c[n], R*s[n], T);
        glVertex3f(R*c[n], R*s[n], B);
    }
    glEnd();
    
    // draw top cap:
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, +1);
    glVertex3f(0, 0,  T);
    for ( size_t n = 0; n <= fin; ++n )
        glVertex3f(R*c[n], R*s[n], T);
    glEnd();
    
    // draw bottom cap:
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, -1);
    glVertex3f(0, 0,  B);
    for ( size_t n = 0; n <= fin; ++n )
        glVertex3f(-R*c[n], R*s[n], B);
    glEnd();
    
#endif
    return true;
}

#else

bool SpaceCylinderZ::draw() const
{
    return false;
}

#endif
