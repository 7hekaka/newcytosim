// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_cylinderZ.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "mecapoint.h"
#include "glossary.h"
#include "meca.h"

/// keyword to allow smooth edges on the Cylinder
#define HAS_SMOOTH_EDGES 1


SpaceCylinderZ::SpaceCylinderZ(SpaceProp const* p)
: Space(p), radius_(0), bot_(0), top_(0), edge_(0)
{
    if ( DIM < 3 )
        throw InvalidParameter("cylinderZ is only valid in 3D");
    update();
}


void SpaceCylinderZ::resize(Glossary& opt)
{
    real rad = radius_, bot = bot_, top = top_, edg = edge_;

    if ( opt.set(rad, "diameter") )
        rad *= 0.5;
    else opt.set(rad, "radius");
    
    if ( rad < 0 )
        throw InvalidParameter("cylinderZ:radius must be >= 0");

    if ( opt.set(top, "length") )
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
        throw InvalidParameter("cylinderZ:bottom must be <= top");
    
#if HAS_SMOOTH_EDGES
    if ( opt.set(edg, "edge") )
    {
        if ( edg < 0 )
            throw InvalidParameter("cylinderZ:edge must be >= 0");
        if ( rad < edg )
            throw InvalidParameter("cylinderZ:edge must be <= radius");
        if ( top - bot < 2 * edg )
            throw InvalidParameter("cylinderZ:edge must be <= length / 2");
    }
#endif
    
    radius_ = rad;
    bot_ = bot;
    top_ = top;
    edge_ = edg;
    update();
}


void SpaceCylinderZ::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-radius_,-radius_, bot_);
    sup.set( radius_, radius_, top_);
}

/**
 Pappus's (2nd) Centroid Theorem: the volume of a planar area of revolution is
 the product of the area A and the length of the path traced by its centroid R,
 i.e., 2 M_PI R.
 
 https://en.wikipedia.org/wiki/Pappus%27s_centroid_theorem
 */
real SpaceCylinderZ::volume() const
{
#if HAS_SMOOTH_EDGES
    const real R0 = radius_ - edge_;
    return M_PI * (( top_ - bot_ - 2 * edge_ ) * radius_ * radius_  //central cylindrical part
                   + 2 * edge_ * R0 * R0                            //top and bottom central cylinders
                   + ( R0 * M_PI + 4/3.0 * edge_ ) * edgeSqr_ );
#else
    return M_PI * ( top_ - bot_ ) * radius_ * radius_;
#endif
}


bool SpaceCylinderZ::inside(Vector const& W) const
{
#if ( DIM > 2 )
    const real RT = W.XX * W.XX + W.YY * W.YY;
# if HAS_SMOOTH_EDGES
    const real R = max_real(0, sqrt(RT)-radius_+edge_);
    const real Z = max_real(0, std::max(bot_+edge_-W.ZZ, W.ZZ-top_+edge_));
    return ( R*R + Z*Z <= edgeSqr_ );
# else
    return (( bot_ <= W.ZZ ) & ( W.ZZ <= top_ ) & ( RT <= radius_ * radius_ ));
# endif
#else
    ABORT_NOW("cylinderZ is only valid in 3D");
    return true;
#endif
}


bool SpaceCylinderZ::allInside(Vector const& W, const real rad) const
{
    assert_true( rad >= 0 );
#if ( DIM > 2 )
    const real RT = W.XX * W.XX + W.YY * W.YY;
# if HAS_SMOOTH_EDGES
    const real E = edge_ + rad;
    const real R = max_real(0, sqrt(RT)-E);
    const real Z = max_real(0, max_real(bot_+E-W.ZZ, W.ZZ-top_+E));
    return ( R*R + Z*Z <= edgeSqr_ );
# else
    return ((bot_ + rad <= W.ZZ ) & ( W.ZZ + rad <= top_ )
            & ( RT <= square(radius_-rad) ));
# endif
#else
    ABORT_NOW("cylinderZ is only valid in 3D");
    return true;
#endif
}


Vector SpaceCylinderZ::randomPlace() const
{
#if HAS_SMOOTH_EDGES
    Vector W;
    do {
        const Vector2 V = Vector2::randB(radius_);
        W.set(V.XX, V.YY, bot_+RNG.preal()*(top_-bot_));
    } while ( !inside(W) );
    return W;
#else
    const Vector2 V = Vector2::randB(radius_);
    return Vector(V.XX, V.YY, bot_+RNG.preal()*(top_-bot_));
#endif
}

//------------------------------------------------------------------------------
Vector SpaceCylinderZ::project(Vector const& W) const
{
    Vector P(W);
#if ( DIM > 2 )
# if HAS_SMOOTH_EDGES
    const real T = top_ - edge_;
    const real B = bot_ + edge_;
    const real R = radius_ - edge_;
#else
    real const& T = top_;
    real const& B = bot_;
    real const& R = radius_;
# endif
    bool in = true;

    if ( T < W.ZZ )
    {
        P.ZZ = T;
        in = false;
    }
    else if ( W.ZZ < B )
    {
        P.ZZ = B;
        in = false;
    }
    
    real n = W.normXY();
    
    if ( n > R )
    {
        n = R / n;
        P.XX = n * W.XX;
        P.YY = n * W.YY;
    }
    else if ( in )
    {
        real ZT = top_ - W.ZZ;
        real ZB = W.ZZ - bot_;
        // check which cap is closer:
        if ( radius_ - n < ZT  &&  radius_ - n < ZB )
        {
            n = radius_ / n;
            P.XX = n * W.XX;
            P.YY = n * W.YY;
        }
        else if ( ZT < ZB )
            P.ZZ = top_;
        else
            P.ZZ = bot_;
        return P;
    }
# if HAS_SMOOTH_EDGES
    //normalize to radius(), and add to p to get the projection
    real dis = edge_ / norm(W-P);
    return dis * ( W - P ) + P;
# endif
#endif
    return P;
}

//------------------------------------------------------------------------------

/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCylinderZ::setInteraction(Vector const& pos, Mecapoint const& pe, Meca& meca, 
                                    real stiff, real R, real B, real T)
{
#if ( DIM >= 3 )
    bool cap = false;
    bool cyl = false;
    real Z;

    // inside cylinder radius_
    if ( 2 * pos.ZZ - B > T )
    {
        Z = T;
        cap = ( pos.ZZ > T );
    }
    else
    {
        Z = B;
        cap = ( pos.ZZ < B );
    }
    
    real dis = pos.XX*pos.XX + pos.YY*pos.YY;
    
    if ( R*R < dis )
    {
        // outside cylinder in XY plane
        cyl = true;
    }
    else if ( ! cap )
    {
        // inside cylinder in XY plane and also inside in Z:
        //if ( abs_real(pos.ZZ-Z) > R - sqrt(dis) )
        if ( dis > square( R - abs_real(pos.ZZ-Z) ) )
            cyl = true;
        else
            cap = true;
    }
    
    if ( cap )
        meca.addPlaneClampZ(pe, Z, stiff);
    
    if ( cyl )
        meca.addCylinderClampZ(pe, R, stiff);
#endif
}


void SpaceCylinderZ::setInteraction(Vector const& pos, Mecapoint const& pe, Meca& meca,
                                    real stiff, real R, real B, real T, real E)
{
#if ( DIM >= 3 )
    real n = pos.normXY();
    bool in = true;
    real Z;
    
    // inside cylinder radius_
    if ( 2 * pos.ZZ - B > T )
    {
        Z = T-E;
        if ( pos.ZZ > T-E )
        {
            in = false;
            if ( n < R-E ) {
                meca.addPlaneClampZ(pe, T, stiff);
                return;
            }
        }
    }
    else
    {
        Z = B+E;
        if ( pos.ZZ < B+E )
        {
            in = false;
            if ( n < R-E ) {
                meca.addPlaneClampZ(pe, B, stiff);
                return;
            }
        }
    }
    
    if ( in )
    {
        if ( n > R-E  ||  R-E-n < abs_real(pos.ZZ-Z) )
            meca.addCylinderClampZ(pe, R, stiff);
        else
            meca.addPlaneClampZ(pe, Z, stiff);
        return;
    }
    
    /*
     we add interaction to the tangent plane, which is okay
     if the edge radius is not too small
     */
    real S = ( R - E ) / n;
    Vector3 rim(S*pos.XX, S*pos.YY, Z);  // on the reduced cylinder rim
    Vector dir = normalize(pos-rim);
    meca.addPlaneClamp(pe, rim+E*dir, dir, stiff);

#endif
}

/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCylinderZ::setInteraction(Vector const& pos, Mecapoint const& pe, Meca& meca, real stiff) const
{
#if HAS_SMOOTH_EDGES
    setInteraction(pos, pe, meca, stiff, radius_, bot_, top_, edge_);
#else
    setInteraction(pos, pe, meca, stiff, radius_, bot_, top_);
#endif
}

/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCylinderZ::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca& meca, real stiff) const
{
    real R = max_real(0, radius_ - rad);
    real T = top_ - rad;
    real B = bot_ + rad;
    
    if ( B > T )
    {
        B = 0.5 * ( top_ + bot_ );
        T = B;
    }
    
#if HAS_SMOOTH_EDGES
    setInteraction(pos, pe, meca, stiff, R, B, T, edge_);
#else
    setInteraction(pos, pe, meca, stiff, R, B, T);
#endif
}


//------------------------------------------------------------------------------

void SpaceCylinderZ::write(Outputter& out) const
{
    out.put_characters("cylinderZ", 16);
    out.writeUInt16(4);
    out.writeFloat(radius_);
    out.writeFloat(bot_);
    out.writeFloat(top_);
    out.writeFloat(edge_);
}


void SpaceCylinderZ::setLengths(const real len[])
{
    radius_ = len[0];
    bot_    = len[1];
    top_    = len[2];
    edge_   = len[3];
}


void SpaceCylinderZ::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    read_data(in, len, "cylinderZ");
    setLengths(len);
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
    GLfloat T = top_;
    GLfloat B = bot_;
    GLfloat R = radius_;
    GLfloat E = edge_;
    GLfloat RE = radius_ - edge_;
    GLfloat TE = top_ - edge_;
    GLfloat BE = bot_ + edge_;
    
    const size_t fin = 256;
    GLfloat c[fin+1], s[fin+1];
    gle::circle(fin, c, s, 1);

    size_t pi_half = fin/4;
    size_t pi_once = fin/2;

    for (size_t i = 0; i < fin; i++)
    {
        GLfloat CU = c[i],   SU = s[i];
        GLfloat CL = c[i+1], SL = s[i+1];

        glBegin(GL_TRIANGLE_STRIP);
        glNormal3f(0, 0, +1);
        glVertex3f(0, 0,  T);
        if ( edge_ > 0 )
        {
            //draw top arc
            for ( size_t j = 0; j <= pi_half; j++ )
            {
                glNormal3f(CU*s[j],        SU*s[j],             c[j]);
                glVertex3f(CU*(RE+E*s[j]), SU*(RE+E*s[j]), TE+E*c[j]);
                glNormal3f(CL*s[j],        SL*s[j],             c[j]);
                glVertex3f(CL*(RE+E*s[j]), SL*(RE+E*s[j]), TE+E*c[j]);
            }
            
            //draw bottom arc
            for ( size_t j = pi_half; j<=pi_once; j++ )
            {
                glNormal3f(CU*s[j],        SU*s[j],             c[j]);
                glVertex3f(CU*(RE+E*s[j]), SU*(RE+E*s[j]), BE+E*c[j]);
                glNormal3f(CL*s[j],        SL*s[j],             c[j]);
                glVertex3f(CL*(RE+E*s[j]), SL*(RE+E*s[j]), BE+E*c[j]);
            }
        }
        else
        {
            glNormal3f(CU,   SU,   0);
            glVertex3f(CU*R, SU*R, T);
            glNormal3f(CL,   SL,   0);
            glVertex3f(CL*R, SL*R, T);
            glNormal3f(CU,   SU,   0);
            glVertex3f(CU*R, SU*R, B);
            glNormal3f(CL,   SL,   0);
            glVertex3f(CL*R, SL*R, B);
        }
        glNormal3f(0, 0, -1);
        glVertex3f(0, 0,  B);
        glEnd();
    }
#endif
    return true;
}

#else

bool SpaceCylinderZ::draw() const
{
    return false;
}

#endif
