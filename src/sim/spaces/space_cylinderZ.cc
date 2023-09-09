// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University
#include "space_cylinderZ.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "mecapoint.h"
#include "glossary.h"
#include "meca.h"


/// keyword to allow smooth edges on the Cylinder
#define SMOOTH_CYLINDER 1


/// surface of a cylinder of given radius, length
static real surface_cylinder(real rad, real len)
{
    // surface elements divided by 2 * M_PI:
    real S0 = rad * len;
    real S1 = square(rad);
    return ( 2 * M_PI ) * ( S0 + S1 );
}

/// surface of a cylinder of given radius, length, with rounded edges
static real surface_smooth_cylinder(real rad, real len, real edg)
{
    real RE = rad - edg;
    real LE = len - 2 * edg;
    real GC = RE + edg * ( 2.0 / M_PI );
    // surface elements divided by 2 * M_PI:
    real S0 = rad * LE;        // cylindrical surface
    real S1 = square(RE);      // top/bottom discs
    real S2 = M_PI * GC * edg; // rounded edges
    return ( 2 * M_PI ) * ( S0 + S1 + S2 );
}


SpaceCylinderZ::SpaceCylinderZ(SpaceProp const* p)
: Space(p), radius_(0), bot_(0), top_(0), edge_(0)
{
    if ( DIM < 3 )
        throw InvalidParameter("cylinderZ is only valid in 3D");
    radius_ = 0;
    bot_ = 0;
    top_ = 0;
    edge_ = 0;
    update();
}


void SpaceCylinderZ::resize(Glossary& opt)
{
    real rad = radius_, edg = edge_;
    real bot = bot_, top = top_;

    if ( opt.set(rad, "diameter") )
        rad *= 0.5;
    else opt.set(rad, "radius");

#if SMOOTH_CYLINDER
    if ( opt.set(edg, "edge") )
    {
        if ( edg < 0 )
            throw InvalidParameter("cylinderZ:edge must be >= 0");
        if ( rad < edg )
            throw InvalidParameter("cylinderZ:edge must be <= radius");
    }
#endif
    
    if ( rad < 0 )
        throw InvalidParameter("cylinderZ:radius must be >= 0");
    
    real suf = 0;
    if ( opt.set(suf, "surface") )
    {
        // The surface is affine with 'len' = S + ( len - 2*edg ) * L;
        real S = surface_smooth_cylinder(rad, 2*edg, edg);
        real L = ( 2 * M_PI ) * rad;
        if ( suf < S )
            throw InvalidParameter("cylinderZ:surface must be > "+std::to_string(S));
        // back calculate length to match the given surface:
        real len = 2 * edg + ( suf - S ) / L;
        //std::clog << "CylinderZ surface " << surface_cylinder(rad, len, edg) << " " << suf << '\n';
        std::clog << "CylinderZ:length=" << len << '\n';
        bot = -0.5 * len;
        top =  0.5 * len;
    }
    else if ( opt.set(top, "length") )
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
    if ( top - bot < 2 * edg )
        throw InvalidParameter("cylinderZ:edge must be <= length / 2");

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


Vector SpaceCylinderZ::place() const
{
#if SMOOTH_CYLINDER
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


Vector SpaceCylinderZ::normalToEdge(Vector const& pos) const
{
#if ( DIM > 2 )
    const real R = std::sqrt(pos.XX * pos.XX + pos.YY * pos.YY);
    real tZ = top_ - pos.ZZ;
    real bZ = pos.ZZ - bot_;
    real dZ = min_real(abs_real(tZ-edge_), abs_real(bZ+edge_));
#if SMOOTH_CYLINDER
    if (( bZ < edge_ ) | ( tZ < edge_ ))
    {
        if ( R < radius_-edge_ )
            return Vector(0, 0, std::copysign(1, bZ-tZ));
        // project on the inner cylinder:
        const real n = min_real(R, radius_-edge_) / R;
        const real X = n * pos.XX;
        const real Y = n * pos.YY;
        const real Z = max_real(min_real(pos.ZZ, top_-edge_), bot_+edge_);
        // this will fail if 'pos' is already inside the inner cylinder:
        return normalize(pos-Vector(X,Y,Z));
    }
#endif
    if ( abs_real(R-radius_) < dZ )
        return Vector(pos.XX/R, pos.YY/R, 0);
    else
        return Vector(0, 0, std::copysign(1, bZ-tZ));
#endif
    return Vector(0, 0, 0);  // intentionally invalid!
}


real SpaceCylinderZ::surface() const
{
#if SMOOTH_CYLINDER
    return surface_smooth_cylinder(radius_, top_-bot_, edge_);
#else
    return surface_cylinder(radius_, top_-bot_);
#endif
}


/**
Pappus's (1st) Centroid Theorem: the volume of a planar area of revolution is
the product of the area A and the length of the path traced by its centroid,
For a half-circle, this is GC = 2 R / M_PI

https://en.wikipedia.org/wiki/Pappus%27s_centroid_theorem
*/
Vector SpaceCylinderZ::placeOnEdge(real) const
{
#if ( DIM > 2 )
#if SMOOTH_CYLINDER
    const real RE = radius_ - edge_;
    const real LE = top_ - bot_ - 2 * edge_;
    const real GC = RE + edge_ * ( 2.0 / M_PI );
    // surface elements divided by 2 * M_PI:
    const real S0 = radius_ * LE;
    const real S1 = square(RE);
    const real S2 = M_PI * GC * edge_;
    const real P = RNG.preal() * ( S0 + S1 + S2 );
    if ( P <= S0 )
    {
        Vector2 XY = Vector2::randU(radius_);
        real Z = bot_ + edge_ + LE * RNG.preal();
        return Vector(XY.XX, XY.YY, Z);
    }
    else if ( P < S0+S1 )
    {
        Vector2 XY = Vector2::randB(RE);
        real Z = RNG.choice(bot_, top_);
        return Vector(XY.XX, XY.YY, Z);
    }
    else
    {
        real Z, R, N;
        do {
            // generate a point inside Cylinder
            R = RE / edge_;
            /* generate R randomly between 0 and 1, with skewed probabilities:
             The probability is 'RE' at 0 and 'RE+edge_' at 1, as needed
             to obtain a uniform volume sampling of the cut cylinder,
             for RE < R < RE + edge_.
             This formula was derived by inverting the cumulative probability */
            R = std::sqrt(square(R)+RNG.preal()*(2*R+1)) - R;
            Z = RNG.sreal();
            // repeat until point is inside Torus:
            N = square(R) + square(Z);
        } while ( N > 1.0 );
        // normalize to get a point on the surface:
        N = edge_ / std::sqrt(N);
        Vector2 XY = Vector2::randU(RE+R*N);
        Z = N * Z + sign_select(Z, bot_+edge_, top_-edge_);
        return Vector(XY.XX, XY.YY, Z);
    }
#else
    // surface elements divided by M_PI * square(radius_)
    const real S0 = top_ - bot_;
    const real S1 = 2;
    const real P = RNG.preal() * ( S0 + S1 );
    if ( P < S0 )
    {
        Vector2 XY = Vector2::randU(radius_);
        real Z = bot_ + ( top_ - bot_ ) * RNG.preal();
        return Vector(XY.XX, XY.YY, Z);
    }
    else
    {
        Vector2 XY = Vector2::randB(radius_);
        real Z = RNG.choice(bot_, top_);
        return Vector(XY.XX, XY.YY, Z);
    }
#endif
#endif
    return Vector(0, 0, 0);  // intentionally invalid!
}


/**
 Pappus's (2nd) Centroid Theorem: the volume of a planar area of revolution is
 the product of the area A and the length of the path traced by its centroid R,
 For a half-disc, this is GC = 4 / 3.0 * R / M_PI
 
 https://en.wikipedia.org/wiki/Pappus%27s_centroid_theorem
 */
real SpaceCylinderZ::volume() const
{
#if SMOOTH_CYLINDER
    const real RE = radius_ - edge_;
    const real LE = top_ - bot_ - 2 * edge_;
    const real GC = RE + 4 / 3.0 * edge_ / M_PI;
    return M_PI * ( LE * square(radius_)       // central cylindrical part
                   + 2 * edge_ * square(RE)    // top and bottom central cylinders
                   + M_PI * GC * square(edge_) );     // revolution
#else
    return M_PI * ( top_ - bot_ ) * square(radius_);
#endif
}


bool SpaceCylinderZ::inside(Vector const& W) const
{
#if ( DIM > 2 )
    const real RT = W.normXYSqr();
# if SMOOTH_CYLINDER
    const real R = max_real(0, std::sqrt(RT)-radius_+edge_);
    const real Z = max_real(0, std::max(bot_+edge_-W.ZZ, W.ZZ-top_+edge_));
    return ( R*R + Z*Z <= edgeSqr_ );
# else
    return (( bot_ <= W.ZZ ) & ( W.ZZ <= top_ ) & ( RT <= square(radius_) ));
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
    const real RT = W.normXYSqr();
# if SMOOTH_CYLINDER
    const real E = edge_ + rad;
    const real R = max_real(0, std::sqrt(RT)-E);
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

//------------------------------------------------------------------------------
Vector SpaceCylinderZ::project(Vector const& W) const
{
    Vector P(W);
#if ( DIM > 2 )
# if SMOOTH_CYLINDER
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
# if SMOOTH_CYLINDER
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
void SpaceCylinderZ::setConfinement(Vector const& pos, Mecapoint const& mp, Meca& meca, 
                                    real stiff, real R, real B, real T)
{
#if ( DIM >= 3 )
    bool cap = false;
    bool cyl = false;
    real Z;

    // in top or bottom half
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
        //if ( abs_real(pos.ZZ-Z) > R - std::sqrt(dis) )
        if ( dis > square( R - abs_real(pos.ZZ-Z) ) )
            cyl = true;
        else
            cap = true;
    }
    
    if ( cap )
        meca.addPlaneClampZ(mp, Z, stiff);
    
    if ( cyl )
        meca.addCylinderClampZ(mp, R, stiff);
#endif
}


void SpaceCylinderZ::setConfinement(Vector const& pos, Mecapoint const& mp, Meca& meca,
                                    real stiff, real R, real B, real T, real E)
{
#if ( DIM >= 3 )
    real n = pos.normXY();
    bool in = true;
    real Z;
    
    if ( 2 * pos.ZZ - B > T )
    {
        // in top half, Z = top end
        Z = T-E;
        if ( pos.ZZ > Z )
        {
            // above top end
            in = false;
            if ( n < R-E ) {
                meca.addPlaneClampZ(mp, T, stiff);
                return;
            }
        }
    }
    else
    {
        //in bottom half, Z = bottom end
        Z = B+E;
        if ( pos.ZZ < Z )
        {
            in = false;
            if ( n < R-E ) {
                meca.addPlaneClampZ(mp, B, stiff);
                return;
            }
        }
    }
    
    if ( in )
    {
        if ( n > R-E  ||  R-E-n < abs_real(pos.ZZ-Z) )
            meca.addCylinderClampZ(mp, R, stiff);
        else
            meca.addPlaneClampZ(mp, Z, stiff);
        return;
    }
    
    /*
     we add interaction to the tangent plane, which is okay
     if the edge radius is not too small
     */
    real S = ( R - E ) / n;
    Vector3 rim(S*pos.XX, S*pos.YY, Z);  // on the reduced cylinder rim
    Vector dir = normalize(pos-rim);
    meca.addPlaneClamp(mp, rim+E*dir, dir, stiff);

#endif
}

/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCylinderZ::setConfinement(Vector const& pos, Mecapoint const& mp, Meca& meca, real stiff) const
{
#if SMOOTH_CYLINDER
    setConfinement(pos, mp, meca, stiff, radius_, bot_, top_, edge_);
#else
    setConfinement(pos, mp, meca, stiff, radius_, bot_, top_);
#endif
}

/**
 This applies the correct forces in the cylindrical and spherical parts.
 */
void SpaceCylinderZ::setConfinement(Vector const& pos, Mecapoint const& mp, real rad, Meca& meca, real stiff) const
{
    real R = max_real(0, radius_ - rad);
    real T = top_ - rad;
    real B = bot_ + rad;
    
    if ( B > T )
    {
        B = 0.5 * ( top_ + bot_ );
        T = B;
    }
    
#if SMOOTH_CYLINDER
    setConfinement(pos, mp, meca, stiff, R, B, T, edge_);
#else
    setConfinement(pos, mp, meca, stiff, R, B, T);
#endif
}


//------------------------------------------------------------------------------

void SpaceCylinderZ::write(Outputter& out) const
{
    writeMarker(out, TAG);
    writeShape(out, "RBTE");
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
    update();
}


void SpaceCylinderZ::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "RBTE");
    setLengths(len);
}

//------------------------------------------------------------------------------
#pragma mark - OpenGL display

#ifdef DISPLAY

#include "gle.h"
#include "gym_flute.h"
#include "gym_draw.h"

void SpaceCylinderZ::draw3D() const
{
    const float T(top_), TE(top_-edge_);
    const float B(bot_), BE(bot_+edge_);
    const float R(radius_), RE(radius_-edge_);
    const float E(edge_);
    
    size_t cnt = 2 * ( gle::pi_once + 3 );
    
    for ( size_t u = 0; u < gle::pi_twice; ++u )
    {
        flute6 * flu = gym::mapBufferV3N3(cnt);
        flute6 * ptr = flu;
        float CU = gle::cos_(u), CL = gle::cos_(u+1);
        float SU = gle::sin_(u), SL = gle::sin_(u+1);

        *ptr++ = {0, 0, T, 0, 0, 1};
#if SMOOTH_CYLINDER
        if ( edge_ > 0 )
        {
            //draw top arc
            for ( size_t j = 0; j <= gle::pi_half; ++j )
            {
                float C = gle::cos_(j), S = gle::sin_(j);
                float RS = RE + E*S;
                ptr[0] = {CU*RS, SU*RS, TE+E*C, CU*S, SU*S, C};
                ptr[1] = {CL*RS, SL*RS, TE+E*C, CL*S, SL*S, C};
                ptr += 2;
            }
            if ( 0 )
            {
                // at pi_half, C = 0 and S = 1
                ptr[0] = {CU*R, SU*R, TE, CU, SU, 0};
                ptr[1] = {CL*R, SL*R, TE, CL, SL, 0};
                ptr[2] = {CU*R, SU*R, BE, CU, SU, 0};
                ptr[3] = {CL*R, SL*R, BE, CL, SL, 0};
                ptr += 4;
            }
            //draw bottom arc
            for ( size_t j = gle::pi_half; j <= gle::pi_once; ++j )
            {
                float C = gle::cos_(j), S = gle::sin_(j);
                float RS = RE + E*S;
                ptr[0] = {CU*RS, SU*RS, BE+E*C, CU*S, SU*S, C};
                ptr[1] = {CL*RS, SL*RS, BE+E*C, CL*S, SL*S, C};
                ptr += 2;
            }
        }
        else
#endif
        {
            ptr[0] = {CU*R, SU*R, T, CU, SU, 0};
            ptr[1] = {CL*R, SL*R, T, CL, SL, 0};
            ptr[2] = {CU*R, SU*R, B, CU, SU, 0};
            ptr[3] = {CL*R, SL*R, B, CL, SL, 0};
            ptr += 4;
        }
        *ptr++ = {0, 0, B, 0, 0, -1};
        gym::unmapBufferV3N3();
        gym::drawTriangleStrip(0, ptr-flu);
        gym::cleanup();
    }
}

#else

void SpaceCylinderZ::draw3D() const {}

#endif
