// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "dim.h"
#include "space_disc.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include "random.h"
#include "meca.h"


SpaceDisc::SpaceDisc(SpaceProp const* p)
: Space(p)
{
    if ( DIM != 3 )
        throw InvalidParameter("disc is only usable in 3D");
    radius_ = 0;
    bot_ = 0;
    top_ = 0;
    mid_ = 0;
}


void SpaceDisc::resize(Glossary& opt)
{
    real rad = radius_;
    real bot = bot_, top = top_;

    if ( opt.set(rad, "diameter") )
        rad *= 0.5;
    else opt.set(rad, "radius");

    if ( rad < 0 )
        throw InvalidParameter("disc:radius must be >= 0");
    
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
    
    radius_ = rad;
    bot_ = bot;
    top_ = top;
    mid_ = ( top_ + bot_ ) * 0.5;
}


void SpaceDisc::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-radius_,-radius_, bot_);
    sup.set( radius_, radius_, top_);
}


real SpaceDisc::volume() const
{
#if ( DIM <= 2 )
    return M_PI * square(radius_);
#else
    return M_PI * ( top_ - bot_ ) * square(radius_);
#endif
}


bool SpaceDisc::inside(Vector const& pos) const
{
#if ( DIM <= 2 )
    return pos.normSqr() <= square(radius_);
#else
    return ( bot_ <= pos.ZZ ) & ( pos.ZZ <= top_ );
#endif
}


bool SpaceDisc::allInside(Vector const& pos, const real rad) const
{
    assert_true( rad >= 0 );
#if ( DIM <= 2 )
    real R = max_real(0, radius_-rad);  // remaining radius
    return pos.normSqr() <= square(R);
#else
    return ( bot_ + rad <= pos.ZZ ) & ( pos.ZZ + rad <= top_ );
#endif
}


Vector SpaceDisc::project(Vector const& pos) const
{
#if ( DIM <= 2 )
    real n = pos.normSqr();
    if ( n > 0 ) {
        return pos * ( radius_ / std::sqrt(n) );
    }
    else {
        //select a random point on the surface
        return radius_ * Vector::randU();
    }
#else
    return Vector(pos.XX, pos.YY, bot_);
#endif
}


/**
 The implementation here is not physically exact, as we only consider the
 radial component, and leads to artifacts in particular for small radius
 @todo Implement correct bounce in Disc/Sphere
 */
Vector SpaceDisc::bounce(Vector const& pos) const
{
#if ( DIM >= 3 )
    real X = pos.XX;
    real Y = pos.YY;
    real n = X*X + Y*Y;
    if ( n > square(radius_) )
    {
        n = 2 * ( radius_ / std::sqrt(n) ) - 1;
        X *= n;
        Y *= n;
    }
    real Z = bounce1(pos.ZZ, bot_, top_-bot_);
    return Vector(X, Y, Z);
#endif
    return pos;
}


Vector SpaceDisc::place() const
{
    const Vector2 V = Vector2::randB(radius_);
    return Vector(V.XX, V.YY, bot_+RNG.preal()*(top_-bot_));
}


Vector SpaceDisc::placeOnEdge(real) const
{
    const Vector2 V = Vector2::randB(radius_);
    return Vector(V.XX, V.YY, bot_);
}


void SpaceDisc::setConfinement(Vector const& pos, Mecapoint const& mp,
                               Meca& meca, real stiff) const
{
# if ( DIM <= 2 )
    meca.addSphereClamp(pos, mp, Vector(0,0,0), radius_, stiff);
#else
    real Z = sign_select(pos.ZZ - mid_, bot_, top_);
    meca.addPlaneClampZ(mp, Z, stiff);
#endif
}


void SpaceDisc::setConfinement(Vector const& pos, Mecapoint const& mp,
                               real rad, Meca& meca, real stiff) const
{
# if ( DIM <= 2 )
    if ( radius_ > rad )
    {
        meca.addSphereClamp(pos, mp, Vector(0,0,0), radius_-rad, stiff);
    }
    else {
        meca.addPointClamp(mp, Vector(0,0,0), stiff);
        std::cerr << "object is too big to fit in SpaceDisc\n";
    }
#else
    real Z = sign_select(pos.ZZ - mid_, bot_+rad, top_-rad);
    meca.addPlaneClampZ(mp, Z, stiff);
#endif
}


void SpaceDisc::write(Outputter& out) const
{
    writeMarker(out, TAG);
    writeShape(out, "RBT");
    out.writeUInt16(4);
    out.writeFloat(radius_);
    out.writeFloat(bot_);
    out.writeFloat(top_);
    out.writeFloat(mid_);
}


void SpaceDisc::setLengths(const real len[8])
{
    radius_ = len[0];
    bot_    = len[1];
    top_    = len[2];
    mid_ = ( top_ + bot_ ) * 0.5;
}

void SpaceDisc::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "RBT");
    setLengths(len);
}


#ifdef DISPLAY

#include "gle.h"
#include "gym_view.h"

void SpaceDisc::draw2D(float width) const
{
    gle::circle(radius_, width);
}

void SpaceDisc::draw3D() const
{
    const float R(radius_);
    const float B(bot_);

    gym::transScale(0, 0, -B, R);
    //gle::tube1();
    gle::disc1();
    //gle::discTop1();
}

#else

void SpaceDisc::draw2D(float) const {}
void SpaceDisc::draw3D() const {}

#endif
