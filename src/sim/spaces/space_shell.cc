// Cytosim -- Spherical shell space
#include "dim.h"
#include "space_shell.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include "random.h"
#include "meca.h"
#include <cmath>

SpaceShell::SpaceShell(SpaceProp const* p)
: Space(p), r_in_(0), r_out_(0)
{
    if ( DIM < 3 )
        throw InvalidParameter("spherical_shell is only valid in 3D");
}


void SpaceShell::resize(Glossary& opt)
{
    real r_in = r_in_, r_out = r_out_;
    real val = 0;
    bool got_inner = false;
    bool got_outer = false;

    if ( opt.set(val, "diameter") )
    {
        r_out = 0.5 * val;
        got_outer = true;
    }
    if ( opt.set(val, "radius") )
    {
        r_out = val;
        got_outer = true;
    }
    if ( opt.set(r_out, "outer") || opt.set(r_out, "radius_outer") || opt.set(r_out, "outer_radius") )
        got_outer = true;
    if ( opt.set(r_in, "inner") || opt.set(r_in, "radius_inner") || opt.set(r_in, "inner_radius") )
        got_inner = true;

    real width = -1;
    bool got_width = opt.set(width, "thickness") || opt.set(width, "width");
    if ( got_width )
    {
        if ( width <= 0 )
            throw InvalidParameter(prop->name()+":thickness must be > 0");
        if ( got_inner && !got_outer )
            r_out = r_in + width;
        else
            r_in = r_out - width;
    }

    if ( !( r_out > r_in && r_in >= 0 ) )
        throw InvalidParameter(prop->name()+": require outer>inner>=0");

    r_in_ = r_in;
    r_out_ = r_out;
}


void SpaceShell::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-r_out_, -r_out_, -r_out_);
    sup.set( r_out_,  r_out_,  r_out_);
}


real SpaceShell::volume() const
{
    return ( 4.0 * M_PI / 3.0 ) * ( cube(r_out_) - cube(r_in_) );
}


real SpaceShell::surface() const
{
    return ( 4.0 * M_PI ) * ( square(r_out_) + square(r_in_) );
}


bool SpaceShell::inside(Vector const& pos) const
{
    const real n = pos.normSqr();
    return ( n >= square(r_in_) - REAL_EPSILON ) && ( n <= square(r_out_) + REAL_EPSILON );
}


bool SpaceShell::allInside(Vector const& pos, const real rad) const
{
    assert_true( rad >= 0 );
    const real R = std::sqrt(pos.normSqr());
    return ( R >= r_in_ + rad ) && ( R <= r_out_ - rad );
}


Vector SpaceShell::place() const
{
    const real u = cube(r_in_) + RNG.preal() * ( cube(r_out_) - cube(r_in_) );
    return Vector::randU(std::cbrt(u));
}


Vector SpaceShell::normalToEdge(Vector const& pos) const
{
    const real n = pos.normSqr();
    if ( n <= REAL_EPSILON )
        return Vector(-1, 0, 0);

    const real R = std::sqrt(n);
    const Vector dir = pos * ( 1.0 / R );
    return ( R - r_in_ <= r_out_ - R ) ? -dir : dir;
}


Vector SpaceShell::placeOnEdge(real) const
{
    const real inner_area = square(r_in_);
    const real outer_area = square(r_out_);
    if ( RNG.preal() * ( inner_area + outer_area ) < inner_area )
        return Vector::randU(r_in_);
    return Vector::randU(r_out_);
}


Vector SpaceShell::project(Vector const& pos) const
{
    const real n = pos.normSqr();
    if ( n <= REAL_EPSILON )
        return Vector::randU(r_in_ > 0 ? r_in_ : r_out_);

    const real R = std::sqrt(n);
    real target = r_out_;
    if ( R < r_in_ )
        target = r_in_;
    else if ( R > r_out_ )
        target = r_out_;
    else
        target = ( R - r_in_ <= r_out_ - R ) ? r_in_ : r_out_;

    return pos * ( target / R );
}


void SpaceShell::clampToRadius(Vector const& pos, Mecapoint const& mp, Meca& meca,
                               const real target, const real stiff) const
{
    if ( pos.normSqr() > REAL_EPSILON )
        meca.addSphereClamp(pos, mp, Vector(0,0,0), target, stiff);
    else
        meca.addPointClamp(mp, Vector(target, 0, 0), stiff);
}


void SpaceShell::setConfinement(Vector const& pos, Mecapoint const& mp,
                                Meca& meca, real stiff) const
{
    const real R = std::sqrt(pos.normSqr());
    real target = r_out_;
    if ( R < r_in_ )
        target = r_in_;
    else if ( R > r_out_ )
        target = r_out_;
    else
        target = ( R - r_in_ <= r_out_ - R ) ? r_in_ : r_out_;
    clampToRadius(pos, mp, meca, target, stiff);
}


void SpaceShell::setConfinement(Vector const& pos, Mecapoint const& mp,
                                real rad, Meca& meca, real stiff) const
{
    assert_true( rad >= 0 );
    real inner = r_in_ + rad;
    real outer = r_out_ - rad;
    if ( inner > outer )
    {
        inner = 0.5 * ( r_in_ + r_out_ );
        outer = inner;
        std::cerr << "object is too big to fit in SpaceShell\n";
    }

    const real R = std::sqrt(pos.normSqr());
    real target = outer;
    if ( R < inner )
        target = inner;
    else if ( R > outer )
        target = outer;
    else
        target = ( R - inner <= outer - R ) ? inner : outer;
    clampToRadius(pos, mp, meca, target, stiff);
}


void SpaceShell::write(Outputter& out) const
{
    writeMarker(out, TAG);
    writeShape(out, "OI");
    out.writeUInt16(2);
    out.writeFloat(r_out_);
    out.writeFloat(r_in_);
}


void SpaceShell::setLengths(const real len[8])
{
    r_out_ = len[0];
    r_in_ = len[1];
    if ( !( r_out_ > r_in_ && r_in_ >= 0 ) )
        throw InvalidParameter("spherical_shell:setLengths: require outer>inner>=0");
}


void SpaceShell::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "OI");
    setLengths(len);
}


#ifdef DISPLAY

#include "gym_draw.h"
#include "gym_flute.h"
#include "gym_flute_dim.h"

namespace {
    static void drawCircleXY(real r, size_t cnt, float w)
    {
        if ( r <= REAL_EPSILON ) return;
        size_t i = 0;
        fluteD* pts = gym::mapBufferVD(cnt + 1);
        for ( size_t k = 0; k <= cnt; ++k )
        {
            const real a = ( 2.0 * M_PI * k ) / cnt;
            pts[i++] = Vector( r * std::cos(a), r * std::sin(a), 0 );
        }
        gym::unmapBufferVD();
        gym::drawLineStrip(w, 0, i);
    }

    static void drawCircleXZ(real r, size_t cnt, float w)
    {
        if ( r <= REAL_EPSILON ) return;
        size_t i = 0;
        fluteD* pts = gym::mapBufferVD(cnt + 1);
        for ( size_t k = 0; k <= cnt; ++k )
        {
            const real a = ( 2.0 * M_PI * k ) / cnt;
            pts[i++] = Vector( r * std::cos(a), 0, r * std::sin(a) );
        }
        gym::unmapBufferVD();
        gym::drawLineStrip(w, 0, i);
    }

    static void drawCircleYZ(real r, size_t cnt, float w)
    {
        if ( r <= REAL_EPSILON ) return;
        size_t i = 0;
        fluteD* pts = gym::mapBufferVD(cnt + 1);
        for ( size_t k = 0; k <= cnt; ++k )
        {
            const real a = ( 2.0 * M_PI * k ) / cnt;
            pts[i++] = Vector( 0, r * std::cos(a), r * std::sin(a) );
        }
        gym::unmapBufferVD();
        gym::drawLineStrip(w, 0, i);
    }
}

void SpaceShell::draw2D(float w) const
{
    drawCircleXY(r_out_, 128, w);
    drawCircleXY(r_in_, 128, w);
}


void SpaceShell::draw3D() const
{
    const size_t cnt = 128;
    const float w = 1.0f;
    drawCircleXY(r_out_, cnt, w);
    drawCircleXZ(r_out_, cnt, w);
    drawCircleYZ(r_out_, cnt, w);
    drawCircleXY(r_in_, cnt, w);
    drawCircleXZ(r_in_, cnt, w);
    drawCircleYZ(r_in_, cnt, w);
}

#else

void SpaceShell::draw2D(float) const {}
void SpaceShell::draw3D() const {}

#endif
