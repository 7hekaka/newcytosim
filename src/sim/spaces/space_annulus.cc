// Cytosim — Annulus / Annular slab
#include "space_annulus.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include <cmath>
#ifdef DISPLAY
#include "gym_draw.h"
#include "gym_flute.h"
#include "gym_flute_dim.h"
#endif

SpaceAnnulus::SpaceAnnulus(SpaceProp const* p) : Space(p)
{
    // 2D and 3D are both valid
    r_in_ = 0;
    r_out_ = 0;
#if (DIM >= 3)
    z_bot_ = 0;
    z_top_ = 0;
#endif
}

void SpaceAnnulus::resize(Glossary& opt)
{
    real r_out = r_out_, r_in = r_in_;
    real rad = -1, width = -1, r_outer = -1, r_inner = -1;

    // accept multiple parameter styles
    if ( opt.set(rad,    "radius") ) { /* handled below */ }
    opt.set(width,       "width");
    opt.set(r_outer,     "radius_outer");
    opt.set(r_inner,     "radius_inner");
    opt.set(r_out,       "outer");
    opt.set(r_in,        "inner");

#if (DIM >= 3)
    opt.set(z_bot_, "bottom");
    opt.set(z_top_, "top");
    if ( z_top_ < z_bot_ ) std::swap(z_bot_, z_top_);
#endif

    if ( r_outer > 0 ) r_out = r_outer;
    if ( r_inner >= 0 ) r_in = r_inner;

    if ( rad > 0 ) {
        r_out = rad;
        if ( width > 0 )
            r_in = std::max<real>(0, r_out - width);
    }

    if ( !(r_out > r_in && r_in >= 0) )
        throw InvalidParameter("annulus: require outer>inner>=0");

    r_out_ = r_out;
    r_in_  = r_in;
}

void SpaceAnnulus::boundaries(Vector& inf, Vector& sup) const
{
#if (DIM >= 3)
    inf.set(-r_out_, -r_out_, z_bot_);
    sup.set( r_out_,  r_out_, z_top_);
#else
    inf.set(-r_out_, -r_out_, 0);
    sup.set( r_out_,  r_out_, 0);
#endif
}

bool SpaceAnnulus::inside(Vector const& W) const
{
    const real r = std::sqrt(W.XX*W.XX + W.YY*W.YY);
#if (DIM >= 3)
    if ( W.ZZ < z_bot_ || W.ZZ > z_top_ ) return false;
#endif
    return (r >= r_in_ - REAL_EPSILON) && (r <= r_out_ + REAL_EPSILON);
}

bool SpaceAnnulus::allInside(Vector const& W, const real rad) const
{
    assert_true( rad >= 0 );
    const real r = std::sqrt(W.XX*W.XX + W.YY*W.YY);
#if (DIM >= 3)
    if ( W.ZZ < z_bot_ + rad || W.ZZ > z_top_ - rad ) return false;
#endif
    return (r >= r_in_ + rad) && (r <= r_out_ - rad);
}

Vector SpaceAnnulus::project(Vector const& W) const
{
    Vector P(W);

    // --- radial clamp/projection in XY ---
    real r = std::sqrt(W.XX*W.XX + W.YY*W.YY);
    real ux = 1, uy = 0;
    if ( r > REAL_EPSILON ) { ux = W.XX / r; uy = W.YY / r; }

    real r_target;
    if ( r < r_in_ )
        r_target = r_in_;
    else if ( r > r_out_ )
        r_target = r_out_;
    else
        r_target = (r - r_in_ <= r_out_ - r) ? r_in_ : r_out_;

    P.XX = ux * r_target;
    P.YY = uy * r_target;

    // --- z clamp (3D slab) ---
#if (DIM >= 3)
    if ( W.ZZ < z_bot_ )      P.ZZ = z_bot_;
    else if ( W.ZZ > z_top_ ) P.ZZ = z_top_;
    else                      P.ZZ = W.ZZ;
#endif
    return P;
}

real SpaceAnnulus::volume() const
{
    const real area = M_PI * ( r_out_*r_out_ - r_in_*r_in_ );
#if (DIM >= 3)
    return area * ( z_top_ - z_bot_ );
#else
    return area; // 2D: area of the annulus (see Space docs)
#endif
}

real SpaceAnnulus::surface() const
{
#if (DIM >= 3)
    const real h = ( z_top_ - z_bot_ );
    const real sides = 2 * M_PI * (r_out_ + r_in_) * h;
    const real caps  = 2 * M_PI * ( r_out_*r_out_ - r_in_*r_in_ );
    return sides + caps;
#else
    return 2 * M_PI * ( r_out_ + r_in_ ); // 2D: total perimeter (two rims)
#endif
}

void SpaceAnnulus::write(Outputter& out) const
{
    writeMarker(out, TAG);
    // Layout tag is only checked against itself when reading (“expected” string).
    // We encode: outer, inner, [bottom, top]
#if (DIM >= 3)
    writeShape(out, "OIBT"); // Outer, Inner, Bottom, Top
    out.writeUInt16(4);
    out.writeFloat(r_out_);
    out.writeFloat(r_in_);
    out.writeFloat(z_bot_);
    out.writeFloat(z_top_);
#else
    writeShape(out, "OI");   // Outer, Inner
    out.writeUInt16(2);
    out.writeFloat(r_out_);
    out.writeFloat(r_in_);
#endif
}

void SpaceAnnulus::setLengths(const real len[8])
{
    // Order: outer, inner, bottom, top
    r_out_ = len[0];
    r_in_  = len[1];
#if (DIM >= 3)
    z_bot_ = len[2];
    z_top_ = len[3];
    if ( z_top_ < z_bot_ ) std::swap(z_bot_, z_top_);
#endif
    if ( !(r_out_ > r_in_ && r_in_ >= 0) )
        throw InvalidParameter("annulus:setLengths: require outer>inner>=0");
}

void SpaceAnnulus::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
#if (DIM >= 3)
    readShape(in, 8, len, "OIBT");
#else
    readShape(in, 8, len, "OI");
#endif
    setLengths(len);
}

// space_annulus.cc
// ...

#ifdef DISPLAY
namespace {
    // helper: draw a circle (line strip) in the XY plane at height z
    static void drawRingAtZ(real z, real r, size_t cnt, float w)
    {
        if ( r <= REAL_EPSILON ) return;               // nothing to draw if degenerate
        size_t i = 0;
        fluteD* pts = gym::mapBufferVD(cnt + 1);       // +1 to close the loop
        for ( size_t k = 0; k <= cnt; ++k )
        {
            real a = ( 2.0 * M_PI * k ) / cnt;
            pts[i++] = Vector( r * std::cos(a), r * std::sin(a), z );
        }
        gym::unmapBufferVD();
        gym::drawLineStrip(w, 0, i);
    }
}

void SpaceAnnulus::draw2D(float w) const
{
    // 2D build: show both rims at z=0
    const size_t cnt = 128;
    drawRingAtZ(0.0, r_out_, cnt, w);
    drawRingAtZ(0.0, r_in_,  cnt, w);
}

void SpaceAnnulus::draw3D() const
{
#if (DIM >= 3)
    const size_t cnt = 128;
    const float  w   = 1.0f;

    // Only top & bottom outlines, both inner and outer:
    drawRingAtZ(z_bot_, r_out_, cnt, w);
    drawRingAtZ(z_bot_, r_in_,  cnt, w);
    drawRingAtZ(z_top_, r_out_, cnt, w);
    drawRingAtZ(z_top_, r_in_,  cnt, w);
#else
    draw2D(1.0f);
#endif
}
#else
void SpaceAnnulus::draw2D(float) const {}
void SpaceAnnulus::draw3D() const {}
#endif