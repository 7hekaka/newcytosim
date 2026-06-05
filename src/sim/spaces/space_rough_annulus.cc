// Cytosim -- annular slab with a static rough inner wall
#include "space_rough_annulus.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#ifdef DISPLAY
#include "gym_draw.h"
#include "gym_flute.h"
#include "gym_flute_dim.h"
#endif

namespace
{
    constexpr unsigned ROUGH_SEED_PACK = 100000;

    static real clamp_coord(real x, real lo, real hi)
    {
        return std::max(lo, std::min(x, hi));
    }

    static real square_real(real x)
    {
        return x * x;
    }

    static void radial_frame(Vector const& W, real& r, real& theta, real& ux, real& uy)
    {
        r = std::sqrt(W.XX*W.XX + W.YY*W.YY);
        if ( r > REAL_EPSILON )
        {
            ux = W.XX / r;
            uy = W.YY / r;
            theta = std::atan2(W.YY, W.XX);
        }
        else
        {
            ux = 1;
            uy = 0;
            theta = 0;
        }
    }

    static Vector radial_point(real ux, real uy, real radius, real z)
    {
        return Vector(ux * radius, uy * radius, z);
    }

    static uint32_t hash_u32(uint32_t x)
    {
        x ^= x >> 16;
        x *= 0x7feb352dU;
        x ^= x >> 15;
        x *= 0x846ca68bU;
        x ^= x >> 16;
        return x;
    }

    static real hash_unit(uint32_t seed, uint32_t channel)
    {
        return real(hash_u32(seed ^ (0x9e3779b9U * (channel + 1U))) & 0x00ffffffU) / real(0x01000000U);
    }

    static unsigned bounded_mode(real mode, unsigned fallback)
    {
        if ( mode < 0 )
            return fallback;
        return unsigned(std::max<real>(0, std::round(mode)));
    }

    static unsigned pack_seed_components(unsigned seed, unsigned components)
    {
        return ROUGH_SEED_PACK * std::min(components, 99U) + seed % ROUGH_SEED_PACK;
    }

    static void unpack_seed_components(real packed, unsigned& seed, unsigned& components)
    {
        const unsigned code = unsigned(std::max<real>(0, std::round(packed)));
        if ( code >= ROUGH_SEED_PACK )
        {
            components = std::max(1U, code / ROUGH_SEED_PACK);
            seed = code % ROUGH_SEED_PACK;
        }
        else
        {
            seed = code;
        }
    }
}


SpaceRoughAnnulus::SpaceRoughAnnulus(SpaceProp const* p) : Space(p)
{
    r_in_ = 0;
    r_out_ = 0;
#if (DIM >= 3)
    z_bot_ = 0;
    z_top_ = 0;
#endif
    amplitude_ = 0;
    theta_mode_ = 11;
    z_mode_ = 7;
    phase_ = 0;
    rough_seed_ = 17;
    rough_components_ = 12;
}


real SpaceRoughAnnulus::innerRadius(real theta, real z) const
{
    if ( amplitude_ <= 0 )
        return r_in_;

#if (DIM >= 3)
    const real height = z_top_ - z_bot_;
    const real z_phase = ( height > REAL_EPSILON ) ? ( 2.0 * M_PI * ( z - z_bot_ ) / height ) : 0;
#else
    const real z_phase = 0;
    (void)z;
#endif

    const unsigned theta_max = std::max(1U, bounded_mode(theta_mode_, 11));
    const unsigned theta_min = std::min(theta_max, 2U);
    const unsigned z_max = bounded_mode(z_mode_, 7);
    const unsigned count = std::max(1U, std::min(rough_components_, 64U));

    real sum = 0;
    real norm = 0;
    for ( unsigned k = 0; k < count; ++k )
    {
        const uint32_t base = rough_seed_ + 0x85ebca6bU * ( k + 1U );
        const unsigned mt = theta_min + unsigned(hash_unit(base, 0) * real(theta_max - theta_min + 1U));
        const unsigned mz = z_max ? unsigned(hash_unit(base, 1) * real(z_max + 1U)) : 0U;
        const real direction = ( hash_unit(base, 2) < 0.5 ) ? -1.0 : 1.0;
        const real component_phase = 2.0 * M_PI * hash_unit(base, 3) + phase_;
        const real weight = ( 0.65 + 0.70 * hash_unit(base, 4) ) / std::sqrt(real(mt + mz + 1U));
        sum += weight * std::sin(real(mt) * theta + direction * real(mz) * z_phase + component_phase);
        norm += weight * weight;
    }

    const real wave = sum / std::sqrt(std::max(norm, REAL_EPSILON));
    return r_in_ + amplitude_ * std::tanh(0.9 * wave);
}


void SpaceRoughAnnulus::validate() const
{
    if ( amplitude_ < 0 )
        throw InvalidParameter("rough_annulus:amplitude must be >= 0");
    if ( rough_components_ < 1 )
        throw InvalidParameter("rough_annulus:rough_components must be >= 1");
    if ( r_in_ - amplitude_ < 0 )
        throw InvalidParameter("rough_annulus: require inner>=amplitude");
    if ( !( r_out_ > r_in_ + amplitude_ && r_in_ >= 0 ) )
        throw InvalidParameter("rough_annulus: require outer>inner+amplitude and inner>=0");
#if (DIM >= 3)
    if ( z_top_ < z_bot_ )
        throw InvalidParameter("rough_annulus:bottom must be <= top");
#endif
}


void SpaceRoughAnnulus::resize(Glossary& opt)
{
    real r_out = r_out_, r_in = r_in_;
    real amp = amplitude_;
    real rad = -1, width = -1, r_outer = -1, r_inner = -1;

    if ( opt.set(rad, "radius") ) { /* handled below */ }
    opt.set(width,   "width");
    opt.set(r_outer, "radius_outer");
    opt.set(r_inner, "radius_inner");
    opt.set(r_out,   "outer");
    opt.set(r_in,    "inner");

    opt.set(amp, "amplitude") || opt.set(amp, "rough_amplitude") || opt.set(amp, "inner_amplitude");
    opt.set(theta_mode_, "theta_mode") || opt.set(theta_mode_, "angular_mode");
    opt.set(z_mode_, "z_mode") || opt.set(z_mode_, "axial_mode");
    opt.set(rough_seed_, "rough_seed") || opt.set(rough_seed_, "seed");
    opt.set(rough_components_, "rough_components") || opt.set(rough_components_, "components");
    opt.set(phase_, "phase");

#if (DIM >= 3)
    opt.set(z_bot_, "bottom");
    opt.set(z_top_, "top");
    if ( z_top_ < z_bot_ )
        std::swap(z_bot_, z_top_);
#endif

    if ( r_outer > 0 )
        r_out = r_outer;
    if ( r_inner >= 0 )
        r_in = r_inner;

    if ( rad > 0 )
    {
        r_out = rad;
        if ( width > 0 )
            r_in = std::max<real>(0, r_out - width);
    }

    r_out_ = r_out;
    r_in_ = r_in;
    amplitude_ = amp;
    validate();
}


void SpaceRoughAnnulus::boundaries(Vector& inf, Vector& sup) const
{
#if (DIM >= 3)
    inf.set(-r_out_, -r_out_, z_bot_);
    sup.set( r_out_,  r_out_, z_top_);
#else
    inf.set(-r_out_, -r_out_, 0);
    sup.set( r_out_,  r_out_, 0);
#endif
}


bool SpaceRoughAnnulus::inside(Vector const& W) const
{
    real r, theta, ux, uy;
    radial_frame(W, r, theta, ux, uy);
#if (DIM >= 3)
    if ( W.ZZ < z_bot_ || W.ZZ > z_top_ )
        return false;
    const real z = W.ZZ;
#else
    const real z = 0;
#endif
    return ( r >= innerRadius(theta, z) - REAL_EPSILON ) && ( r <= r_out_ + REAL_EPSILON );
}


bool SpaceRoughAnnulus::allInside(Vector const& W, const real rad) const
{
    assert_true( rad >= 0 );
    if ( !inside(W) )
        return false;
#if (DIM >= 3)
    if ( W.ZZ < z_bot_ + rad || W.ZZ > z_top_ - rad )
        return false;
#endif
    return ( project(W) - W ).normSqr() >= rad * rad;
}


Vector SpaceRoughAnnulus::project(Vector const& W) const
{
    real r, theta, ux, uy;
    radial_frame(W, r, theta, ux, uy);

#if (DIM >= 3)
    const real z_side = clamp_coord(W.ZZ, z_bot_, z_top_);
#else
    const real z_side = 0;
#endif
    Vector best = radial_point(ux, uy, innerRadius(theta, z_side), z_side);
    real best_d = ( W - best ).normSqr();

    Vector candidate = radial_point(ux, uy, r_out_, z_side);
    real dist = ( W - candidate ).normSqr();
    if ( dist < best_d )
    {
        best = candidate;
        best_d = dist;
    }

#if (DIM >= 3)
    const real r_bot = clamp_coord(r, innerRadius(theta, z_bot_), r_out_);
    candidate = radial_point(ux, uy, r_bot, z_bot_);
    dist = ( W - candidate ).normSqr();
    if ( dist < best_d )
    {
        best = candidate;
        best_d = dist;
    }

    const real r_top = clamp_coord(r, innerRadius(theta, z_top_), r_out_);
    candidate = radial_point(ux, uy, r_top, z_top_);
    dist = ( W - candidate ).normSqr();
    if ( dist < best_d )
        best = candidate;
#endif

    return best;
}


real SpaceRoughAnnulus::averageInnerRadiusSqr(size_t theta_cnt, size_t z_cnt) const
{
    real sum = 0;
    size_t cnt = 0;
    for ( size_t i = 0; i < theta_cnt; ++i )
    {
        const real theta = 2.0 * M_PI * ( real(i) + 0.5 ) / real(theta_cnt);
        for ( size_t j = 0; j < z_cnt; ++j )
        {
#if (DIM >= 3)
            const real z = z_bot_ + ( z_top_ - z_bot_ ) * ( real(j) + 0.5 ) / real(z_cnt);
#else
            const real z = 0;
#endif
            sum += square_real(innerRadius(theta, z));
            ++cnt;
        }
    }
    return sum / real(cnt);
}


real SpaceRoughAnnulus::volume() const
{
    const real area = M_PI * ( r_out_ * r_out_ - averageInnerRadiusSqr(128, 64) );
#if (DIM >= 3)
    return area * ( z_top_ - z_bot_ );
#else
    return area;
#endif
}


real SpaceRoughAnnulus::surface() const
{
    const size_t theta_cnt = 128;
    const real dtheta = 2.0 * M_PI / real(theta_cnt);

#if (DIM >= 3)
    const size_t z_cnt = 64;
    const real dz = ( z_top_ - z_bot_ ) / real(z_cnt);
    real inner_side = 0;
    for ( size_t i = 0; i < theta_cnt; ++i )
    {
        const real theta = 2.0 * M_PI * ( real(i) + 0.5 ) / real(theta_cnt);
        for ( size_t j = 0; j < z_cnt; ++j )
        {
            const real z = z_bot_ + ( z_top_ - z_bot_ ) * ( real(j) + 0.5 ) / real(z_cnt);
            const real r = innerRadius(theta, z);
            const real dr_dtheta = ( innerRadius(theta + dtheta, z) - innerRadius(theta - dtheta, z) ) / ( 2.0 * dtheta );
            const real dr_dz = ( innerRadius(theta, z + dz) - innerRadius(theta, z - dz) ) / ( 2.0 * dz );
            inner_side += std::sqrt(r*r * (1.0 + dr_dz*dr_dz) + dr_dtheta*dr_dtheta) * dtheta * dz;
        }
    }
    const real height = z_top_ - z_bot_;
    const real outer_side = 2.0 * M_PI * r_out_ * height;
    const real caps = 2.0 * M_PI * ( r_out_ * r_out_ - averageInnerRadiusSqr(theta_cnt, z_cnt) );
    return outer_side + inner_side + caps;
#else
    real inner_perimeter = 0;
    for ( size_t i = 0; i < theta_cnt; ++i )
    {
        const real theta = 2.0 * M_PI * ( real(i) + 0.5 ) / real(theta_cnt);
        const real r = innerRadius(theta, 0);
        const real dr_dtheta = ( innerRadius(theta + dtheta, 0) - innerRadius(theta - dtheta, 0) ) / ( 2.0 * dtheta );
        inner_perimeter += std::sqrt(r*r + dr_dtheta*dr_dtheta) * dtheta;
    }
    return 2.0 * M_PI * r_out_ + inner_perimeter;
#endif
}


void SpaceRoughAnnulus::write(Outputter& out) const
{
    writeMarker(out, TAG);
#if (DIM >= 3)
    writeShape(out, "OIBTAMZS");
    out.writeUInt16(8);
    out.writeFloat(r_out_);
    out.writeFloat(r_in_);
    out.writeFloat(z_bot_);
    out.writeFloat(z_top_);
    out.writeFloat(amplitude_);
    out.writeFloat(theta_mode_);
    out.writeFloat(z_mode_);
    out.writeFloat(real(pack_seed_components(rough_seed_, rough_components_)));
#else
    writeShape(out, "OIAMS");
    out.writeUInt16(5);
    out.writeFloat(r_out_);
    out.writeFloat(r_in_);
    out.writeFloat(amplitude_);
    out.writeFloat(theta_mode_);
    out.writeFloat(real(pack_seed_components(rough_seed_, rough_components_)));
#endif
}


void SpaceRoughAnnulus::setLengths(const real len[8])
{
    r_out_ = len[0];
    r_in_ = len[1];
#if (DIM >= 3)
    z_bot_ = len[2];
    z_top_ = len[3];
    amplitude_ = len[4];
    theta_mode_ = len[5];
    z_mode_ = len[6];
    phase_ = 0;
    unpack_seed_components(len[7], rough_seed_, rough_components_);
    if ( z_top_ < z_bot_ )
        std::swap(z_bot_, z_top_);
#else
    amplitude_ = len[2];
    theta_mode_ = len[3];
    phase_ = 0;
    unpack_seed_components(len[4], rough_seed_, rough_components_);
#endif
    validate();
}


void SpaceRoughAnnulus::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
#if (DIM >= 3)
    len[7] = pack_seed_components(rough_seed_, rough_components_);
    readShape(in, 8, len, "OIBTAMZS");
#else
    len[4] = pack_seed_components(rough_seed_, rough_components_);
    readShape(in, 8, len, "OIAMS");
#endif
    setLengths(len);
}


#ifdef DISPLAY
namespace
{
    static flute6 vertex_normal(real x, real y, real z, real nx, real ny, real nz)
    {
        return flute6(float(x), float(y), float(z), float(nx), float(ny), float(nz));
    }

    static void normalize(real& x, real& y, real& z)
    {
        const real n = std::sqrt(x*x + y*y + z*z);
        if ( n > REAL_EPSILON )
        {
            x /= n;
            y /= n;
            z /= n;
        }
        else
        {
            x = 1;
            y = 0;
            z = 0;
        }
    }

    static void drawRoughRingAtZ(SpaceRoughAnnulus const& spc, real z, size_t cnt, float w)
    {
        size_t i = 0;
        fluteD* pts = gym::mapBufferVD(cnt + 1);
        for ( size_t k = 0; k <= cnt; ++k )
        {
            const real a = ( 2.0 * M_PI * k ) / cnt;
            const real r = spc.innerRadius(a, z);
            pts[i++] = Vector(r * std::cos(a), r * std::sin(a), z);
        }
        gym::unmapBufferVD();
        gym::drawLineStrip(w, 0, i);
    }

    static void drawSmoothRingAtZ(real z, real r, size_t cnt, float w)
    {
        if ( r <= REAL_EPSILON )
            return;
        size_t i = 0;
        fluteD* pts = gym::mapBufferVD(cnt + 1);
        for ( size_t k = 0; k <= cnt; ++k )
        {
            const real a = ( 2.0 * M_PI * k ) / cnt;
            pts[i++] = Vector(r * std::cos(a), r * std::sin(a), z);
        }
        gym::unmapBufferVD();
        gym::drawLineStrip(w, 0, i);
    }

    static void innerWallNormal(SpaceRoughAnnulus const& spc, real theta, real z, real dtheta, real dz,
                                real& nx, real& ny, real& nz)
    {
        const real r = spc.innerRadius(theta, z);
        const real dr_dtheta = ( spc.innerRadius(theta + dtheta, z) - spc.innerRadius(theta - dtheta, z) ) / ( 2.0 * dtheta );
        const real dr_dz = ( spc.innerRadius(theta, z + dz) - spc.innerRadius(theta, z - dz) ) / ( 2.0 * dz );
        const real c = std::cos(theta);
        const real s = std::sin(theta);

        const real ptheta_x = dr_dtheta * c - r * s;
        const real ptheta_y = dr_dtheta * s + r * c;
        const real pz_x = dr_dz * c;
        const real pz_y = dr_dz * s;

        // ptheta x pz points toward increasing radius. The inner boundary normal
        // for the annular volume points toward the central hole, so flip it.
        nx = -ptheta_y;
        ny =  ptheta_x;
        nz = -( ptheta_x * pz_y - ptheta_y * pz_x );
        normalize(nx, ny, nz);
    }

    static void drawOuterWall(real z_bot, real z_top, real radius, size_t cnt)
    {
        if ( radius <= REAL_EPSILON )
            return;

        size_t i = 0;
        flute6* pts = gym::mapBufferV3N3(2 * (cnt + 1));
        for ( size_t k = 0; k <= cnt; ++k )
        {
            const real a = ( 2.0 * M_PI * k ) / cnt;
            const real c = std::cos(a);
            const real s = std::sin(a);
            pts[i++] = vertex_normal(radius * c, radius * s, z_top, c, s, 0);
            pts[i++] = vertex_normal(radius * c, radius * s, z_bot, c, s, 0);
        }
        gym::unmapBufferV3N3();
        gym::drawTriangleStrip(0, i);
        gym::cleanupVN();
    }

    static void drawInnerWall(SpaceRoughAnnulus const& spc, real z_bot, real z_top, size_t theta_cnt, size_t z_cnt)
    {
        const real height = z_top - z_bot;
        const real dtheta = 2.0 * M_PI / real(theta_cnt);
        const real dz = std::max(height / real(std::max<size_t>(z_cnt, 1)), REAL_EPSILON);

        for ( size_t j = 0; j < z_cnt; ++j )
        {
            const real z0 = z_bot + height * real(j) / real(z_cnt);
            const real z1 = z_bot + height * real(j + 1) / real(z_cnt);
            size_t i = 0;
            flute6* pts = gym::mapBufferV3N3(2 * (theta_cnt + 1));
            for ( size_t k = 0; k <= theta_cnt; ++k )
            {
                const real a = ( 2.0 * M_PI * k ) / theta_cnt;
                const real c = std::cos(a);
                const real s = std::sin(a);
                real nx, ny, nz;

                const real r0 = spc.innerRadius(a, z0);
                innerWallNormal(spc, a, z0, dtheta, dz, nx, ny, nz);
                pts[i++] = vertex_normal(r0 * c, r0 * s, z0, nx, ny, nz);

                const real r1 = spc.innerRadius(a, z1);
                innerWallNormal(spc, a, z1, dtheta, dz, nx, ny, nz);
                pts[i++] = vertex_normal(r1 * c, r1 * s, z1, nx, ny, nz);
            }
            gym::unmapBufferV3N3();
            gym::drawTriangleStrip(0, i);
            gym::cleanupVN();
        }
    }

    static void drawTopCap(SpaceRoughAnnulus const& spc, real z, real outer, size_t cnt)
    {
        size_t i = 0;
        flute6* pts = gym::mapBufferV3N3(2 * (cnt + 1));
        for ( size_t k = 0; k <= cnt; ++k )
        {
            const real a = ( 2.0 * M_PI * k ) / cnt;
            const real c = std::cos(a);
            const real s = std::sin(a);
            const real inner = spc.innerRadius(a, z);
            pts[i++] = vertex_normal(inner * c, inner * s, z, 0, 0, 1);
            pts[i++] = vertex_normal(outer * c, outer * s, z, 0, 0, 1);
        }
        gym::unmapBufferV3N3();
        gym::drawTriangleStrip(0, i);
        gym::cleanupVN();
    }

    static void drawBottomCap(SpaceRoughAnnulus const& spc, real z, real outer, size_t cnt)
    {
        size_t i = 0;
        flute6* pts = gym::mapBufferV3N3(2 * (cnt + 1));
        for ( size_t k = 0; k <= cnt; ++k )
        {
            const real a = ( 2.0 * M_PI * k ) / cnt;
            const real c = std::cos(a);
            const real s = std::sin(a);
            const real inner = spc.innerRadius(a, z);
            pts[i++] = vertex_normal(outer * c, outer * s, z, 0, 0, -1);
            pts[i++] = vertex_normal(inner * c, inner * s, z, 0, 0, -1);
        }
        gym::unmapBufferV3N3();
        gym::drawTriangleStrip(0, i);
        gym::cleanupVN();
    }
}


void SpaceRoughAnnulus::draw2D(float w) const
{
    const size_t cnt = 128;
    drawSmoothRingAtZ(0.0, r_out_, cnt, w);
    drawRoughRingAtZ(*this, 0.0, cnt, w);
}


void SpaceRoughAnnulus::draw3D() const
{
#if (DIM >= 3)
    const size_t theta_cnt = 160;
    const size_t z_cnt = 32;

    gym::color_both(0.62f, 0.66f, 0.70f, 0.32f);
    drawOuterWall(z_bot_, z_top_, r_out_, theta_cnt);
    drawInnerWall(*this, z_bot_, z_top_, theta_cnt, z_cnt);
    drawTopCap(*this, z_top_, r_out_, theta_cnt);
    drawBottomCap(*this, z_bot_, r_out_, theta_cnt);
#else
    draw2D(1.0f);
#endif
}
#else
void SpaceRoughAnnulus::draw2D(float) const {}
void SpaceRoughAnnulus::draw3D() const {}
#endif
