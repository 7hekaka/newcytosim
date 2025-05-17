// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_polygon.h"
#include "exceptions.h"
#include "mecapoint.h"
#include "glossary.h"
#include "polygon.h"
#include "meca.h"
#include <fstream>


SpacePolygon::SpacePolygon(SpaceProp const* p)
: Space(p)
{
    inf_.reset();
    sup_.reset();
    surface_ = 0;
    height_ = 0;
    
    if ( DIM == 1 )
        throw InvalidParameter("polygon is not usable in 1D");
}


SpacePolygon::~SpacePolygon()
{
}

//------------------------------------------------------------------------------
#pragma mark - I/O

/**
 recalculate bounding box, volume
 and points offsets that are used to project
 */
void SpacePolygon::resize(Glossary& opt)
{
    unsigned ord = 6;
    std::string file;
    
    if ( opt.set(file, "file") )
        poly_.read(file);
    else if ( !prop->dimensions.empty() )
        poly_.read(prop->dimensions);
    else if ( opt.has_key("points") )
    {
        // specify vertices directly:
        unsigned nbp = (unsigned)opt.num_values("points");
        poly_.allocate(nbp);
        for ( unsigned p = 0; p < nbp; ++p )
        {
            Vector2 vec(0,0);
            if ( ! opt.set(vec, "points", p) )
                throw InvalidParameter("polygon:points must be a list of comma-separated points: X Y, X Y, X Y, etc.");
            poly_.setPoint(p, vec.XX, vec.YY);
        }
    }
    else if ( opt.set(ord, "order") )
    {
        real rad = 1, ang = 0;
        opt.set(rad, "radius");
        opt.set(ang, "angle");
        poly_.set(ord, rad, ang);
    }
    else
        return;
    
    real x;
    if ( opt.set(x, "scale") )
        poly_.transform(x, x, 0, 0);

    if ( opt.set(x, "inflate") )
        poly_.inflate(x);
    
#if ( DIM == 3 )
    x = height_;
    if ( opt.set(x, "height") )
        x *= 0.5;
    if ( x < 0 )
        throw InvalidParameter("polygon:height must be >= 0");
    height_ = x;
#endif

    update();
}


void SpacePolygon::update()
{
    surface_ = poly_.surface();
    if ( surface_ < 0 )
    {
        //std::clog << "flipping clockwise polygon `" << file << "'" << '\n';
        poly_.flip();
        surface_ = poly_.surface();
    }
    assert_true( surface_ > 0 );
    
    if ( poly_.complete(1e-6) )
        throw InvalidParameter("unfit polygon: consecutive points may overlap");

    real box[4];
    poly_.find_extremes(box);
    inf_.set(box[0], box[2], -height_);
    sup_.set(box[1], box[3],  height_);
}


bool SpacePolygon::inside(Vector const& W) const
{
#if ( DIM > 2 )
    if ( abs_real(W.ZZ) > height_ )
        return false;
#endif
#if ( DIM > 1 )
    return poly_.inside(W.XX, W.YY, 1);
#else
    return false;
#endif
}


Vector SpacePolygon::place() const
{
    if ( surface_ <= 0 )
        throw InvalidParameter("cannot pick point inside polygon of null surface");
    return Space::place();
}


Vector SpacePolygon::project(Vector const& W) const
{
    Vector P(W);
#if ( DIM == 2 )
    
    unsigned hit;
    poly_.project(W.XX, W.YY, P.XX, P.YY, hit);
    
#elif ( DIM > 2 )
    
    if ( abs_real(W.ZZ) > height_ )
    {
        if ( poly_.inside(W.XX, W.YY, 1) )
        {
            // too high or too low in the Z axis, but inside XY
            P.XX = W.XX;
            P.YY = W.YY;
        }
        else
        {
            // outside in Z and XY
            unsigned hit;
            poly_.project(W.XX, W.YY, P.XX, P.YY, hit);
        }
        P.ZZ = std::copysign(height_, W.ZZ);
    }
    else
    {
        unsigned hit;
        poly_.project(W.XX, W.YY, P.XX, P.YY, hit);
        if ( poly_.inside(W.XX, W.YY, 1) )
        {
            // inside in the Z axis and the XY polygon:
            // to the polygonal edge in XY plane:
            real HH = (W.XX-P.XX)*(W.XX-P.XX) + (W.YY-P.YY)*(W.YY-P.YY);
            // to the top/bottom plates:
            real V = height_ - abs_real(W.ZZ);
            // compare distances
            if ( V * V < HH )
                return Vector(W.XX, W.YY, std::copysign(height_, W.ZZ));
        }
        P.ZZ = W.ZZ;
    }
    
#endif
    return P;
}

//------------------------------------------------------------------------------
#pragma mark - setConfinement

/**
 The current procedure tests the vertices of fibers against the segments of the polygon.
 This fails for non-convext polygon since the re-entrant corners can intersect the Fiber.
 
 @todo Also project re-entrant polygon corners on the segments of the Fiber.
 */
void SpacePolygon::setConfinement(Vector const& pos, Mecapoint const& mp,
                                  Meca& meca, real stiff) const
{    
#if ( DIM > 1 )
    
    unsigned hit;
    real pX, pY;
    int edg = poly_.project(pos.XX, pos.YY, pX, pY, hit);
    real nX = -poly_.pts_[hit].dy;
    real nY =  poly_.pts_[hit].dx;
    
#if ( DIM > 2 )
    bool in = poly_.inside(pos.XX, pos.YY, 1);

    if ( abs_real(pos.ZZ) >= height_ )
    {
        meca.addPlaneClampZ(mp, std::copysign(height_, pos.ZZ), stiff);
        if ( in ) return;
    }
    else if ( in )
    {
        // Compare distance to top/bottom plate:
        real V = height_ - abs_real(pos.ZZ);
        // and distance to polygonal edge in XY plane:
        real HH = (pos.XX-pX)*(pos.XX-pX) + (pos.YY-pY)*(pos.YY-pY);
        
        if ( V * V < HH )
        {
            meca.addPlaneClampZ(mp, std::copysign(height_, pos.ZZ), stiff);
            return;
        }
    }
#endif

    if ( edg )
        meca.addPlaneClamp(mp, Vector(pX,pY,0), Vector(nX,nY,0), stiff);
    else
        meca.addPointClampXY(mp, Vector(pX,pY,0), stiff);
#endif
}


void SpacePolygon::setConfinement(Vector const& pos, Mecapoint const& mp,
                                  real rad, Meca& meca, real stiff) const
{
    //setConfinement(pos, mp, meca, stiff);
    std::cerr << "unfinished SpacePolygon::setConfinement(with radius)\n";
}

#include "fiber_segment.h"
#include "fiber_set.h"

void SpacePolygon::setInteractions(Meca& meca, Simul const&) const
{
#if ( 0 )
    /// WORK IN PROGRESS
    Polygon::Point2D const* pts = poly_.pts_;
    const int n_pik = 2;
    const int inx[n_pik] = { 0, 100 };
    Vector pik[n_pik];
    
    for ( int i = 0; i < n_pik; ++i )
        pik[i].set(pts[inx[i]].xx, pts[inx[i]].yy, 0);
    
    for ( Fiber * fib=fibers.first(); fib; fib=fib->next() )
    {
        real ls = fib->segmentation();
        for ( index_t seg = 0; seg < fib->nbSegments() ; ++seg )
        {
            FiberSegment loc(fib, seg);
            for ( int i = 0; i < n_pik; ++i )
            {
                real dis;
                real abs = loc.projectPoint(pik[i], abs, dis);
                if ( 0 <= abs  &&  abs < ls )
                {
                    if ( !inside(loc.pos(abs)) || !inside(loc.pos1()) || !inside(loc.pos2()) )
                        meca.addPointClamp(Interpolation(loc, abs), pik[i], 100);
                }
            }
        }
    }
#endif
}

//------------------------------------------------------------------------------
#pragma mark - I/O

void SpacePolygon::write(Outputter& out) const
{
    writeMarker(out, TAG);
    writeShape(out, "L");
    out.writeUInt16(2);
    out.writeFloat(height_);
    out.writeFloat(0.f);
}


void SpacePolygon::setLengths(const real len[8])
{
    height_ = len[0];
}


void SpacePolygon::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "L");
    setLengths(len);
    if ( !prop->dimensions.empty() )
    {
        poly_.read(prop->dimensions);
        update();
    }
}


//------------------------------------------------------------------------------
#pragma mark - Display

#ifdef DISPLAY

#include "gle.h"
#include "gym_flute.h"
#include "gym_draw.h"

void SpacePolygon::drawPolygon(float lines, float points) const
{
    const index_t nbp = poly_.nbPoints();
    Polygon::Point2D const* pts = poly_.pts_;
    flute2 * flt = gym::mapBufferV2(nbp+1);
    for ( index_t n = 0; n <= nbp; ++n )
        flt[n].set(pts[n].xx, pts[n].yy);
    gym::unmapBufferV2();
    
    glEnable(GL_STENCIL_TEST);
    glClearStencil(1);
    glClear(GL_STENCIL_BUFFER_BIT);
    glStencilFunc(GL_EQUAL, 1, ~0U);
    glStencilOp(GL_KEEP, GL_ZERO, GL_ZERO);
    if ( lines > 0 )
    {
        gym::drawLineStrip(lines, 0, nbp+1);
    }
    if ( points > 0 )
    {
        gym::drawPoints(points, 0, nbp);
    }
    glClear(GL_STENCIL_BUFFER_BIT);
    glDisable(GL_STENCIL_TEST);
}

/*
void SpacePolygon::drawPolygonPoints() const
{
    const index_t nbp = poly_.nbPoints();
    Polygon::Point2D const* pts = poly_.pts_;
    // indicate index of each point:
    char tmp[32];
    for ( size_t n = 0; n < nbp; ++n )
    {
        snprintf(tmp, sizeof(tmp), "%lu", n);
        gym::translate_ref(pts[n].xx, pts[n].yy, height_);
        fgStrokeString(0, 0, pixelSize, 1, tmp, 1);
    }
}
*/

void SpacePolygon::draw3D() const
{
    const float H(-height_);
    const index_t nbp = poly_.nbPoints();
    Polygon::Point2D const* pts = poly_.data();
    flute3 * flt = gym::mapBufferV3(2*nbp+2);
    for ( index_t i = 0; i <= nbp; ++i )
    {
        float X(pts[i].xx), Y(pts[i].yy);
        flt[2*i  ] = { X, Y, -H };
        flt[2*i+1] = { X, Y,  H };
    }
    gym::unmapBufferV3();
    // display sides
    gym::drawTriangleStrip(0, 2*nbp+2);

    float lines = 2;
    if ( lines )
    {
        // display bottom
        gym::rebindBufferV3(2, 0);
        gym::drawLineStrip(lines, 0, nbp+1);
        // display top
        gym::rebindBufferV3(2, 1);
        gym::drawLineStrip(lines, 0, nbp+1);
    }
    gym::cleanupV();
}

#else

void SpacePolygon::drawPolygon(float, float) const {}
void SpacePolygon::draw3D() const {}

#endif
