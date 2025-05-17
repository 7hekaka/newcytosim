// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_polygonZ.h"
#include "exceptions.h"
#include "mecapoint.h"
#include "glossary.h"
#include "polygon.h"
#include "random.h"
#include "meca.h"
#include <fstream>


SpacePolygonZ::SpacePolygonZ(SpaceProp const* p)
: Space(p)
{
    volume_ = 0;
    
    if ( DIM < 3 )
        throw InvalidParameter("polygonZ is only usable in 3D");
}


void SpacePolygonZ::resize(Glossary& opt)
{
    std::string file;
    opt.set(file, "file");

    poly_.read(file);
    
#if ( DIM > 2 )
    Vector vec;
    if ( opt.set(vec, "translate") )
        poly_.transform(1.0, 1.0, vec.XX, vec.YY);

    real len;
    if ( opt.set(len, "inflate") && len > 0 )
        poly_.inflate(len);
#endif

    update();
}


SpacePolygonZ::~SpacePolygonZ()
{
}


/**
 The volume is estimated with a Monte-Carlo approach,
 but taking into account that the volume is axisymmetric.
 The result should be more precise than Space::estimateVolume()
 */
real SpacePolygonZ::estimateVolumeZ(size_t cnt) const
{
    real box[4];
    poly_.find_extremes(box);
    
    real H = box[3] - box[2];

    real Ri = box[0];
    real Ro = box[1];
    if ( Ri < 0 )
    {
        box[0] = 0;
        Ri = 0;
    }
    real W = Ro - Ri;
    
    real in = 0, out = 0;
    for ( size_t i = 0; i < cnt; ++i )
    {
        real x = box[0] + W * RNG.preal();
        real y = box[2] + H * RNG.preal();
        if ( poly_.inside(x, y, 1) )
            in  += x;
        else
            out += x;
    }
    
    return M_PI * H * ( Ro*Ro - Ri*Ri ) * in / ( in + out );
}


/**
 recalculate bounding box, volume
 and points offsets that are used to project
 */
void SpacePolygonZ::update()
{
    if ( poly_.surface() < 0 )
    {
        //std::clog << "flipping clockwise polygon `" << file << "'" << '\n';
        poly_.flip();
    }

    if ( poly_.complete(1e-6) )
        throw InvalidParameter("unfit polygon: consecutive points may overlap");

    real box[4];
    poly_.find_extremes(box);
    inf_.set(-box[1],-box[1], box[2]);
    sup_.set( box[1], box[1], box[3]);

    volume_ = estimateVolumeZ(1<<17);
}


bool SpacePolygonZ::inside(Vector const& W) const
{
    return poly_.inside(W.normXY(), W[2], 1);
}


Vector SpacePolygonZ::project(Vector const& W) const
{
    real P, Z, R = W.normXY();
    unsigned hit;
    poly_.project(R, W.z(), P, Z, hit);
    
    real n = P / R;
    return Vector( W.XX * n, W.y() * n, Z);
}


/**
 The current procedure tests the vertices of fibers against the segments of the polygon.
 This fails for non-convex polygon since the re-entrant corners can intersect the Fiber.
 
 @todo Also project re-entrant polygon corners on the segments of the Fiber.
 */
void SpacePolygonZ::setConfinement(Vector const& pos, Mecapoint const& mp,
                                   Meca& meca, real stiff) const
{
    //Space::setConfinement(pos, mp, meca, stiff); return;
#if ( DIM > 2 )
    real P, R = pos.normXY();
    unsigned hit;
    
    Vector prj;

    int edg = poly_.project(R, pos.ZZ, P, prj.ZZ, hit);
    real nX = -poly_.pts_[hit].dy;
    real nY =  poly_.pts_[hit].dx;

    prj.XX = pos.XX * P / R;
    prj.YY = pos.YY * P / R;
    
    if ( edg )
    {
        Vector dir( nX * pos.XX / R, nX * pos.YY / R, nY );
        meca.addPlaneClamp(mp, prj, dir, stiff);
    }
    else
    {
        Vector dir( -pos.YY / R, pos.XX / R, 0 );
        meca.addLineClamp(mp, prj, dir, stiff);
    }
#endif
}
                       

void SpacePolygonZ::setConfinement(Vector const& pos, Mecapoint const& mp,
                                   real rad, Meca& meca, real stiff) const
{
    //setConfinement(pos, mp, meca, stiff);
    std::cerr << "unfinished SpacePolygonZ::setConfinement(with radius)\n";
}


void SpacePolygonZ::setInteractions(Meca& meca, Simul const&) const
{
    /// @todo add interactions between fibers and reentrant corners!
#if ( 0 )
    real stiffness = 0;
    Vector dir(0,0,1);
    
    for (Fiber * fib=fibers.first(); fib; fib=fib->next() )
    {
        for ( index_t s = 0; s < fib->nbSegments() ; ++s )
        {
            //project on point on segment
            if ( 0 <= abs  &&  abs < 1 )
                ;// meca.addCylinderClampX(Interpolation(seg, abs), Vector(0,0,0), neck, 100)
        }
    }
#endif
}

//------------------------------------------------------------------------------
#pragma mark - OpenGL display

#ifdef DISPLAY

#include "gle.h"
#include "gym_flute.h"
#include "gym_draw.h"

void SpacePolygonZ::draw3D() const
{
    const unsigned inc = 1;
    const unsigned nbp = poly_.nbPoints();
    Polygon::Point2D const* pts = poly_.pts_;
    
    for ( unsigned n = 0; n < nbp; n++ )
    {
        // do not display special edges
        if ( pts[n].spot )
            continue;
        
        float R1(pts[n].xx), R2(pts[n+1].xx);
        float Z1(pts[n].yy), Z2(pts[n+1].yy);
        float nX(pts[n].dy), nY(-pts[n].dx);
        
        if (( R1 >= 0 ) & ( R2 >= 0 ))
        {
            flute6 * flu = gym::mapBufferV3N3(2+2*gle::pi_4half);
            flute6 * ptr = flu;
            for ( unsigned j = 0; j <= gle::pi_4half; j += inc )
            {
                float C = gle::cos_(j), S = gle::sin_(j);
                ptr[0] = {R2*C, R2*S, Z2, nX*C, nX*S, nY};
                ptr[1] = {R1*C, R1*S, Z1, nX*C, nX*S, nY};
                ptr += 2;
            }
            gym::unmapBufferV3N3();
            gym::drawTriangleStrip(0, ptr-flu);
        }
    }
    gym::cleanupVN();
}


//display rings around:
void SpacePolygonZ::drawRings(float width) const
{
    const unsigned inc = 1;
    const unsigned nbp = poly_.nbPoints();
    Polygon::Point2D const* pts = poly_.pts_;
    
    for ( unsigned n = 0; n < nbp; n++ )
    {
        float R(pts[n].xx);
        float Z(pts[n].yy);
        if ( R > 0 )
        {
            float nX(pts[n].dy), nY(-pts[n].dx);
            flute6 * flu = gym::mapBufferV3N3(2+gle::pi_4half);
            for ( unsigned j = 0; j <= gle::pi_4half; j += inc )
            {
                float C = gle::cos_(j), S = gle::sin_(j);
                flu[j] = {R*C, R*S, Z, nX*C, nX*S, nY};
            }
            gym::unmapBufferV3N3();
            gym::drawLineStrip(width, 0, gle::pi_4half+1);
        }
    }
    gym::cleanupVN();
}


#else

void SpacePolygonZ::draw3D() const
{
}

#endif
