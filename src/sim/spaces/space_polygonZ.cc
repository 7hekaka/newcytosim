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
        poly_.translate(vec.XX, vec.YY);

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

    if ( poly_.complete(REAL_EPSILON) )
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
    size_t hit;
    poly_.project(R, W.z(), P, Z, hit);
    
    real n = P / R;
    return Vector( W.XX * n, W.y() * n, Z);
}


/**
 The current procedure tests the vertices of fibers against the segments of the polygon.
 This fails for non-convex polygon since the re-entrant corners can intersect the fibers.
 
 @todo Also project re-entrant polygon corners on the segments of the Fiber.
 */
void SpacePolygonZ::setInteraction(Vector const& pos, Mecapoint const& pe, Meca& meca, real stiff) const
{
    //Space::setInteraction(pos, pe, meca, stiff); return;
#if ( DIM > 2 )
    real P, R = pos.normXY();
    size_t hit;
    
    Vector prj;

    int edg = poly_.project(R, pos.ZZ, P, prj.ZZ, hit);
    real nX = -poly_.pts_[hit].dy;
    real nY =  poly_.pts_[hit].dx;

    prj.XX = pos.XX * P / R;
    prj.YY = pos.YY * P / R;
    
    if ( edg )
    {
        Vector dir( nX * pos.XX / R, nX * pos.YY / R, nY );
        meca.addPlaneClamp(pe, prj, dir, stiff);
    }
    else
    {
        Vector dir( -pos.YY / R, pos.XX / R, 0 );
        meca.addLineClamp(pe, prj, dir, stiff);
    }
#endif
}
                       

void SpacePolygonZ::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca& meca, real stiff) const
{
    //setInteraction(pos, pe, meca, stiff);
    std::cerr << "unfinished SpacePolygonZ::setInteractions(with radius)\n";
}


void SpacePolygonZ::setInteractions(Meca& meca) const
{
    /// @todo add interactions between fibers and reentrant corners!
#if ( 0 )
    real stiffness = 0;
    Vector dir(0,0,1);
    
    for (Fiber * fib=fibers.first(); fib; fib=fib->next() )
    {
        for ( size_t s = 0; s < fib->nbSegments() ; ++s )
        {
            //project on point on segment
            if ( 0 <= abs  &&  abs < 1 )
                ;// meca.addCylinderClampX(Interpolation(seg, abs), Vector(0,0,0), neck, 100)
        }
    }
#endif
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"
#include "gle_flute.h"

void SpacePolygonZ::draw3D() const
{
    const size_t npts = poly_.nbPoints();
    Polygon::Point2D const* pts = poly_.pts_;
    
    for ( size_t n = 0; n < npts; n++ )
    {
        // do not display special edges
        if ( pts[n].spot )
            continue;
        
        GLfloat R1(pts[n].xx), R2(pts[n+1].xx);
        GLfloat Z1(pts[n].yy), Z2(pts[n+1].yy);
        GLfloat nX(pts[n].dy), nY(-pts[n].dx);
        
        if (( R1 >= 0 ) & ( R2 >= 0 ))
        {
            fluteVN * flu = gle::mapVertexNormalBuffer(2+2*gle::pi_twice);
            size_t i = 0;
            for ( size_t j = 0; j <= gle::pi_twice; ++j )
            {
                float C = gle::cos_(j), S = gle::sin_(j);
                flu[i++] = fluteVN{R2*C, R2*S, Z2, nX*C, nX*S, nY};
                flu[i++] = fluteVN{R1*C, R1*S, Z1, nX*C, nX*S, nY};
            }
            gle::unmapVertexNormalBuffer();
            glDrawArrays(GL_TRIANGLE_STRIP, 0, i);
        }
    }
    glDisableClientState(GL_NORMAL_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}


//display rings around:
void SpacePolygonZ::drawRings() const
{
    const size_t npts = poly_.nbPoints();
    Polygon::Point2D const* pts = poly_.pts_;
    
    glLineWidth(0.5);
    for ( size_t n = 0; n < npts; n++ )
    {
        GLfloat R = pts[n].xx;
        GLfloat Z = pts[n].yy;
        GLfloat nX(pts[n].dy), nY(-pts[n].dx);
        if ( R > 0 )
        {
            fluteVN * flu = gle::mapVertexNormalBuffer(2+gle::pi_twice);
            size_t i = 0;
            for ( size_t j = 0; j <= gle::pi_twice; ++j )
            {
                float C = gle::cos_(j), S = gle::sin_(j);
                flu[i++] = fluteVN{R*C, R*S, Z, nX*C, nX*S, nY};
            }
            gle::unmapVertexNormalBuffer();
            glDrawArrays(GL_LINE_LOOP, 0, i);
        }
    }
}


#else

void SpacePolygonZ::draw3D() const
{
}

#endif
