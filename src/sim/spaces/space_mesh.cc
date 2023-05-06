// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Copyright 2019, Cambridge University
#include "dim.h"
#include "space_mesh.h"
#include "exceptions.h"
#include "mecapoint.h"
#include "glossary.h"
#include "meca.h"
#include <fstream>

#define COMPLETE_SPACE_MESH 0


SpaceMesh::SpaceMesh(SpaceProp const* p)
: Space(p)
{
    if ( DIM < 3 )
        throw InvalidParameter("mesh is only usable in 3D");
}


SpaceMesh::~SpaceMesh()
{
}


/**
 recalculate bounding box, volume
 and points offsets that are used to project
 */
void SpaceMesh::resize(Glossary& opt)
{
#if COMPLETE_SPACE_MESH
    std::string file;
    
    if ( opt.set(file, "file") )
        mesh_.read(file);
    
    real x;
    if ( opt.set(x, "scale") )
        mesh_.scale(x, x);

    if ( opt.set(x, "inflate") )
        mesh_.inflate(x);

    volume_ = mesh_.volume();
    if ( volume_ < 0 )
        throw InvalidParameter("mesh volume is < 0");

    real box[4];
    mesh_.find_extremes(box);
    inf_.set(box[0], box[2], 0);
    sup_.set(box[1], box[3], 0);
#endif
}


bool SpaceMesh::inside(Vector const& W) const
{
#if COMPLETE_SPACE_MESH
    return mesh_.inside(W.XX, W.YY, W.ZZ);
#endif
    return true;
}


Vector SpaceMesh::project(Vector const& W) const
{
#if COMPLETE_SPACE_MESH
    return mesh_.project(W);
#endif
    return Vector(0,0,0);
}


/**
 The current procedure tests the vertices of fibers against the segments of the polygon.
 This fails for non-convext polygon since the re-entrant corners can intersect the Fiber.
 
 @todo Also project re-entrant polygon corners on the segments of the Fiber.
 */
void SpaceMesh::setConfinement(Vector const& pos, Mecapoint const& mp, Meca& meca, real stiff) const
{    
    std::cerr << "unfinished SpaceMesh::setConfinement\n";
}


void SpaceMesh::setConfinement(Vector const& pos, Mecapoint const& mp, real rad, Meca& meca, real stiff) const
{
    std::cerr << "unfinished SpaceMesh::setConfinement(with radius)\n";
}


void SpaceMesh::setInteractions(Meca& meca, Simul const&) const
{
}

//------------------------------------------------------------------------------

void SpaceMesh::write(Outputter& out) const
{
}


void SpaceMesh::read(Inputter& in, Simul&, ObjectTag)
{
}

//------------------------------------------------------------------------------
#pragma mark - OpenGL display

#ifdef DISPLAY
#include "gle.h"

void SpaceMesh::draw3D() const
{
}

#else

void SpaceMesh::draw3D() const {}

#endif
