// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Copyright 2019, Cambridge University
#include "dim.h"
#include "space_mesh.h"
#include "mecapoint.h"
#include "exceptions.h"
#include "glossary.h"
#include "meca.h"
#include <fstream>


SpaceMesh::SpaceMesh(SpaceProp const* p, Glossary& opt)
: Space(p)
{
    if ( DIM < 3 )
        throw InvalidParameter("mesh is only usable in 3D");
    
#if COMPLETE_SPACE_MESH
    mMesh.read(p->shape_spec);

    real x;
    if ( opt.set(x, "scale") )
        mMesh.scale(x, x);

    if ( opt.set(x, "inflate") )
        mMesh.inflate(x);
#endif
    resize();
}


SpaceMesh::~SpaceMesh()
{
}


/**
 recalculate bounding box, volume
 and points offsets that are used to project
 */
void SpaceMesh::resize()
{
#if COMPLETE_SPACE_MESH
    mVolume = mMesh.volume();
    if ( mVolume < 0 )
        throw InvalidParameter("mesh volume is < 0");

    real box[4];
    mMesh.find_extremes(box);
    mInf.set(box[0], box[2], 0);
    mSup.set(box[1], box[3], 0);
#endif
}


bool SpaceMesh::inside(Vector const& w) const
{
#if COMPLETE_SPACE_MESH
    return mMesh.inside(w.XX, w.YY, w.ZZ);
#endif
    return true;
}


Vector SpaceMesh::project(Vector const& w) const
{
#if COMPLETE_SPACE_MESH
    return mMesh.project(w);
#endif
    return Vector(0,0,0);
}


/**
 The current procedure tests the vertices of fibers against the segments of the polygon.
 This fails for non-convext polygon since the re-entrant corners can intersect the fibers.
 
 @todo Also project re-entrant polygon corners on the segments of the Fiber.
 */
void SpaceMesh::setInteraction(Vector const& pos, Mecapoint const& pe, Meca & meca, real stiff) const
{    
      std::cerr << "unfinished SpaceMesh::setInteraction\n";
}


void SpaceMesh::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca & meca, real stiff) const
{
    std::cerr << "unfinished SpaceMesh::setInteraction(with radius)\n";
}


void SpaceMesh::setInteractions(Meca & meca, FiberSet const& fibers) const
{
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"

bool SpaceMesh::draw() const
{
    return true;
}

#else

bool SpaceMesh::draw() const
{
    return false;
}

#endif
