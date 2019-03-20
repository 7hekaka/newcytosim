// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_MESH_H
#define SPACE_MESH_H

#include "space.h"

/// a volume defined by a triangular mesh
/**
 Space `mesh` implements a generic volume in 3D.
 The mesh is read from a file.

    mesh file_name

 @ingroup SpaceGroup
*/
class SpaceMesh : public Space
{
private:
    
    /// The 3D mesh object
    //Mesh              mMesh;
    
    /// pre-calculated bounding box derived from mMesh
    Vector            mInf, mSup;
    
    /// Volume calculated from mMesh
    real              mVolume;

public:
        
    /// constructor
    SpaceMesh(const SpaceProp *, Glossary&);
    
    /// destructor
    ~SpaceMesh();
    
    /// return bounding box in `inf` and `sup`
    void        boundaries(Vector& inf, Vector& sup) const { inf=mInf; sup=mSup; }
    
    /// the volume inside
    real        volume() const { return mVolume; }
    
    /// true if the point is inside the Space
    bool        inside(Vector const&) const;
    
    /// set `proj` as the point on the edge that is closest to `point`
    Vector      project(Vector const& pos) const;

    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, Meca &, real stiff) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, real rad, Meca &, real stiff) const;
    
    /// add interactions between fibers and reentrant corners
    void        setInteractions(Meca &, FiberSet const&) const;

    /// update since length have changed
    void        resize();
    
    /// OpenGL display function; returns true if successful
    bool        draw() const;
    
};

#endif

