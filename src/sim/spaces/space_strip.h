// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_STRIP_H
#define SPACE_STRIP_H

#include "space.h"
#include "modulo.h"

///a rectangular Space with partial periodic boundary conditions
/**
 Space `periodic` implements periodic boundary condition in all but the last dimension.
 The volume only has edge in the last dimension, and otherwise wraps on itself.
 The last dimension is Y in 2D and Z in 3D.
 
    strip sizeX sizeY sizeZ
 
 With:
 - sizeX = half-width along X
 - sizeY = half-width along Y
 - sizeZ = half-width along Z
 .
 
 */
class SpaceStrip : public Space
{
public:
    
    /// creator
    SpaceStrip(const SpaceProp*);

    /// check number and validity of specified lengths
    void       resize();
    
    /// initialize Modulo Object
    void       setModulo(Modulo&) const;
    
    /// return bounding box in `inf` and `sup`
    void       boundaries(Vector& inf, Vector& sup) const;
    
    /// the volume inside
    real       volume()           const;
    
    /// true if the point is inside the Space
    bool       inside(Vector const&) const;
    
    /// set `proj` as the point on the edge that is closest to `point`
    Vector     project(Vector const& pos) const;

    
    /// apply a force directed towards the edge of the Space
    void       setInteraction(Vector const& pos, Mecapoint const&, Meca &, real stiff) const;
    
    /// apply a force directed towards the edge of the Space
    void       setInteraction(Vector const& pos, Mecapoint const&, real rad, Meca &, real stiff) const;
    
    
    /// OpenGL display function; returns true if successful
    bool       draw() const;
    
};

#endif

