// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_PERIODIC_H
#define SPACE_PERIODIC_H

#include "space.h"
#include "modulo.h"

/// a rectangular Space with periodic boundary conditions
/**
 Space `periodic` implements periodic boundary condition in all dimensions.
 The volume has no edge and wraps on itself.
 
    periodic sizeX sizeY sizeZ
 
 With:
 - sizeX = half-width along X
 - sizeY = half-width along Y
 - sizeZ = half-width along Z
 .
 
 */
class SpacePeriodic : public Space
{
    
public:
    
    /// creator
    SpacePeriodic(const SpaceProp*);

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
    
    /// OpenGL display function; returns true if successful
    bool       draw() const;
    
};

#endif

