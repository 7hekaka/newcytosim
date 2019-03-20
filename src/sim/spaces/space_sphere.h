// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_SPHERE_H
#define SPACE_SPHERE_H

#include "space.h"

/// sphere centered at the origin.
/**
 Space `sphere` is a sphere centered around the origin
 
    sphere radius
 
With:
 - radius = radius of the sphere
 .
 
 @ingroup SpaceGroup
 */

class SpaceSphere : public Space
{
protected:
    
    /// the radius of the sphere (set by mLength[0])
    real & radius;
    
    /// square of the radius (set by mLengthSqr[0])
    real & radiusSqr;
    
public:
    
    /// constructor
    SpaceSphere(const SpaceProp*);

    /// check number and validity of specified lengths
    void        resize() { Space::checkLengths(1, false); }

    /// return bounding box in `inf` and `sup`
    void        boundaries(Vector& inf, Vector& sup) const;
    
    /// the volume inside
    real        volume() const;
    
    /// true if the point is inside the Space
    bool        inside(Vector const&) const;
    
    /// a random position inside the volume
    Vector      randomPlace() const { return Vector::randB(radius); }
    
    /// direct normal direction calculation
    Vector      normalToEdge(Vector const& pos) const { return normalize(pos); }

    /// set `proj` as the point on the edge that is closest to `point`
    Vector      project(Vector const& pos) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, Meca &, real stiff) const;

    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, real rad, Meca &, real stiff) const;
    
    /// OpenGL display function; returns true if successful
    bool        draw() const;
    
};

#endif

