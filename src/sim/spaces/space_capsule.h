// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_CAPSULE_H
#define SPACE_CAPSULE_H

#include "space.h"

/// a spherocylinder (cylinder capped with hemispheres)
/**
 Space `capsule` is cylinder ending with hemispheres (a spherocylinder)

    capsule length radius 
 
 With: 
 - length: half the length of the central cylinder
 - radius: the radius of the hemisphere and central cylinder
 .

 @ingroup SpaceGroup
 */
class SpaceCapsule : public Space
{    
    /// apply a force directed towards the edge of the Space
    static void setInteraction(Vector const& pos, Mecapoint const&, Meca &, real stiff, real len, real rad);

private:
    
    /// half the length of the central cylinder (alias to mLength[0])
    real &      length;
    
    /// the radius of the hemisphere (alias to mLength[1])
    real &      radius;
    
    /// the square of the radius (alias to mLengthSqr[1])
    real &      radiusSqr;

public:
        
    /// creator
    SpaceCapsule(const SpaceProp*);
    
    /// check number and validity of specified lengths
    void        resize() { Space::checkLengths(2, false); }

    /// return bounding box in `inf` and `sup`
    void        boundaries(Vector& inf, Vector& sup) const;

    /// the volume inside
    real        volume() const;
    
    /// true if the point is inside the Space
    bool        inside(Vector const&) const;
    
    /// true if the bead is inside the Space
    bool        allInside(Vector const&, real rad) const;
    
    /// a random position inside the volume
    Vector      randomPlace() const;

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
