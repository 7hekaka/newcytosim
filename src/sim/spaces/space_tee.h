// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_TEE_H
#define SPACE_TEE_H

#include "space.h"

/// a Capsule with a cylindrical arm perpendicular to it
/**
 Space `tee` is a capsule with a cylindrical arms in a perpendicular direction.
 
    tee length radius arm_position arm_length 

 With:
 - length = half length of central cylinder
 - radius = radius/width of central and arm cylinders, radius of caps
 - arm_position = position of perpendicular arm on the cylinder
 - arm_length = length of perpendicular arm
 .
 

 This is an OLD Space: 
 The re-entrant corners at the base of the arm are not
 properly considered in setInteraction().
 
 @todo Update SpaceTee::setInteraction() if you want to use this SpaceTee.
*/
class SpaceTee : public Space
{
private:
    
    ///the half length of the central cylinder
    real tLength;
    
    ///the length of the perpendicular part on the cylinder
    real tArmLength;
    
    ///the position of the perpendicular part on the cylinder
    real tJunction;
    
    ///the radius of the caps, and square of it
    real tWidth,  tWidthSq;
    
    ///project on base cylinder, return distance
    real projectOnBase(const Vector w, Vector& p)  const;
    
    ///project on side-arm, return distance
    real projectOnArm(const Vector w, Vector& p)   const;
    
    ///project on intersection line
    void projectOnInter(const Vector w, Vector& p) const;
    
public:
        
    ///constructor
    SpaceTee(const SpaceProp*);
   
    /// check number and validity of specified lengths
    void        resize();
    
    /// return bounding box in `inf` and `sup`
    void        boundaries(Vector& inf, Vector& sup) const;
    
    /// the volume inside
    real        volume() const;
    
    /// true if the point is inside the Space
    bool        inside(Vector const&) const;
    
    /// set `proj` as the point on the edge that is closest to `point`
    Vector      project(Vector const& pos) const;
    
    /// OpenGL display function; returns true if successful
    bool        draw() const;
    
};

#endif
