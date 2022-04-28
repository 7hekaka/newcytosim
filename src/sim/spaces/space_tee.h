// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_TEE_H
#define SPACE_TEE_H

#include "space.h"

/// a Capsule with a cylindrical arm perpendicular to it
/**
 Space `tee` is a capsule with a cylindrical arms in a perpendicular direction.

 Parameters:
     - length = half length of central cylinder
     - width = radius/width of central and arm cylinders, radius of caps
     - junction = position of perpendicular arm on the cylinder
     - arm_length = length of perpendicular arm
     .
 
 This is an OLD Space: 
 The re-entrant corners at the base of the arm are not
 properly considered in setConfinement().
 
 @todo Update SpaceTee::setConfinement() if you want to use this SpaceTee.
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
    real tRadius,  tRadiusSq;
    
    /// update derived lengths
    void update() { tRadiusSq = square(tRadius); }
    
    ///project on base cylinder, return distance
    real projectOnBase(const Vector&, Vector& p)  const;
    
    ///project on side-arm, return distance
    real projectOnArm(const Vector&, Vector& p)   const;
    
    ///project on intersection line
    void projectOnInter(const Vector&, Vector& p) const;
    
public:
        
    ///constructor
    SpaceTee(SpaceProp const*);
   
    /// change dimensions
    void resize(Glossary& opt);

    /// return bounding box in `inf` and `sup`
    void boundaries(Vector& inf, Vector& sup) const;
    
    /// the volume inside
    real volume() const;
    
    /// true if the point is inside the Space
    bool inside(Vector const&) const;
    
    /// return point on the edge that is closest to `pos`
    Vector project(Vector const& pos) const;

    
    /// write to file
    void write(Outputter&) const;

    /// get dimensions from array `len`
    void setLengths(const real len[8]);
    
    /// read from file
    void read(Inputter&, Simul&, ObjectTag);

    
    /// OpenGL display function
    void draw2D(float) const;
    
    /// OpenGL display function
    void draw3D() const;
};

#endif
