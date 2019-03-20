// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_CYLINDER_H
#define SPACE_CYLINDER_H

#include "space.h"

///a cylinder of axis X
/**
 Space `cylinder` is radial symmetric along the X-axis.
 The crosssection in the YZ plane is a disc.
 It is terminated by flat discs at `X = +/- length/2`.
 For spherical caps, see `capsule`.
 
    cylinder length radius

 With:
 - length = half-length of the cylinder along X
 - radius = radius of the cylinder
 .

 @ingroup SpaceGroup
 */
class SpaceCylinder : public Space
{    
    /// apply a force directed towards the edge of the Space
    static void setInteraction(Vector const& pos, Mecapoint const&, Meca &, real stiff, real len, real rad);

private:
    
    /// half the length of the central cylinder (alias to mLength[0])
    real &      length;
    
    /// the radius of the cylinder (alias to mLength[1])
    real &      radius;
    
public:
        
    ///creator
    SpaceCylinder(const SpaceProp*);
    
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

