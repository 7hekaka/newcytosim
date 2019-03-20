// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_CYLINDERZ_H
#define SPACE_CYLINDERZ_H

#include "space.h"

///a cylinder of axis Z
/**
 Space `cylinderZ` is radial symmetric along the Z-axis.
 The crosssection in the XY plane is a disc.

    cylinderZ radius bottom top
 
 With:
 - radius = radius of cylinder
 - bottom = smallest Z
 - top = highest Z
 .

 @ingroup SpaceGroup
 */
class SpaceCylinderZ : public Space
{    
    /// apply a force directed towards the edge of the Space
    static void setInteraction(Vector const& pos, Mecapoint const&, Meca &, real stiff, real, real, real);

private:
    
    /// the radius of the cylinder (alias to mLength[0])
    real &      radius;
    
    /// position in Z of the bottom limit (alias to mLength[1])
    real &      bot;
    
    /// position in Z of the top limit (alias to mLength[2])
    real &      top;
    
public:
        
    ///creator
    SpaceCylinderZ(const SpaceProp*);
    
    /// check number and validity of specified lengths
    void        resize();
    
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

