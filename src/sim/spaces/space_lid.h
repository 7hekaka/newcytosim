// Cytosim 3.0 -  Copyright Francois Nedelec et al.  EMBL 2007-2013

#ifndef SPACE_LID_H
#define SPACE_LID_H

#include "space.h"
#include "modulo.h"

/**
 a rectangular Space with partial periodic boundary conditions similar to SpaceStrip
 The space can grow with a constant velocity or a velocity depending on total force on the membrane
 
    strip sizeX sizeY sizeZ
 
 With:
 - sizeX = half-width along X
 - sizeY = half-width along Y
 - sizeZ = half-width along Z
 .
 
 Author: Antonio Z. Politi
 */
class SpaceLid : public Space
{
private:
    
    /// the position of the top lid (an alias to mLength[3])
    real &  ceiling;
    
    /// force during last time step
    mutable real force;
    
public:
    
    /// creator
    SpaceLid(const SpaceProp*);
    
    /// check number and validity of specified lengths
    void       resize();
    
    /// initialize Modulo Object
    void       setModulo(Modulo&) const;

    /// true if the Space is periodic in dimension ii
    bool       isPeriodic(int ii) const { return ( ii < DIM-1 ); }
    
    /// maximum extension along each axis
    Vector     extension()        const;
    
    /// the volume inside
    real       volume()           const;
    
    /// true if the point is inside the Space
    bool       inside(Vector const&) const;
    
    /// true if a sphere (\a center, \a radius) is entirely insde this Space
    bool       allInside(Vector const&, real rad) const;
    
    /// true if a sphere (\a center, \a radius) is entirely outside this Space
    bool       allOutside(Vector const&, real rad) const;

    /// project point on the closest edge of the Space
    Vector     project(Vector const& pos) const;
    
    
    /// apply a force directed towards the edge of the Space
    void       setInteraction(Vector const& pos, Mecapoint const&, Meca &, real stiff) const;
    
    /// apply a force directed towards the edge of the Space
    void       setInteraction(Vector const& pos, Mecapoint const&, real rad, Meca &, real stiff) const;

    
    /// add interactions to a Meca
    void       setInteractions(Meca &, FiberSet const&) const;
    
    /// one Monte-Carlo simulation step
    void       step();
    
    /// near the top edge
    Vector     randomPlaceNearEdge(real radius, unsigned long) const;

    /// OpenGL display function; returns true if successful
    bool       draw() const;
    
};

#endif

