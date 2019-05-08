// Cytosim 3.0 -  Copyright Francois Nedelec et al.  EMBL 2007-2013

#ifndef SPACE_LID_H
#define SPACE_LID_H

#include "space.h"
#include "modulo.h"

/**
 SpaceLid is a rectangular Space with partial periodic boundary conditions in
 all except the last dimension. The top surface can move with a constant velocity,
 or with a velocity proportional to the force it experience.
 
 Parameters:
     - length = along X, Y and Z
     - ceiling = position of the top boundary
     .
 
 Author: Antonio Z. Politi

 @ingroup SpaceGroup
*/
/// Periodic boundary conditions in all but the last dimension
class SpaceLid : public Space
{
private:
    
    /// dimensions
    real  length_[3];
    
    /// the position of the top lid
    real  top_;

    /// force during last time step
    mutable real force_;
    
public:
    
    /// creator
    SpaceLid(SpaceProp const*);
    
    /// update geometry
    void        resize(Glossary& opt);
    
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
    void        setInteractions(Meca &, FiberSet const&) const;
    
    /// one Monte-Carlo simulation step
    void        step();
    
    /// near the top edge
    Vector      randomPlaceNearEdge(real radius, unsigned long) const;

    /// OpenGL display function; returns true if successful
    bool        draw() const;
    
    /// write to file
    void        write(Outputter&) const;

    /// get dimensions from array `len`
    void        setLengths(const real len[8]);
    
    /// read from file
    void        read(Inputter&, Simul&, ObjectTag);

};

#endif

