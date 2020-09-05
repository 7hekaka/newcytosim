// Cytosim 3.0 -  Copyright Francois Nedelec et al.  EMBL 2007-2013

#ifndef SPACE_LID_H
#define SPACE_LID_H

#include "space.h"
#include "modulo.h"
#include "space_dynamic_prop.h"

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
    
    /// half the length in each dimension
    real  length_[3];

    /// the position of the top lid
    real  top_;

    /// force during last time step
    mutable real force_;
  
    /// Object to handle periodic boundary conditions
    Modulo modulo_;

public:
    
    /// creator
    SpaceLid(SpaceDynamicProp const*);
    
    /// properties
    const SpaceDynamicProp* prop;
    
    /// change dimensions
    void       resize(Glossary& opt);
    
    /// return Modulo Object
    Modulo const* getModulo() const { return &modulo_; }
    
    /// match sizes of Modulo object
    void       update();

    /// true if the Space is periodic in dimension ii
    bool       isPeriodic(int ii) const { return ( ii < DIM-1 ); }
    
    /// return bounding box in `inf` and `sup`
    void       boundaries(Vector& inf, Vector& sup) const;

    /// the volume inside
    real       volume() const;
    
    /// true if the point is inside the Space
    bool       inside(Vector const&) const;
    
    /// true if a sphere (\a center, \a radius) is entirely insde this Space
    bool       allInside(Vector const&, real rad) const;
    
    /// true if a sphere (\a center, \a radius) is entirely outside this Space
    bool       allOutside(Vector const&, real rad) const;

    /// project point on the closest edge of the Space
    Vector     project(Vector const& pos) const;
    
    /// equivalent to 'Modulo::fold'
    void       bounce(Vector&) const;

    
    /// apply a force directed towards the edge of the Space
    void       setInteraction(Vector const& pos, Mecapoint const&, Meca&, real stiff) const;
    
    /// apply a force directed towards the edge of the Space
    void       setInteraction(Vector const& pos, Mecapoint const&, real rad, Meca&, real stiff) const;

    
    /// add interactions to a Meca
    void       setInteractions(Meca&) const;
    
    /// one Monte-Carlo simulation step
    void       step();
    
    /// near the top edge
    Vector     randomPlaceOnEdge(real) const;

    /// OpenGL display function; returns true if successful
    bool       draw() const;
    
    /// write to file
    void       write(Outputter&) const;

    /// get dimensions from array `len`
    void       setLengths(const real len[8]);
    
    /// read from file
    void       read(Inputter&, Simul&, ObjectTag);

};

#endif

