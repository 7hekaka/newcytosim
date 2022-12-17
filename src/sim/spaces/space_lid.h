// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#ifndef SPACE_LID_H
#define SPACE_LID_H

#include "space.h"
#include "modulo.h"
#include "space_dynamic_prop.h"


///a rectangular Space with partial periodic boundary conditions, and mobile boundary
/**
 SpaceLid is a rectangular Space with partial periodic boundary conditions in
 all except the last dimension, like SpaceStrip.
 In addition, the top surface can move, depending on the force it experiences.
 
 Parameters:
 
     - length = extent in X, and Y in 3D
     - bottom = lower limit in Z
     - top    = position of mobile upper limit in Z

 Author: Antonio Z. Politi (2013)

 @ingroup SpaceGroup
*/
class SpaceLid : public Space
{
private:
    
    /// half the lenth in each dimension
    real half_[2];
    
    /// lower position of the bottom limit: Y in 2D and Z in 3D
    real bot_;
    
    /// upper position of the top limit: Y in 2D and Z in 3D
    real top_;
    
    /// Object to handle periodic boundary conditions
    Modulo modulo_;
    
    /// force applied on top boundary in last time step
    mutable real force_;

public:
    
    /// creator
    SpaceLid(SpaceDynamicProp const*);
    
    /// Property
    SpaceDynamicProp const* prop() const { return static_cast<SpaceDynamicProp const*>(Space::prop); }

    /// change dimensions
    void resize(Glossary& opt);
    
    /// match sizes of Modulo object
    void update();

    /// return Modulo Object
    Modulo const* getModulo() const { return &modulo_; }
    
    /// return bounding box in `inf` and `sup`
    void boundaries(Vector& inf, Vector& sup) const;
    
    /// near the top edge
    Vector placeOnEdge(real) const;

    /// the volume inside
    real volume() const;
    
    /// the surface area of the edge
    real surface() const;

    /// true if the point is inside the Space
    bool inside(Vector const&) const;
    
    /// true if a sphere (\a center, \a radius) is entirely inside this Space
    bool allInside(Vector const&, real rad) const;
    
    /// true if a sphere (\a center, \a radius) is entirely outside this Space
    bool allOutside(Vector const&, real rad) const;

    /// project point on the closest edge of the Space
    Vector project(Vector const& pos) const;
    
    /// return a position inside, resulting from bouncing off on the edges of the Space
    Vector bounce(Vector const&) const;

    
    /// apply a force directed towards the edge of the Space
    void setConfinement(Vector const& pos, Mecapoint const&, Meca&, real stiff) const;
    
    /// apply a force directed towards the edge of the Space
    void setConfinement(Vector const& pos, Mecapoint const&, real rad, Meca&, real stiff) const;
    
    /// add interactions to a Meca
    void setInteractions(Meca&, Simul const&) const;
    
    /// one Monte-Carlo simulation step
    void step();
    
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

