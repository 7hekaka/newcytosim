// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#ifndef SPACE_STRIP_H
#define SPACE_STRIP_H

#include "space.h"
#include "modulo.h"

///a rectangular Space with partial periodic boundary conditions
/**
 SpaceStrip implements periodic boundary conditions in all dimensions except the
 last one: Y in 2D and Z in 3D. The volume is a cuboid, with edges in the last
 dimension, while other dimensions (X and Y in 3D) wrap periodically.
 
 Parameters:
 
     - length = extent in X, and Y in 3D
     - bottom = lower limit in Z
     - top    = upper limit in Z
 
 To display a periodic Space, use simul:display parameter 'tile'.
 @ingroup SpaceGroup
 */
class SpaceStrip : public Space
{
private:
    
    /// half the lenth in each dimension
    real half_[2];
    
    /// lower position of the bottom limit: Y in 2D and Z in 3D
    real bot_;
    
    /// upper position of the top limit: Y in 2D and Z in 3D
    real top_;
    
    /// derived quantity: mid_ = ( top_ + bot_ ) / 2
    real mid_;
    
    /// an object used for cooking on which a top lid fits
    real pot_;

    /// Object to handle periodic boundary conditions
    Modulo modulo_;

    /// if true, only the bottom plate is used for projection
    int no_top_;
    
public:
    
    /// creator
    SpaceStrip(SpaceProp const*);
    
    /// change dimensions
    void resize(Glossary& opt);
    
    /// match sizes of Modulo object
    void update();

    /// return Modulo Object
    Modulo const* getModulo() const { return &modulo_; }
    
    /// return bounding box in `inf` and `sup`
    void boundaries(Vector& inf, Vector& sup) const;
    
    /// thickness in Z
    real thickness() const { return ( top_ - bot_ ); }

    /// the volume inside
    real volume() const;
    
    /// the surface area of the edge
    real surface() const;

    /// near the top or bottom edge
    Vector placeOnEdge(real) const;
    
    /// a vector aligned along Z
    Vector normalToEdge(Vector const&) const;

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

