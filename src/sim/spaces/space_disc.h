// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef SPACE_DISC_H
#define SPACE_DISC_H

#include "space.h"
#include "space_prop.h"

/// A disc centered at the origin at adjustable Z position
/**
 Space `disc` is a disc centered around the origin and parallel to the XY plane.
 It is only usable in 3D, where it offers forces only on the bottom plate.
 
 Parameters:
     - radius = radius of the disc
     - bottom = lower Z
     - top = upper Z

 @ingroup SpaceGroup
 
 FJN, Cambridge 5.10.2021
 */

class SpaceDisc : public Space
{
private:
    
    /// the radius of the disc
    real radius_;
    
    /// position in Z of the bottom limit
    real bot_;
    
    /// position in Z of the top limit
    real top_;

    /// derived quantity: mid_ = ( top_ + bot_ ) / 2
    real mid_;
    
public:
    
    /// constructor
    SpaceDisc(SpaceProp const*);

    /// change dimensions
    void resize(Glossary& opt);
 
    /// the volume inside
    real volume() const;
    
    /// true if the point is inside the Space
    bool inside(Vector const&) const;
    
    /// true if the bead of radius `rad` is inside the Space
    bool allInside(Vector const&, real rad) const;

    /// return bounding box in `inf` and `sup`
    void boundaries(Vector& inf, Vector& sup) const;

    /// a random position in volume
    Vector place() const;
    
    /// direct normal direction calculation
    Vector normalToEdge(Vector const& pos) const { return Vector(0, 0, -1); }
    
    /// a random position on the bottom disc
    Vector placeOnEdge(real) const;

    /// return point on the edge that is closest to `pos`
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

