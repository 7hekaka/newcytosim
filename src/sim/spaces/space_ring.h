// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#ifndef SPACE_RING_H
#define SPACE_RING_H

#include "space.h"

/// a cylinder of axis X, but without caps
/**
 Space `ring` is radial symmetric around the X-axis.
 It is like a cylinder except that the two end discs are not part of the surface.
 project() will always project on the curvy surface of the cylinder.

 Parameters:
     - length = total extent of the cylinder in X
     - radius = radius of the cylinder
     .

 @ingroup SpaceGroup
 */
class SpaceRing : public Space
{    
private:
    
    /// half the length of the cylinder
    real half_;
    
    /// the radius of the ring
    real radius_;

public:
        
    ///creator
    SpaceRing(SpaceProp const*);
    
    /// change dimensions
    void resize(Glossary& opt);

    /// return bounding box in `inf` and `sup`
    void boundaries(Vector& inf, Vector& sup) const;
    
    /// the volume inside
    real volume() const { return 2 * M_PI * half_ * square(radius_); }

    /// surface area of the boundary
    real surface() const { return 4 * M_PI * half_ * radius_; }

    /// true if the point is inside the Space
    bool inside(Vector const&) const;
    
    /// true if the bead of radius `rad` is inside the Space
    bool allInside(Vector const&, real rad) const;
    
    /// a random position inside the volume
    Vector place() const;

    /// return point on the edge that is closest to `pos`
    Vector project(Vector const& pos) const;
    
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
    void draw3D() const;
    
    /// OpenGL display function
    void draw2D(float) const { draw3D(); }
};

#endif

