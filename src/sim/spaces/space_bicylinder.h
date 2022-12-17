// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.
#ifndef SPACE_BICYLINDER_H
#define SPACE_BICYLINDER_H

#include "space.h"

///the intersection of two cylinders
/**
 Space `bicylinder` is defined as the intersection of two cylinders with
 the same radius: one with axis X, and the other with axis Y.
 Seen from the top Z direction, it appears square,
 but seen from the X or Y direction, crossections are circular.
 https://en.wikipedia.org/wiki/Steinmetz_solid
 
 Parameters:
     - radius = radius of the cylinders
     .

 @ingroup SpaceGroup
 */
class SpaceBicylinder : public Space
{    
    /// apply a force directed towards the edge of the Space
    static void setConfinement(Vector const& pos, Mecapoint const&, Meca&, real stiff, real rad);

private:
    
    /// the radius of the cylinder
    real radius_;
    
public:
        
    ///creator
    SpaceBicylinder(SpaceProp const*);

    /// change dimensions
    void resize(Glossary& opt);
 
    /// return bounding box in `inf` and `sup`
    void boundaries(Vector& inf, Vector& sup) const;
    
    /// radius
    real thickness() const { return 2*radius_; }

    /// the volume inside
    real volume() const { return 16 / 3.0 * cube(radius_); }

    /// the area of the edge surface
    real surface() const { return 16 * square(radius_); }

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
};

#endif

