// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#ifndef SPACE_CYLINDERP_H
#define SPACE_CYLINDERP_H

#include "space.h"
#include "modulo.h"

/// A cylinder of axis X with periodic boundary conditions in X
/**
 Space `cylinderP` is a cylinder with periodic boundary conditions
 along the X-axis. It has no ends and loops on itself like a torus,
 but without the curvature.

 Parameters:
     - length = length of the cylinder in X
     - radius = radius of the cylinder
     .

 To display a periodic Space, use simul:display parameter 'tile'.
 @ingroup SpaceGroup
 */
class SpaceCylinderP : public Space
{
private:
    
    /// half the length of the central cylinder
    real half_;
    
    /// the radius of the cylinder
    real radius_;

    /// Object to handle periodic boundary conditions
    Modulo modulo_;
    
public:
        
    ///creator
    SpaceCylinderP(SpaceProp const*);

    /// change dimensions
    void resize(Glossary& opt);
 
    /// return Modulo Object
    Modulo const* getModulo() const { return &modulo_; }
    
    /// match sizes of Modulo object
    void update();

    /// return bounding box in `inf` and `sup`
    void boundaries(Vector& inf, Vector& sup) const;
    
    /// radius
    real thickness() const { return 2*radius_; }
    
    /// the volume inside
    real volume() const { return 2 * M_PI * half_ * square(radius_); }

    /// surface area of the boundary
    real surface() const { return 2 * M_PI * ( 2 * half_ + radius_ ) * radius_; }

    /// true if the point is inside the Space
    bool inside(Vector const&) const;

    /// true if the bead of radius `rad` is inside the Space
    bool allInside(Vector const&, real rad) const;
    
    /// a random position inside the volume
    Vector place() const;
    
    /// direct normal direction calculation
    Vector normalToEdge(Vector const&) const;
    
    /// direct surface placement
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
    void draw3D() const;
};

#endif

