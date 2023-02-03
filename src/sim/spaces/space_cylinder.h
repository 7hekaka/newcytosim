// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#ifndef SPACE_CYLINDER_H
#define SPACE_CYLINDER_H

#include "space.h"

///a cylinder of axis X
/**
 Space `cylinder` is radial symmetric along the X-axis.
 The cross section in the YZ plane is a disc.
 It is terminated by flat discs at `X = +/- length/2`.
 For spherical caps, see `capsule`.
 
 Parameters:
     - length = length of the cylinder in X
     - radius = radius of the cylinder
     .

 @ingroup SpaceGroup
 */
class SpaceCylinder : public Space
{    
    /// apply a force directed towards the edge of the Space
    static void setConfinement(Vector const& pos, Mecapoint const&, Meca&, real stiff, real len, real rad);

private:
    
    /// half the length of the central cylinder
    real half_;
    
    /// the radius of the cylinder
    real radius_;
    
public:
        
    ///creator
    SpaceCylinder(SpaceProp const*);

    /// change dimensions
    void resize(Glossary& opt);
 
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
    
    /// return point on the edge that is closest to `pos`
    Vector project(Vector const& pos) const;

    /// apply a force directed towards the edge of the Space
    void setConfinement(Vector const& pos, Mecapoint const& mp, Meca& meca, real stiff) const
    {
        setConfinement(pos, mp, meca, stiff, half_, radius_);
    }

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

