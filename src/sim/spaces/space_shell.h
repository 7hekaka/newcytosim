// Cytosim -- Spherical shell space
#ifndef SPACE_SHELL_H
#define SPACE_SHELL_H

#include "space.h"

/// spherical shell centered at the origin.
/**
 Space `spherical_shell` confines objects between an inner and an outer sphere.

 Parameters:
   - outer = outer radius
   - inner = inner radius
   - radius_outer, outer_radius = alternative names for outer
   - radius_inner, inner_radius = alternative names for inner
   - thickness or width = shell thickness when one radius is already specified
 */
class SpaceShell : public Space
{
private:
    real r_in_;
    real r_out_;

    /// apply a radial clamp to target radius
    void clampToRadius(Vector const&, Mecapoint const&, Meca&, real target, real stiff) const;

public:
    /// constructor
    SpaceShell(SpaceProp const*);

    /// change dimensions
    void resize(Glossary& opt);

    /// return bounding box in `inf` and `sup`
    void boundaries(Vector& inf, Vector& sup) const;

    /// outer diameter
    real thickness() const { return 2*r_out_; }

    /// volume inside the shell
    real volume() const;

    /// area of inner plus outer surfaces
    real surface() const;

    /// true if the point is inside the shell
    bool inside(Vector const&) const;

    /// true if a bead of radius `rad` is entirely inside
    bool allInside(Vector const&, real rad) const;

    /// a random position inside the shell volume
    Vector place() const;

    /// direct normal direction calculation
    Vector normalToEdge(Vector const&) const;

    /// direct surface placement
    Vector placeOnEdge(real) const;

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
    void draw2D(float) const override;

    /// OpenGL display function
    void draw3D() const override;
};

#endif
