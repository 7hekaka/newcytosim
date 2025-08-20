// Cytosim — Annulus / Annular slab
#ifndef SPACE_ANNULUS_H
#define SPACE_ANNULUS_H

#include "space.h"

/// an annulus (2D) or annular slab (3D)
/**
 Parameters:
   - outer = outer radius
   - inner = inner radius  (must satisfy 0 <= inner < outer)
   - radius + width = alternative to set outer/inner
   - radius_outer, radius_inner = alternative names
   - bottom, top = z-extent (only used if DIM >= 3)

 In 2D, volume() returns the area; in 3D, the geometric volume. (See Space doc.) 
*/
class SpaceAnnulus : public Space
{
private:
    real r_in_  = 0;
    real r_out_ = 0;
#if (DIM >= 3)
    real z_bot_ = 0;
    real z_top_ = 0;
#endif

public:
    /// creator
    SpaceAnnulus(SpaceProp const*);

    /// change dimensions
    void resize(Glossary& opt);

    /// return bounding box in `inf` and `sup`
    void boundaries(Vector& inf, Vector& sup) const;

    /// used for piston effect in some contexts (outer diameter)
    real thickness() const { return 2*r_out_; }

    /// 3D: volume; 2D: area
    real volume() const;

    /// 3D: surface area of boundary; 2D: perimeter length
    real surface() const;

    /// true if the point is inside the Space
    bool inside(Vector const&) const;

    /// true if a bead of radius `rad` is entirely inside
    bool allInside(Vector const&, real rad) const;

    /// closest point on the boundary (inner/outer rim and z-caps)
    Vector project(Vector const& pos) const;

    /// write to file
    void write(Outputter&) const;

    /// get dimensions from array `len`
    void setLengths(const real len[8]);

    /// read from file
    void read(Inputter&, Simul&, ObjectTag);

    void draw2D(float w) const override;
    void draw3D() const override;

};

#endif
