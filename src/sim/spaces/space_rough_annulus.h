// Cytosim -- annular slab with a static rough inner wall
#ifndef SPACE_ROUGH_ANNULUS_H
#define SPACE_ROUGH_ANNULUS_H

#include "space.h"
#include <cstddef>

/// annulus with a rigid, spatially varying inner radius
/**
 Parameters:
   - outer = outer radius
   - inner = mean inner radius
   - amplitude / rough_amplitude = maximum inner-wall radial variation
   - theta_mode = maximum angular roughness mode
   - z_mode = maximum axial roughness mode (3D only)
   - rough_seed = deterministic seed for the roughness field
   - rough_components = number of Fourier components in the roughness field
   - phase = optional phase offset for deterministic variants
   - bottom, top = z-extent (only used if DIM >= 3)

 The outer wall remains smooth. The inner wall is rigid but imperfect:
 its local radius varies as a seeded irregular Fourier field in theta and z.
*/
class SpaceRoughAnnulus : public Space
{
private:
    real r_in_ = 0;
    real r_out_ = 0;
#if (DIM >= 3)
    real z_bot_ = 0;
    real z_top_ = 0;
#endif
    real amplitude_ = 0;
    real theta_mode_ = 11;
    real z_mode_ = 7;
    real phase_ = 0;
    unsigned rough_seed_ = 17;
    unsigned rough_components_ = 12;

    real averageInnerRadiusSqr(size_t theta_cnt, size_t z_cnt) const;
    void validate() const;

public:
    /// local inner-wall radius at angular coordinate `theta` and height `z`
    real innerRadius(real theta, real z) const;

    /// creator
    SpaceRoughAnnulus(SpaceProp const*);

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

    /// closest point on the boundary
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
