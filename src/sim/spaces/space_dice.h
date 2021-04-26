// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
#ifndef SPACE_DICE_H
#define SPACE_DICE_H

#include "space.h"

#define ADVANCED_DICE_INTERACTIONS 1

/// A rectangle ( or a cube ) with rounded edges. 
/**
 Space `dice` is a cube with smooth edges.

 It is build by expanding a cube by a distance `radius` in all directions.
 Mathematically, a point is inside the `dice` if it is at most at distance
 `radius` from the inner cube obtained by subtracting `radius` to the sizes.
 The dice is thus included in the rectangular space of similar size.

 Parameters:
     - length = total extent along X, Y and Z
     - radius = rounding radius of edges
     .

 Note: Dice::setInteraction() relies on project(), and numerical instabilities
 may arise in particular if `radius << size`, because determining the tangent
 plane to a point becomes imprecise.
 
 @ingroup SpaceGroup
 */
class SpaceDice : public Space
{
private:
    
    /// half the lenth in each dimension
    real half_[4];
    
    /// the radius by which the corners are smoothed
    real edge_;
    
    /// the square of the radius
    real edgeSqr_;
    
    /// calculate edgeSqr_
    void update() { edgeSqr_ = square(edge_); }
    
    /// apply a force directed towards the edge of the Space
    void setInteraction(Vector const& pos, Mecapoint const&, Meca&, real, const real[], real) const;

public:
    
    /// constructor
    SpaceDice(SpaceProp const*);

    /// change dimensions
    void resize(Glossary& opt);
 
    /// return bounding box in `inf` and `sup`
    void boundaries(Vector& inf, Vector& sup) const;
    
    /// the volume inside
    real volume() const;
    
    /// the area of the edge surface
    real surface() const;

    /// true if the point is inside the Space
    bool inside(Vector const&) const;
    
    /// true if the bead is inside the Space
    bool allInside(Vector const&, real rad) const;

    /// set `proj` as the point on the edge that is closest to `point`
    Vector project(Vector const& pos) const;
    
#if ADVANCED_DICE_INTERACTIONS
    /// apply a force directed towards the edge of the Space
    void setInteraction(Vector const& pos, Mecapoint const&, Meca&, real stiff) const;
    
    /// apply a force directed towards the edge of the Space
    void setInteraction(Vector const& pos, Mecapoint const&, real rad, Meca&, real stiff) const;
#endif

    /// write to file
    void write(Outputter&) const;

    /// get dimensions from array `len`
    void setLengths(const real len[8]);
    
    /// read from file
    void read(Inputter&, Simul&, ObjectTag);
    
    /// OpenGL display function
    void draw2D() const;

    /// OpenGL display function
    void draw3D() const;
};

#endif
