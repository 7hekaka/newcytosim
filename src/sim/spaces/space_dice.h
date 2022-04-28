// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
#ifndef SPACE_DICE_H
#define SPACE_DICE_H

#include "space.h"

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
    static void setConfinement(Vector const& pos, Mecapoint const&, Meca&, real, const real[], real);

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
    
    /// true if the bead of radius `rad` is inside the Space
    bool allInside(Vector const&, real rad) const;

    /// return point on the edge that is closest to `pos`
    Vector project(Vector const& pos) const;
    
    /// apply a force directed towards the edge of the Space
    void setConfinement(Vector const&, Mecapoint const&, Meca&, real stiff) const;
    
    /// apply a force directed towards the edge of the Space
    void setConfinement(Vector const&, Mecapoint const&, real rad, Meca&, real stiff) const;

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
