// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
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

    dice sizeX sizeY sizeZ radius

 With:
 - sizeX = half-width along X
 - sizeY = half-width along Y
 - sizeZ = half-width along Z
 - radius = rounding radius of edges
 .

 Note: Dice::setInteraction() relies on project(), and numerical instabilities
 may arise in particular if `radius << size`, because determining the tangent
 plane to a point becomes imprecise.
*/
class SpaceDice : public Space
{
public:
    
    /// the radius by which the corners are smoothed (alias to mLength[3])
    real &      radius;
    
    /// the square of the radius
    real &      radiusSqr;
    
public:
        
    /// constructor
    SpaceDice(const SpaceProp*);
    
    /// check number and validity of specified lengths
    void        resize();
    
    /// return bounding box in `inf` and `sup`
    void        boundaries(Vector& inf, Vector& sup) const;
    
    /// the volume inside
    real        volume() const;
    
    /// true if the point is inside the Space
    bool        inside(Vector const&) const;
    
    /// set `proj` as the point on the edge that is closest to `point`
    Vector      project(Vector const& pos) const;

    
    /// OpenGL display function; returns true if successful
    bool        draw() const;
    
};

#endif
