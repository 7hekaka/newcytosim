// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_SQUARE_H
#define SPACE_SQUARE_H

#include "space.h"

///a rectangular region
/**
 Space `square` is a 2D or 3D rectangular volume.
 
    square sizeX sizeY sizeZ

 With:
 - sizeX = half-width along X
 - sizeY = half-width along Y
 - sizeZ = half-width along Z
 .


 @ingroup SpaceGroup
 */
class SpaceSquare : public Space
{
private:
    
    /// apply a force directed towards the edge of the Space
    static void setInteraction(const real pos[], Mecapoint const&, Meca &, real stiff, const real dim[]);
    
public:
    
    ///creator
    SpaceSquare(const SpaceProp*);
    
    /// check number and validity of specified lengths
    void        resize() { Space::checkLengths(DIM, false); }
    
    /// return bounding box in `inf` and `sup`
    void        boundaries(Vector& inf, Vector& sup) const;
    
    /// the volume inside
    real        volume() const;
    
    /// true if the point is inside the Space
    bool        inside(Vector const&) const;

    /// true if a sphere (center, radius) fits in the space, edges included
    bool        allInside(Vector const&, real rad) const;
    
    /// set `proj` as the point on the edge that is closest to `point`
    Vector      project(Vector const& pos) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, Meca &, real stiff) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, real rad, Meca &, real stiff) const;

    
    /// OpenGL display function; returns true if successful
    bool        draw() const;
    
};

#endif


