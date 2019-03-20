// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef SPACE_BANANA_H
#define SPACE_BANANA_H

#include "space.h"

/// a bent cylinder of constant diameter terminated by hemispheric caps
/**
 Space `banana` is comprised from a section of a torus,
 terminated by two hemispheres. It is defined by three parameters:
 
    banana length width radius
 
 With:
 - `length` = the overall length minus 2*width
 - `width`  = the diameter of the torus in its crosssections.
 - `radius` = the main radius of the torus, which defines curvature
 .
 
 This class was first conceived by Dietrich Foethke, to simulate S. pombe.
 
 */
class SpaceBanana : public Space
{
private:
    
    /// dimensions
    real  bLength;
    real  bWidth, bWidthSqr;
    real  bRadius;

    /// angle covered by torus section
    real bAngle;
    
    /// X and Y coordinates of the right end
    real bEnd[2];

    /// coordinates of the center of the torus
    Vector bCenter;
    
    /// project on the backbone circle
    Vector project0(Vector const& pos) const;

public:
        
    /// constructor
    SpaceBanana(const SpaceProp* p);
        
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
