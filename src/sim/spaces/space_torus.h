// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef SPACE_TORUS_H
#define SPACE_TORUS_H

#include "space.h"

///a torus of constant diameter centered on the origin
/**
 Space `torus` is defined by two parameters: 
 
    torus radius width
 
 With:
 - `radius` = the main radius of the torus
 - `width`  = the diameter of the torus in its crosssections.
 .
 
 */
class SpaceTorus : public Space
{
private:
    
    /// main radius
    real  bRadius;
    
    /// thickness
    real  bWidth, bWidthSqr;
    
    /// project on the backbone
    Vector project0(Vector const& pos) const;
    
public:
        
    /// constructor
    SpaceTorus(const SpaceProp* p);
        
    /// this is called if any length has been changed
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
