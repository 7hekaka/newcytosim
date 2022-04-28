// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#ifndef SPACE_TORUS_H
#define SPACE_TORUS_H

#include "space.h"

///a torus of constant diameter centered on the origin
/**
 Space `torus` is defined by two parameters: 
 
 Parameters:
     - `curvature` = the (large) radius of the torus centerline
     - `width` = the diameter of the torus in its cross sections.
     .

 @ingroup SpaceGroup
 */
class SpaceTorus : public Space
{
private:
    
    /// main radius of curvature
    real bCurve;
    
    /// thickness
    real bRadius, bRadiusSqr;
    
    /// set bRadiusSqr
    void update() { bRadiusSqr = square(bRadius); }
    
    /// project on the backbone
    Vector backbone(Vector const& pos) const;
    
public:
 
    /// constructor
    SpaceTorus(SpaceProp const*);
 
    /// change dimensions
    void resize(Glossary& opt);

    /// return bounding box in `inf` and `sup`
    void boundaries(Vector& inf, Vector& sup) const;
    
    /// radius
    real thickness() const { return 2*bRadius; }

    /// the volume inside
    real volume() const;
    
    /// true if the point is inside the Space
    bool inside(Vector const&) const;
    
    /// return point on the edge that is closest to `pos`
    Vector project(Vector const& pos) const;
    
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
