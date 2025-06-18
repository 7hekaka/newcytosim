// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University
#ifndef SPACE_CYLINDERZ_H
#define SPACE_CYLINDERZ_H

#include "space.h"

///a cylinder of axis Z
/**
 Space `cylinderZ` is radial symmetric along the Z axis.
 The cross section in the XY plane is a disc.

 Parameters:
     - radius = radius of cylinder
     - bottom = lower limit in Z
     - top    = upper limit in Z
     - top_edge = smoothing radius at the top edge
     - bottom_edge = smoothing radius at the bottom edge
     .

 Added edge smoothing 09.12.2019
 @ingroup SpaceGroup
 */
class SpaceCylinderZ : public Space
{    
    /// apply a force directed towards the edge of the Space
    static void setConfinement(Vector const& pos, Mecapoint const&, Meca&, real stiff, real, real, real);
    /// apply a force directed towards the edge of the Space
    static void setConfinement(Vector const& pos, Mecapoint const&, Meca&, real stiff, real, real, real, real, real);

private:
    
    /// the radius of the cylinder
    real radius_;
    
    /// position in Z of the bottom limit
    real bot_;
    
    /// position in Z of the top limit
    real top_;
    
    /// the radius of smoothing of the edges
    real top_edge_;
    
    /// the radius of smoothing of the edges
    real bot_edge_;
    
    /// mid height in Z:
    real mid_;
    
    /// calculate derived parameters
    void update() { mid_ = 0.5 * ( top_ + bot_ ); }

public:
        
    ///creator
    SpaceCylinderZ(SpaceProp const*);

    /// change dimensions
    void resize(Glossary& opt);
 
    /// return bounding box in `inf` and `sup`
    void boundaries(Vector& inf, Vector& sup) const;
    
    /// radius
    real thickness() const { return 2*radius_; }
    
    /// direct normal direction calculation
    Vector normalToEdge(Vector const& pos) const;
    
    /// direct surface placement
    Vector placeOnEdge(real) const;

    /// the volume inside
    real volume() const;
    
    /// the area of the edge surface
    real surface() const;

    /// true if the point is inside the Space
    bool inside(Vector const&) const;
    
    /// true if the bead of radius `rad` is inside the Space
    bool allInside(Vector const&, real rad) const;
    
    /// a random position inside the volume
    Vector place() const;
    
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
    void draw3D() const;
};

#endif

