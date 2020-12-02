// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_CAPSULE_H
#define SPACE_CAPSULE_H

#include "space.h"

/// a spherocylinder (cylinder capped with hemispheres)
/**
 Space `capsule` is cylinder ending with hemispheres (a spherocylinder)

 Parameters:
     - length: total length in X
     - radius: the radius of the hemisphere and central cylinder
     .

 @ingroup SpaceGroup
 */
class SpaceCapsule : public Space
{    
    /// apply a force directed towards the edge of the Space
    static void setInteraction(Vector const& pos, Mecapoint const&, Meca&, real stiff, real len, real rad);

private:
    
    /// half the distance between the centers of the two hemispherical caps
    real  length_;
    
    /// the radius of the hemispheres caps
    real  radius_;

public:
        
    /// creator
    SpaceCapsule(SpaceProp const*);
    
    /// change dimensions
    void        resize(Glossary& opt);

    /// return bounding box in `inf` and `sup`
    void        boundaries(Vector& inf, Vector& sup) const;
    
    /// radius
    real        thickness() const { return radius_; }

    /// the volume inside
    real        volume() const;
    
    /// the area of the edge surface
    real        surface() const;

    /// true if the point is inside the Space
    bool        inside(Vector const&) const;
    
    /// true if the bead is inside the Space
    bool        allInside(Vector const&, real rad) const;
    
    /// a random position inside the volume
    Vector      randomPlace() const;
    
    /// direct normal direction calculation
    Vector      normalToEdge(Vector const&) const;
    
    /// direct surface placement
    Vector      randomPlaceOnEdge(real) const;

    /// set `proj` as the point on the edge that is closest to `point`
    Vector      project(Vector const& pos) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, Meca&, real stiff) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, real rad, Meca&, real stiff) const;
    
    /// write to file
    void        write(Outputter&) const;

    /// get dimensions from array `len`
    void        setLengths(const real len[8]);
    
    /// read from file
    void        read(Inputter&, Simul&, ObjectTag);
    
    /// OpenGL display function
    void        draw2D() const;
    
    /// OpenGL display function
    void        draw3D() const;
};

#endif
