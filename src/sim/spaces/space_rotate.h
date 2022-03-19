// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_ROTATE_H
#define SPACE_ROTATE_H

#include "space.h"

/// Rotate a Space 90 degrees around the Y axis
/**
 Space `rotate` rotates another Space.
 This simply rotates any vector argument before calling the Space functions,
 and eventually inversely transforms the result.
 
 @todo SpaceRotate is unfinished: setConfinement() not implemented
 */
class SpaceRotate : public Space
{
private:
        
    ///Space in which objects are excluded
    Space const* mSpace;
    
    /// forward transform of a vector
    Vector forward(Vector const&) const;
    
    /// backward transform of a vector
    Vector backward(Vector const&) const;
    
public:
        
    ///creator
    SpaceRotate(Space *);
    
    ///destructor
    ~SpaceRotate();
    
    /// return bounding box in `inf` and `sup`
    void boundaries(Vector& inf, Vector& sup) const;
    
    /// volume is unchanged
    real volume() const { return mSpace->volume(); }
    
    /// true if the point is inside the Space
    bool inside(Vector const&) const;
    
    /// true if the bead of radius `rad` is inside the Space
    bool allInside(Vector const&, real rad) const;
    
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
