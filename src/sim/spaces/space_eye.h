// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef SPACE_EYE_H
#define SPACE_EYE_H

#include "space.h"

/// The intersection of two discs, or two spheres (not implemented yet). 
/**
 Space `eye` is the intersection of two discs. 
 The 3D version is not implemented.
 
 @code
    eye width height angle radius_curvature
 @endcode

 With:
 - width = maximum width
 - height = height
 - angle = angle
 - radius_curvature = radius of curvature
 .


 The discontinuties at the intersection of the two circles are not
 properly considered in setInteraction().
 
 @todo Update SpaceEye::setInteraction() if you want to use this SpaceEye.
 @todo Rename SpaceEye SpaceTwinSpheres
*/
class SpaceEye : public Space
{
public:
    
    ///creator
    SpaceEye(const SpaceProp*);
        
    /// check number and validity of specified lengths
    void        resize();
    
    /// maximum extension along each axis
    Vector      extension() const;
   
    /// the volume inside
    real        volume() const;
    
    /// true if the point is inside the Space
    bool        inside(const real point[]) const;
    
    /// project point on the closest edge of the Space
    void        project(const real point[], real proj[]) const;
    
    /// OpenGL display function, return true is display was done
    bool        display() const;
    
};

#endif

