// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.
#ifndef SPACE_ELLIPSE_H
#define SPACE_ELLIPSE_H

#include "dim.h"
#include "space.h"

/// enables code to shortcut the projection in case two of the ellipse axes are equal
#define ELLIPSE_HAS_SPHEROID 0

/// ellipse in 2D, ellipsoid or spheroid in 3D 
/**
 The ellipse/ellipsoid is aligned with the principal axes X, Y and Z.
 
 Parameters:
     - length = total length in X, Y and Z directions
     .
 
 The projection of one point on the surface of the ellipse is done numerically. 
 In 3D this is the only solution if the 3 axes have different length.
 setConfinement() relies on project() and thus uses the tangent plane at the
 projection point to approximate the confinement force.

 @ingroup SpaceGroup
 */

class SpaceEllipse : public Space
{
protected:

    /// dimensions
    real radius_[3];
  
    /// dimensions squared
    real radiusSqr_[3];

#if ELLIPSE_HAS_SPHEROID
    /// indicates that two axes are equal
    int spheroid_;
#endif
    
    /// amount added to radius in DISPLAY
    real thickness_;
    
    /// update derived lengths
    void update();
    
public:
        
    /// creator
    SpaceEllipse(SpaceProp const*);
        
    /// change dimensions
    void resize(Glossary& opt);
 
    /// return bounding box in `inf` and `sup`
    void boundaries(Vector& inf, Vector& sup) const;
    
    /// direct normal direction calculation
    Vector normalToEdge(Vector const& pos) const;
    
    /// the volume inside
    real volume() const;
    
    /// the area of the edge surface
    real surface() const;

    /// true if the point is inside the Space
    bool inside(Vector const&) const;
    
    /// return point on the edge that is closest to `pos`
    Vector1 project1D(Vector1 const&) const;
    
    /// return point on the edge that is closest to `pos`
    Vector2 project2D(Vector2 const&) const;
    
    /// return point on the edge that is closest to `pos`
    Vector3 project3D(Vector3 const&) const;
    
    /// return point on the edge that is closest to `pos`
    Vector project(Vector const& pos) const
    {
#if ( DIM == 1 )
        return project1D(pos);
#elif ( DIM == 2 )
        return project2D(pos);
#else
        return project3D(pos);
#endif
    }
    

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

