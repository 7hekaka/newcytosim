// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_POLYGON_H
#define SPACE_POLYGON_H

#include "space.h"
#include "polygon.h"

/// a polygonal convex region in space
/**
 Space `polygon` implements a polygon. It works best for convex polygon.
 In 3D, the thickness in Z can be specified to describe a generalized 
 cylinder of axis Z, that has the 2D polygon as cross-section.
 
 Parameters:
     - file: name of file with polygon data
     - height: height of polygon in Z
    .

 Alternatively:
     - order : number of sides
     - radius : distance from center
     - angle : rotation offset in radian
     .
 
 Example:
 
     change cell { order=4; radius=13; angle=0.7853; }
 
 @ingroup SpaceGroup
 @todo add SpacePolygon::setConfinement() for re-entrant corners
*/
class SpacePolygon : public Space
{
private:
        
    /// pre-calculated bounding box derived from poly_
    Vector inf_, sup_;

    /// The 2D polygon object
    Polygon poly_;
    
    /// Surface of polygon
    real surface_;
    
    /// half the total height in Z
    real height_;

    /// update derived lengths
    void update();

public:
        
    /// constructor
    SpacePolygon(const SpaceProp *);
    
    /// destructor
    ~SpacePolygon();
    
    /// change dimensions
    void resize(Glossary& opt);
    
    /// return bounding box in `inf` and `sup`
    void boundaries(Vector& inf, Vector& sup) const { inf=inf_; sup=sup_; }
    
    /// the volume inside
    real volume() const { return ( DIM>2 ? 2*height_ : 1 ) * surface_; }
    
    /// true if the point is inside the Space
    bool inside(Vector const&) const;
    
    /// a random position inside the volume
    Vector place() const;

    /// return point on the edge that is closest to `pos`
    Vector project(Vector const& pos) const;

    /// apply a force directed towards the edge of the Space
    void setConfinement(Vector const& pos, Mecapoint const&, Meca&, real stiff) const;
    
    /// apply a force directed towards the edge of the Space
    void setConfinement(Vector const& pos, Mecapoint const&, real rad, Meca&, real stiff) const;
    
    /// add interactions between fibers and reentrant corners
    void setInteractions(Meca&, Simul const&) const;
    
    /// write to file
    void write(Outputter&) const;

    /// get dimensions from array `len`
    void setLengths(const real len[8]);

    /// read from file
    void read(Inputter&, Simul&, ObjectTag);

    /// OpenGL display function
    void drawPolygon(float, float) const;
    /// OpenGL display function
    void draw2D(float width) const { drawPolygon(width, width); }
    /// OpenGL display function
    void draw3D() const;
    /// draw polygon points
    void drawPolygonPoints() const;
};

#endif

