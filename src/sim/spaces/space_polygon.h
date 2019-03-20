// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_POLYGON_H
#define SPACE_POLYGON_H

#include "space.h"
#include "polygon.h"

/// a polygonal convex region in space
/**
 Space `polygon` implements a polygon. It works best for convex polygon.
 In 3D, and additional HEIGHT can be specified to describe a generalized 
 cylinder of axis Z, that has the 2D polygon as cross-section.
 
 The coordinates of the polygon are read from a file.

    polygon file_name HEIGHT

 @ingroup SpaceGroup
 @todo add SpacePolygon::setInteraction() for re-entrant corners
*/
class SpacePolygon : public Space
{
private:
    
    /// The 2D polygon object
    Polygon           mPoly;
        
    /// pre-calculated bounding box derived from mPoly
    Vector            mInf, mSup;
    
    /// Volume calculated from mPoly
    real              mVolume;
    
    /// half the total height in Z (set by mLength[0])
    real              mHeight;

public:
        
    /// constructor
    SpacePolygon(const SpaceProp *, Glossary&);
    
    /// destructor
    ~SpacePolygon();
    
    /// return bounding box in `inf` and `sup`
    void        boundaries(Vector& inf, Vector& sup) const { inf=mInf; sup=mSup; }
    
    /// the volume inside
    real        volume() const { return mVolume; }
    
    /// true if the point is inside the Space
    bool        inside(Vector const&) const;
    
    /// set `proj` as the point on the edge that is closest to `point`
    Vector      project(Vector const& pos) const;

    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, Meca &, real stiff) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, real rad, Meca &, real stiff) const;
    
    /// add interactions between fibers and reentrant corners
    void        setInteractions(Meca &, FiberSet const&) const;

    /// update since length have changed
    void        resize();
    
    /// OpenGL display function; returns true if successful
    bool        draw() const;
    
};

#endif

