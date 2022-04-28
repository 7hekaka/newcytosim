// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#ifndef SPACE_PERIODIC_H
#define SPACE_PERIODIC_H

#include "space.h"
#include "modulo.h"

/// a rectangular Space with periodic boundary conditions
/**
 Space `periodic` implements periodic boundary condition in all dimensions.
 The volume has no edge and wraps on itself.

 Parameters:
     - length = total extent in X, Y and Z
     .
 
 To display a periodic Space, use simul:display parameter 'tile'.
 @ingroup SpaceGroup
 */
class SpacePeriodic : public Space
{
    
    /// half the lenth in each dimension
    real half_[4];

    /// Object to handle periodic boundary conditions
    Modulo modulo_;

public:
    
    /// creator
    SpacePeriodic(SpaceProp const*);

    /// change dimensions
    void resize(Glossary& opt);
    
    /// return Modulo Object
    Modulo const* getModulo() const { return &modulo_; }
    
    /// match sizes of Modulo object
    void update();

    /// return bounding box in `inf` and `sup`
    void boundaries(Vector& inf, Vector& sup) const;
    
    /// the volume inside
    real volume() const;
    
    /// true if the point is inside the Space
    bool inside(Vector const&) const;
    
    /// return point on the edge that is closest to `pos`
    Vector project(Vector const& pos) const;
    
    /// return a position inside, resulting from bouncing off on the edges of the Space
    Vector bounce(Vector const&) const;
    
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

