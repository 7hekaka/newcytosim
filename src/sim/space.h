// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
#ifndef SPACE_H
#define SPACE_H

#include <string>
#include <vector>

#include "cymdef.h"
#include "real.h"
#include "vector.h"
#include "object.h"
#include "modulo.h"
#include "space_prop.h"


class Mecapoint;
class Interpolation;
class FiberSet;
class Modulo;
class Simul;
class Meca;

//------------------------------------------------------------------------------

/// Defines the spatial constrains in cytosim
/**
The Space defines a few important functions:\n
 - volume(), which returns the volume contained inside the boudaries,
 - inside(x), which tells if a position `x` is inside the space or not,
 - project(x,p), which calculates `p`, the closest point to `x` on the edge of the space.
 .
The edges are considered to be inside.
*/
class Space : public Object
{
protected:
    
    /// write shape on 16 characters
    static void writeShape(Outputter&, std::string const&);

    /// read shape string from file
    static void readShape(Inputter&, std::string const& expected);
    
    /// read lengths from file into array `len` of size `n_len`
    static size_t readLengths(Inputter&, size_t n_len, real* len);
    
    /// read shape string and lengths
    static size_t readShape(Inputter&, size_t n_len, real* len, std::string const&);
    
    /// bounce Z between boundaries at B and B+W, returning image within [B, B+W]
    static real bounce1(real Z, real const& B, real W)
    {
        Z = ( Z - B ) / W;
        int i = (int)std::floor(Z);
        W = std::copysign(W, (i&1)?-1:1);
        return B + W * ( Z - ((i+1)&~1) );
    }

public:
    
    /// parameters
    SpaceProp const* prop;
    
    /// constructor
    Space(SpaceProp const*);
    
    /// destructor
    virtual ~Space();
    
    //------------------------------ BASIC -------------------------------------
    
    /// this is called if any length has been changed
    virtual void resize(Glossary& opt) {};

    /// initialize Modulo if this Space has some periodic dimensions
    virtual Modulo const* getModulo() const { return nullptr; }
    
    /// radius used for piston effect (and defined only for certain shapes)
    virtual real thickness() const { return 0; }

    //------------------------------ OBJECT ------------------------------------
    
    /// the volume inside in 3D, or the surface area in 2D
    virtual real volume() const { return -1; }
    
    /// the surface area of the boundary in 3D
    virtual real surface() const { return -1; }

    /// return the bounds for the coordinates of the points inside the Space
    /**
     By setting values for its two arguments:
         `inf` as [ min(X), min(Y), min(Z) ]
         `sup` as [ max(X), max(Y), max(Z) ]
     for any point (X, Y, Z) contained inside the Space, the function will define
     a cuboid aligned with the main axes, containing the entire volume.
     */
    virtual void boundaries(Vector& inf, Vector& sup) const { inf.set(-1,-1,-1); sup.set(1,1,1); }
    
    /// true if `point` is inside or on the edge of this Space
    virtual bool inside(Vector const&) const { return true; }
    
    /// return point on the edge that is closest to `pos`
    /*
     If the edge is a smooth surface, this should correspond to the usual orthogonal projection.
     */
    virtual Vector project(Vector const& pos) const { ABORT_NOW("base Space is unbounded"); };
    
    /// apply a force directed towards the edge of this Space, for a point located at `pos`
    virtual void setConfinement(Vector const& pos, Mecapoint const&, Meca&, real stiff) const;
    
    /// apply a force directed towards the edge of this Space deflated by `radius`
    virtual void setConfinement(Vector const& pos, Mecapoint const&, real rad, Meca&, real stiff) const;
    
#if ( 0 )
    /// apply a force directed towards the edge of this Space
    virtual void setConfinement(Vector const&, Interpolation const&, Meca&, real stiff) const;

    /// apply a force directed towards the edge of this Space
    virtual void setConfinement(Interpolation const&, Meca&, real stiff, Confinement conf) const;
#endif
    
    /// true if all points of the sphere (`center`, `radius`) are inside this Space
    virtual bool allInside(Vector const&, real rad) const;
    
    /// true if no point of the sphere (`center`, `radius`) is inside this Space
    virtual bool allOutside(Vector const&, real rad) const;
    
    //--------------- FUNCTIONS THAT CAN BE CALCULATED--------------------------
    
    /// returns the maximum absolute value of any coordinate
    real maxExtension() const;

    /// true if `point` is outside this Space ( defined as !inside(point) )
    bool outside(Vector const& pos)  const { return ! inside(pos); }
    
    /// return projection of `point` on edges deflated by `rad`
    Vector projectDeflated(Vector const&, real rad) const;
    
    /// estimate Volume using a crude Monte-Carlo method with `cnt` calls to Space::inside()
    real estimateVolume(size_t cnt) const;
    
    /// bring a position back inside, as if it bounced off the edges of the Space
    Vector bounceOnEdges(Vector const&) const;

    /// return a position inside, resulting from bouncing off on the edges of the Space
    /** This is also used for periodic boundary conditions*/
    virtual Vector bounce(Vector const&) const;

    
    /// the square of the distance to the edge of this Space
    real distanceToEdgeSqr(Vector const&) const;
    
    /// the distance to the edge, always positive
    real distanceToEdge(Vector const& pos) const { return std::sqrt(distanceToEdgeSqr(pos)); }
    
    /// the distance to the edge, positive if `point` is outside, and negative if inside
    real signedDistanceToEdge(Vector const&) const;
    
    /// return a random position located inside and at most at distance `rad` from the edge
    Vector placeNearEdge(real rad, size_t max_trials) const;
    
    /// a crude method returning a random position located on the surface
    Vector onSurface(real rad, size_t max_trials) const;

    //------------- DERIVED FUNCTIONS THAT CAN BE OVERWRITTEN ------------------

    /// a random position inside the volume, uniformly distributed in the volume
    virtual Vector place() const;

    /// a Vector perpendicular to the space edge at `pos`, directed outward
    virtual Vector normalToEdge(Vector const& pos) const;
    
    /// a random position located inside and within distance `rad` from the surface
    virtual Vector placeNearEdge(real rad) const { return placeNearEdge(rad, 1<<20); }
    
    /// a random position located on the surface of the Space, broadly distributed
    virtual Vector placeOnEdge(real rad) const { return onSurface(rad, 1<<20); }

    //------------------------------ SIMULATION --------------------------------
    
    /// perform one Monte-Carlo simulation step
    virtual void step() {}
    
    /// add interactions to a Meca
    virtual void setInteractions(Meca&, Simul const&) const {}

    //------------------------------ READ/WRITE --------------------------------
    
    /// a unique character identifying the class
    static const ObjectTag TAG = 'e';
    
    /// return unique character identifying the class
    ObjectTag tag() const { return TAG; }
    
    /// return associated Property
    Property const * property() const { return prop; }
    
    //--------------------------------------------------------------------------

    /// a static_cast<> of Object::next()
    Space * next() const { return static_cast<Space*>(next_); }
    
    /// a static_cast<> of Object::prev()
    Space * prev() const { return static_cast<Space*>(prev_); }
    
    //--------------------------------------------------------------------------

    /// write dimensions to file
    virtual void write(Outputter&) const;

    /// read dimensions from file
    virtual void read(Inputter&, Simul&, ObjectTag);
    
    /// get dimensions from array `len`
    virtual void setLengths(const real len[8]) {}

    /// convert pointer to Space* if possible
    static Space* toSpace(Object * obj)
    {
        if ( obj  &&  obj->tag() == TAG )
            return static_cast<Space*>(obj);
        return nullptr;
    }

    //------------------------------ DISPLAY -----------------------------------
    
    /// Default 2D display, tracing the outline of a section of the Volume
    void drawSection(int dim, real pos, size_t cnt, float width) const;

    /// outline the surface using lines, return true if drawn
    virtual void draw2D(float) const {}
    
    /// draw surface of the volume, using triangles, return true if drawn
    virtual void draw3D() const {}

};

#endif

