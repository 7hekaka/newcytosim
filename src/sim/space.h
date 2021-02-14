// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_H
#define SPACE_H

#include <string>
#include <vector>

#include "sim.h"
#include "real.h"
#include "vector.h"
#include "object.h"
#include "common.h"
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
    
    /// read numbers from file into array `len` of size `n_len`
    static void readShape(Inputter&, size_t n_len, real* len, std::string const&);

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
    virtual real   volume() const { return 1; }

    /// return the bounds for the coordinates of the points inside the Space
    /**
     set inf as [ min(X), min(Y), min(Z) ]
     and sup as [ max(X), max(Y), max(Z) ]
     for any point (X, Y, Z) contained inside the Space.
     
     It thus defines a cuboid aligned with the main axes, and containing the entire volume.
     */
    virtual void   boundaries(Vector& inf, Vector& sup) const { inf.set(-1,-1,-1); sup.set(1,1,1); }
    
    /// true if `point` is inside or on the edge of this Space
    virtual bool   inside(Vector const&) const { return true; }
    
    /// set `proj` as the point on the edge that is closest to `point`
    /*
     If the edge is a smooth surface, this should correspond to the usual orthogonal projection.
     */
    virtual Vector project(Vector const& pos) const { ABORT_NOW("base Space is unbounded"); };
    
    /// apply a force directed towards the edge of this Space, for a point located at `pos`
    virtual void   setInteraction(Vector const& pos, Mecapoint const&, Meca&, real stiff) const;
    
    /// apply a force directed towards the edge of this Space deflated by `radius`
    virtual void   setInteraction(Vector const& pos, Mecapoint const&, real rad, Meca&, real stiff) const;
    
#if ( 0 )
    /// apply a force directed towards the edge of this Space
    virtual void   setInteraction(Vector const&, Interpolation const&, Meca&, real stiff) const;

    /// apply a force directed towards the edge of this Space
    virtual void   setInteraction(Interpolation const&, Meca&, real stiff, Confinement conf) const;
#endif
    
    /// true if all points of the sphere (`center`, `radius`) are inside this Space
    virtual bool   allInside(Vector const&, real rad) const;
    
    /// true if no point of the sphere (`center`, `radius`) is inside this Space
    virtual bool   allOutside(Vector const&, real rad) const;
    
    //--------------- FUNCTIONS THAT CAN BE CALCULATED--------------------------
    
    /// returns the maximum absolute value of any coordinate
    real           max_extension() const;

    /// true if `point` is outside this Space ( defined as !inside(point) )
    bool           outside(Vector const& pos)  const { return ! inside(pos); }
    
    /// project `point` on this Space deflated by `radius`, putting the result in `proj`
    Vector         projectDeflated(Vector const&, real rad) const;
    
    /// estimate Volume using a crude Monte-Carlo method with `cnt` calls to Space::inside()
    real           estimateVolume(size_t cnt) const;
    
    
    /// bring a position back inside, as if it bounced off the edges of the Space
    void           bounceOnEdges(Vector&) const;

    /// bring a position back inside, as if it bounced off the edges of the Space
    /** This is also used for periodic boundary conditions*/
    virtual void   bounce(Vector&) const;

    
    /// the square of the distance to the edge of this Space
    real           distanceToEdgeSqr(Vector const&) const;
    
    /// the distance to the edge, always positive
    real           distanceToEdge(Vector const& pos) const { return std::sqrt(distanceToEdgeSqr(pos)); }
    
    /// the distance to the edge, positive if `point` is outside, and negative if inside
    real           signedDistanceToEdge(Vector const&) const;
    
    /// calculate a random position located inside and at most at distance `rad` from the edge
    Vector         randomPlaceNearEdge(real rad, size_t nb_trials) const;
    
    /// calculate a random position located on the edge
    Vector         randomPlaceOnEdge(real rad, size_t nb_trials) const;

    //------------- DERIVED FUNCTIONS THAT CAN BE OVERWRITTEN ------------------

    /// a random position inside the volume, uniformly distributed in the volume
    virtual Vector randomPlace() const;

    /// a Vector perpendicular to the space edge at `point`, directed towards the outside
    virtual Vector normalToEdge(Vector const& pos) const;
    
    /// a random position located on the edge of the Space, uniformly distributed on the surface
    virtual Vector randomPlaceNearEdge(real rad) const { return randomPlaceNearEdge(rad, 1<<20); }
    
    /// a random position located on the edge of the Space, uniformly distributed on the surface
    virtual Vector randomPlaceOnEdge(real rad) const { return randomPlaceOnEdge(rad, 1<<20); }

    //------------------------------ SIMULATION --------------------------------
    
    /// one Monte-Carlo simulation step
    virtual void   step() {}
    
    /// add interactions to a Meca
    virtual void   setInteractions(Meca&) const {}

    //------------------------------ READ/WRITE --------------------------------
    
    /// a unique character identifying the class
    static const ObjectTag TAG = 'e';
    
    /// return unique character identifying the class
    ObjectTag      tag() const { return TAG; }
    
    /// return associated Property
    Property const* property() const { return prop; }
    
    /// returns the name of the Property
    std::string    name() const { return prop->name(); }
    
    /// a static_cast<> of Object::next()
    Space*         next() const { return static_cast<Space*>(nextO); }
    
    /// a static_cast<> of Object::prev()
    Space*         prev() const { return static_cast<Space*>(prevO); }
    
    /// write shape on 16 characters
    static void    writeShape(Outputter&, std::string const&);
    
    /// write to file
    virtual void   write(Outputter&) const;

    /// read from file
    virtual void   read(Inputter&, Simul&, ObjectTag);
    
    /// get dimensions from array `len`
    virtual void   setLengths(const real len[8]) {}
    
    /// print descriptive quantities to stream
    virtual void   report(std::ostream&) const {}

    //------------------------------ DISPLAY -----------------------------------
    
    /// Default 2D display, tracing the outline of a section of the Volume
    void           drawSection(int dim, real pos, size_t cnt) const;

    /// outline the surface using lines, return true if drawn
    virtual void   draw2D() const {}
    
    /// draw surface of the volume, using triangles, return true if drawn
    virtual void   draw3D() const {}

};

#endif

