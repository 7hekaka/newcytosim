// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef BEAD_H
#define BEAD_H

#include "dim.h"
#include "array.h"
#include "object.h"
#include "mecable.h"
#include "bead_prop.h"

class Meca;
class Single;
class SingleProp;

/// The Bead is constructed from a vertex and a radius
/**
 The Bead is the simplest Mecable.
 It represents a spherical object using: 
 - a position vector,
 - a radius.
 .
 The orientational degrees of freedom are neglected.
 Translation follows Stokes's law.
 A Single can be attached in the center of the Bead.
 
 For more elaborate models, see Sphere and Solid.
*/
class Bead : public Mecable
{
private:
    
    /// radius
    real paRadius;

    /// the total drag coefficient for translation
    real paDrag;
    
public:
    
    /// Property
    BeadProp const* prop;
    
    /// create following specifications
    Bead(BeadProp const*, Vector pos, real rad);

    /// destructor
    virtual ~Bead();
    
    //--------------------------------------------------------------------------
    
    /// Bead only accepts translation
    int mobile() const { return 1; }

    /// return the position in space of the object
    Vector pos() const { return Vector(pPos); }

    /// return the position in space of the object (virtual function)
    Vector position() const { return Vector(pPos); }
    
    /// move the object position ( position += given vector )
    void translate(Vector const& x) { x.add_to(pPos); }
    
    /// set the object position ( position = given vector )
    void setPosition(Vector const& x) { x.store(pPos); }

    //--------------------------------------------------------------------------
        
    /// the radius of the Bead
    real radius() const { return paRadius; }
    
    /// the volume of the bead
    real radiusSqr() const { return paRadius * paRadius; }
    
    /// the volume of the Bead
    real volume() const { return ((M_PI*4.0/3.0) * paRadius) * (paRadius * paRadius); }

    /// set the radius of the Bead
    void resize(real R) { assert_true(R>0); paRadius = R; setDragCoefficient(); }
    
    //--------------------------------------------------------------------------
    
    /// sets the mobility
    void setDragCoefficient();
    
    /// the total drag-coefficient of object (force = drag * speed)
    real dragCoefficient() const { return paDrag; }
    
    /// The mobility of a model vertex ( speed = mobility * point_force )
    real pointMobility() const { return 1 / paDrag; }

    /// the mobility is already set in resize()
    void prepareMecable() { }
    
    /// calculates the speed of points in Y, for the forces given in X
    void projectForces(const real* X, real* Y) const;
    
    /// add contribution of Brownian forces
    real addBrownianForces(real const* rnd, real, real* rhs) const;

    /// add the interactions due to confinement
    void setInteractions(Meca&) const;
    
    //--------------------------------------------------------------------------

    /// a static_cast<> of Object::next()
    Bead * next() const { return static_cast<Bead*>(nextO); }
    
    /// a static_cast<> of Object::prev()
    Bead * prev() const { return static_cast<Bead*>(prevO); }
    
    //--------------------------------------------------------------------------

    /// a unique character identifying the class
    static const ObjectTag TAG = 'b';
    
    /// return unique character identifying the class
    ObjectTag tag() const { return TAG; }
    
    /// return associated Property
    Property const* property() const { return prop; }
    
    /// convert pointer to Fiber* if the conversion seems valid; returns 0 otherwise
    static Bead* toBead(Object * obj)
    {
        if ( obj  &&  obj->tag() == TAG )
            return static_cast<Bead*>(obj);
        return nullptr;
    }
    
    /// convert pointer to Fiber* if the conversion seems valid; returns 0 otherwise
    static Bead const* toBead(Object const* obj)
    {
        if ( obj  &&  obj->tag() == TAG )
            return static_cast<Bead const*>(obj);
        return nullptr;
    }

    //--------------------------------------------------------------------------

    /// read from file
    void read(Inputter&, Simul&, ObjectTag);
    
    /// write to file
    void write(Outputter&) const;

};

#endif
