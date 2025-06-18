// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#ifndef MOVABLE_H
#define MOVABLE_H

#include "vector.h"
#include "isometry.h"

class Space;
class Modulo;


/// Can be moved and rotated in space
/**
Movable provides a common interface, for Object that can be moved or rotated.
The actual operations need to be implemented by redefining the virtual functions:
 
    if ( mobile() == 0 ):
        the object has no position defined

    if ( mobile() == 1 ):
        position() and translate() are implemented

    if ( mobile() == 2 ):
        rotate() is implemented

    if ( mobile() == 3 ):
        position() and translate() are implemented
        rotate() is implemented

 To support periodic boundary conditions, foldPosition() should be defined.
 */
class Movable
{
public:
    
    /// constructor
    Movable() {}
    
    /// destructor
    ~Movable() {}
    
    /// true if object can be translated (default=false)
    /**
     mobile() returns a bit field:
     
         ( mobile() & 1 ) indicates if the object can be translated.
         ( mobile() & 2 ) indicates if the object can be rotated.
     
     Thus,
     
         if ( mobile() & 1 ):
             position() and translate() should be implemented
     
         if ( mobile() & 2 ):
             rotate() should be implemented
     */
    virtual int mobile() const { return 0; }
    
    /// return the spatial position of the Object
    virtual Vector position() const { return Vector(0.0,0.0,0.0); }
    
    /// return a vector representing the orientation of the Object
    virtual Vector direction() const { return Vector(1.0,0.0,0.0); }

    /// move Object to specified position
    void setPosition(Vector const& X) { translate( X - position() ); }

    /// move Object ( position += given vector )
    virtual void translate(Vector const&) { ABORT_NOW("undefined Movable::translate() called"); }
    
    /// rotate Object around the Origin
    virtual void rotate(Rotation const&) { ABORT_NOW("undefined Movable::rotate() called"); }

    /// return translation derived from applying rotation around the Origin
    Vector translation(Rotation const& R) const
    {
        Vector pos = position();
        return R * pos - pos;
    }

    /// rotate Object around its current position, using translate() and rotate()
    void revolve(Rotation const& R)
    {
        Vector G = position();
        translate(-G);
        rotate(R);
        translate(G);
    }
    
    /// Apply isometry
    void move(Isometry const& iso)
    {
        switch ( mobile() )
        {
            case 1: translate(iso.mov+translation(iso.rot)); break;
            case 2: rotate(iso.rot); break;
            case 3: rotate(iso.rot); translate(iso.mov); break;
        }
    }

    /// bring object to centered image using periodic boundary conditions
    void foldPosition(Modulo const*) {}

};

#endif

