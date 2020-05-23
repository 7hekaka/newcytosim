// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef NUCLEUS_H
#define NUCLEUS_H

#include "nucleus_prop.h"
#include "organizer.h"
#include "sphere.h"
#include "fiber.h"

//------------------------------------------------------------------------------
/// Organizer built around a Sphere
/**
A Nucleus attaches Fibers to a Sphere:\n
 - organized(0) is the Sphere
 - organized(n) for n>0 is a Fiber attached to the sphere
 - prop->stiffness is the stiffness of the link.
 .

 @ingroup OrganizerGroup
 */
class Nucleus : public Organizer
{
public:
        
    /// Sphere on which the Nucleus is built
    Sphere*       nuSphere;

    /// Properties for the Nucleus
    NucleusProp const* prop;
        
    //------------------- construction and destruction -------------------------
    /// constructor
    Nucleus(NucleusProp const* p) : nuSphere(nullptr), prop(p) { }
    
    /// destructor
    virtual      ~Nucleus();

    /// create a Nucleus and requested associated Objects
    ObjectList    build(Glossary&, Simul&);
    
    //------------------- simulation -------------------------------------------    

    /// Stochastic simulation
    void          step();
    
    ///add interactions for this object to a Meca
    void          setInteractions(Meca&) const;

    //------------------- querying the nucleus ---------------------------------    
    
    ///position of center of gravity (returns the center of the sphere)
    Vector        position() const { return nuSphere->position(); }
    
    ///the Sphere on which the nucleus is built
    Sphere *      sphere()   const { return nuSphere; }
    
    /// i-th fiber attached to the nucleus
    Fiber *       fiber(size_t i) const { return static_cast<Fiber*>(organized(i)); }
    
    
    /// retrieve links end-points for display
    bool          getLink(size_t, Vector&, Vector&) const;
    
    /// return Sphere
    Mecable*      core() const { return nuSphere; }
    
    /// display parameters
    PointDisp const* disp() const { if ( nuSphere ) return nuSphere->prop->disp; return nullptr; }

    //------------------------------ read/write --------------------------------
    
    /// a unique character identifying the class
    static const ObjectTag TAG = 'n';
    
    /// return unique character identifying the class
    ObjectTag       tag() const { return TAG; }
    
    /// return associated Property
    Property const* property() const { return prop; }
    
    /// read from IO
    void            read(Inputter&, Simul&, ObjectTag);
    
    /// write to IO
    void            write(Outputter&) const;

};


#endif

