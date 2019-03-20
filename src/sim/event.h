// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef EVENT_H
#define EVENT_H

#include "assert_macro.h"
#include "event_prop.h"
#include "object.h"

class Meca;
class Simul;

/// Performs actions on the simulation
/** 
*/
class Event: public Object
{
    
    friend class EventSet;
    
private:
    
    EventProp const* prop;
    
public:

    /// default constructor
    Event(EventProp const* p) { prop = p; }
    
    /// destructor
    virtual      ~Event();
    
    
    /// a unique character identifying the class
    static const ObjectTag TAG = 't';

    /// an ASCII character identifying the class of this object
    ObjectTag            tag() const { return TAG; }
    
    /// Property associated with the Object
    Property const* property() const { return prop; }

    //--------------------------------------------------------------------------
    
    /// prepare for simulation
    virtual void          prepare() {}
    
    /// monte-carlo simulation step
    virtual void          step() {}
    
    /// add interactions to the Meca
    virtual void          setInteractions(Meca &) const {}
    
    //--------------------------------------------------------------------------
    
    /// a static_cast<> of Node::next()
    Event *   next()  const  { return static_cast<Event*>(nNext); }
    
    /// a static_cast<> of Node::prev()
    Event *   prev()  const  { return static_cast<Event*>(nPrev); }
    
    //--------------------------------------------------------------------------

    /// read
    void          read(Inputter&, Simul&, ObjectTag);
    
    /// write
    void          write(Outputter&) const;
};



#endif
