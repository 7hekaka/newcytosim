// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#ifndef EVENT_H
#define EVENT_H

#include "object.h"

class Meca;
class Simul;
class Glossary;

/// an Event acts on the simulation world by executing code
/**
 An Event is a class that can perform some action in the simulation world,
 specified as a line of code interpreted by Cytosim.
 This can be used to add or remove objects or change parameter values.
 
 The firing time can be specified to occur:
     - only once at a given time by setting the parameter `time`,
     - at regular intervals by setting `delay`,
     - at stochastic time by setting `rate`.
 .
 
 It is a special class in the sense that is not associated with a Property,
 and can be created with 'new' without a 'set' beforehand.
 
 Events are not saved to trajectory files.
*/
class Event: public Object
{
    /// clear member variables
    void clear();
    
public:
    
    /**
     @defgroup EventPar Parameters of Event
     @ingroup Parameters
     These are the parameters for an Event
     @{
     */

    /// code to be executed
    std::string activity;
    
    /// true if event will fire multiple times
    bool recurrent;

    /// true if event can fire multiple time within the same time step
    bool multiplexed;
    
    /// rate of occurence of firing events
    double rate;
    
    /// delay in unit time between firing events (used if `rate` is not set)
    double delay;
    
    /// time after which event will not fire
    double stop;
    
    ///@}
    
    /// time of next event
    double nextTime;
    
public:

    /// default constructor
    Event() { clear(); }
    
    /// constructor
    Event(double time, Glossary&);

    /// destructor
    virtual ~Event();
    
    /// set next firing time
    void fire_always_after(double time);

    /// set next firing time
    void fire_once_at(double time);

    /// recalculate next firing time, given current time
    void reload(double time);
    
    /// a unique character identifying the class
    static const ObjectTag TAG = 'q';

    /// an ASCII character identifying the class of this object
    ObjectTag tag() const { return TAG; }

    /// returns 0, since Event have no Property
    Property const* property() const { return nullptr; }

    //--------------------------------------------------------------------------
    
    /// Stochastic simulation
    void step(Simul&);
    
    /// add interactions to a Meca
    void setInteractions(Meca&) const {}
    
    
    //--------------------------------------------------------------------------

    /// a static_cast<> of Object::next()
    Event * next() const { return static_cast<Event*>(next_); }
    
    /// a static_cast<> of Object::prev()
    Event * prev() const { return static_cast<Event*>(prev_); }
    
    //--------------------------------------------------------------------------

    /// read
    void read(Inputter&, Simul&, ObjectTag);
    
    /// write
    void write(Outputter&) const;
};


#endif
