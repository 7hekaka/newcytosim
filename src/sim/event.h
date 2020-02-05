// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef EVENT_H
#define EVENT_H

#include "object.h"

class Meca;
class Simul;
class Glossary;

/// an Event performs action on the simulation by executing code
/**
 An Event is a class that can perform some action in the simulation world,
 specified as a line of code interpreted by Cytosim.
 This can be used for example to add or remove objects.
 
 The firing time can be specified to occur:
     - only once at a given time with `time`,
     - at regular interval  with `delay`,
     - at stochastic time with `rate`.
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

    /// rate of occurence of firing events
    real        rate;
    
    /// delay in unit time between firing events (used if `rate` is not set)
    real        delay;

    /// true if event will fire multiple times
    bool        recurrent;
    
    ///@}
    
    /// time of next event
    real        nextTime;
    
public:

    /// default constructor
    Event() { clear(); }
    
    /// constructor
    Event(real time, Glossary&);

    /// destructor
    virtual ~Event();

    /// set next firing time
    void fire_once(real time);

    /// recalculate next firing time, given current time
    void reset(real time);
    
    /// a unique character identifying the class
    static const ObjectTag TAG = 'q';

    /// an ASCII character identifying the class of this object
    ObjectTag       tag() const { return TAG; }

    /// returns 0, since Event have no Property
    Property const* property() const { return nullptr; }

    //--------------------------------------------------------------------------
    
    /// monte-carlo simulation step
    void      step(Simul&);
    
    /// add interactions to a Meca
    void      setInteractions(Meca&) const {}
    
    
    /// a static_cast<> of Node::next()
    Event *   next()  const  { return static_cast<Event*>(nNext); }
    
    /// a static_cast<> of Node::prev()
    Event *   prev()  const  { return static_cast<Event*>(nPrev); }
    

    /// read
    void      read(Inputter&, Simul&, ObjectTag);
    
    /// write
    void      write(Outputter&) const;
};


#endif
