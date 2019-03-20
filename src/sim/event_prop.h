// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef EVENT_PROP_H
#define EVENT_PROP_H

#include "real.h"
#include "property.h"


class Event;


/// Property for an Event
/**
 @ingroup Properties
 
 */
class EventProp : public Property
{
    
    friend class Event;
    
public:
    
    /**
     @defgroup EventPar Parameters of Event
     @ingroup Parameters
     @{
     */
    
    /// code to be executed
    std::string   code;
    
    /// rate at which code is executed
    real          rate;
    
    /// @}
    

public:
    
    /// constructor
    EventProp(const std::string& n) : Property(n)  { clear(); }
    
    /// destructor
    ~EventProp() { }
    
    /// create an Event with this property
    Event * newEvent(Glossary * opt = nullptr) const;

    /// identifies the property
    std::string category() const { return "event"; }
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// check and derive parameters
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new EventProp(*this); }

    /// write all values
    void write_values(std::ostream&) const;
    
};

#endif

