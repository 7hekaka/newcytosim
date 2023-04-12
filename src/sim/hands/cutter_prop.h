// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University.

#ifndef CUTTER_PROP_H
#define CUTTER_PROP_H

#include "hand_prop.h"

/// values for the 'cutter:selective'
enum { CUT_ANY_FIBER = 0, CUT_IF_BRIDGE = 1, CUT_TOP_FIBER = 2 };


/// Additional Property for Cutter
/**
 @ingroup Properties
 */
class CutterProp : public HandProp
{
    friend class Cutter;
    
public:
    
    /**
     @defgroup CutterPar Parameters of Cutter
     @ingroup Parameters
     Inherits @ref HandPar.
     @{
     */
    
    /// different modes of operation
    /**
     Cutting can be restricted to certain configurations:
     selective = 0: cutting occurs always on any type of fiber
     selective = 1: cut only if the cutter is part of a doubly bound couple
     selective = 2: cut only the topmost fiber (highest Z) if doubly bound
     */
    int selective;
    
    /// unidimensional diffusion coefficient while bound to a Fiber
    real line_diffusion;

    /// rate of cutting event
    real cutting_rate;
    
    /// max distance from the minus end for cutting
    real cutting_range;
    
    /// dynamic state of newly created Fiber ends
    /**
      This defines the dynamic state of the new ends that are created by a cut:
      - new_end_state[0] is for the new plus end,
      - new_end_state[1] is for the new minus end
      .
     */
    state_t new_end_state[2];
    
    /// @}
    
    real line_diffusion_dt, movability_dt;
    
private:
    
    real cutting_rate_dt;
    
public:
    
    /// constructor
    CutterProp(const std::string& n) : HandProp(n)  { clear(); }
    
    /// destructor
    ~CutterProp() { }
    
    /// return a Hand with this property
    virtual Hand * newHand(HandMonitor*) const;
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// compute values derived from the parameters
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new CutterProp(*this); }

    /// write all values
    void write_values(std::ostream&) const;
    
};

#endif

