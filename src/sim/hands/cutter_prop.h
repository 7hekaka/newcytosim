// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University.

#ifndef CUTTER_PROP_H
#define CUTTER_PROP_H

#include "hand_prop.h"


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
    
    /// unidimensional diffusion coefficient while bound to a Fiber
    real line_diffusion;

    /// rate of cutting event
    real cutting_rate;

    /// max distance from the minus end for cutting
    real cutting_range;
    
    /// amount of polymer removed upon cutting
    real cut_width;

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
    
    float cutting_rate_dt;
    
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
    
    /// perform additional tests for the validity of parameters, given the elasticity
    void checkStiffness(real stiff, real len, real mul, real kT) const;

    /// return a carbon copy of object
    Property* clone() const { return new CutterProp(*this); }

    /// write all values
    void write_values(std::ostream&) const;
    
};

#endif

