// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#ifndef NUCLEATOR_PROP_H
#define NUCLEATOR_PROP_H

#include "hand_prop.h"

class FiberProp;
class FiberSet;


/// additional Property for Nucleator
/**
 @ingroup Properties
 */
class NucleatorProp : public HandProp
{
public:
    
    friend class Nucleator;
    
    /// indicates a specificity
    enum BranchSpecificity
    {
        BRANCH_PARALLEL = 0,
        BRANCH_MOSTLY_PARALLEL,
        BRANCH_ANTIPARALLEL,
        BRANCH_RANDOM,
        BRANCH_SPECIFIED
    };

public:
    
    /**
     @defgroup NucleatorPar Parameters of Nucleator
     @ingroup Parameters
     Inherits @ref HandPar.
     @{
     */
    
    /// rate for nucleation (also known as `nucleate[0]`)
    real nucleation_rate;

    /// name of fiber that is nucleated (also known as `nucleate[1]`)
    std::string fiber_type;
    
    /// specifications of a new fiber (also known as `nucleate[2]`)
    /**
     Options for the newly created Fiber may be specified here:
     see @ref FiberGroup.
     */
    std::string fiber_spec;
    
    /// angle, in radian, of newly made fiber, relative to mother filament
    real branch_angle;
    
    /// if specified, the space on the surface of which new fiber will be created
    std::string nucleate_in_plane;

    /// factor to limit the number of fibers nucleated in total
    float nucleation_limit;
    
    /// specificity for nucleation on the side of existing filament
    int branch_direction;

    /// specifies if the nucleator attaches to fibers that it creates [none, minus_end, plus_end]
    /**
     This option controls if the nucleator will be attached (or not) to a fiber that it created.
     Possible values for `hold_end`:
     - off
     - plus_end
     - minus_end
     .
     Note that a nucleator remains inactive as long as it is bound to a fiber.
     Thus, setting `hold_end = minus_end` in combination with a detachment rate of zero
     will limit nucleation to one fiber at a time, because in that case the nucleator
     can only be released when the fiber is deleted, for example after a catastrophe.
     (default value is `minus_end`)
     */
    FiberEnd hold_end;

    /// option to track a specified end [none, minus_end, plus_end]
    /**
     If `track_end` is set to `plus_end` or `minus_end`, the hand will stay always
     positionned at the given fiber end, even if this end is growing or shrinking.
     Possible values:
     - off
     - plus_end
     - minus_end
     .
     */
    FiberEnd track_end;
    
    /// if true, set the Dynamic State of the nearest filament end to STATE_RED upon detachment
    int addictive;
    
    /// the dynamic state of the fiber tip after detachment
    state_t addictive_state;
    
    /// @}
    
    /// type of fiber to be nucleated
    FiberProp const* fiber_class;
    
    /// if true, new fiber is created in the surface of the Space
    Space const* nucleate_space;

private:
    
    float nucleation_rate_dt;
    
public:
    
    /// constructor
    NucleatorProp(const std::string& n) : HandProp(n)  { clear(); }
    
    /// destructor
    ~NucleatorProp() { }
    
    /// return a Hand with this property
    virtual Hand * newHand(HandMonitor*) const;
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// compute values derived from the parameters
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new NucleatorProp(*this); }

    /// write all values
    void write_values(std::ostream&) const;
   
};

#endif

