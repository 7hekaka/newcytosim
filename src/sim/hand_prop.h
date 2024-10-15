// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef HAND_PROP
#define HAND_PROP

#include "cymdef.h"
#include "real.h"
#include "property.h"
#include "fiber.h" // for FiberLattice


/// enables "bind_only_free_end" to limit binding of Hands to Fibers
#define NEW_BIND_ONLY_FREE_END 0

/**
 Set below POOL_UNATTACHED = N to pool N successive Hand's attachment trials.
 The effective diffusion coefficient and binding rates are increased to compensate.
 This can be advantageous since FiberGrid::paintGrid() is spared for N-1 steps
 as the attachment algorithm runs only once every N steps.

 Define POOL_UNATTACHED as 1 to disable the feature, and 4, 8 or 16 to enable it.
 */
#define POOL_UNATTACHED 1

class Hand;
class HandMonitor;
class PointDisp;

/// Property for Hand
/**
 @ingroup Properties
*/
class HandProp : public Property
{
    friend class Hand;
    
public:
      
    /// return one of the Property derived from HandProp
    static HandProp * newProperty(const std::string& n, Glossary&);
    
public:
    
    /**
     @defgroup HandPar Parameters of Hand
     @ingroup Parameters
     @{
     */
    
    /// rate of attachment when Hand is within `binding_range` (also known as `binding[0]`)
    /**
     This has units of 1/second.
     The molecular binding_rate of conventional kinesin is 4.7 +/- 2.4 /s:
         Leduc et al. PNAS 2004 vol. 101 no. 49 17096-17101
         http://dx.doi.org/10.1073/pnas.0406598101 \n
         http://www.pnas.org/content/101/49/17096.abstract
     
     <em>default value = 0</em>
     */
    real binding_rate;
    
    
    /// maximum distance at which the Hand can bind (also known as `binding[1]`)
    real binding_range;
    
    
    /// This bitfield can be set to restrict binding to certain type of Fiber
    /**
     The binding to a fiber is allowed only if the keys of the Hand and Fiber match.
     The test uses a BITWISE_AND operation of the two keys:

         if ( fiber:binding_key BITWISE_AND hand:binding_key )
             allowed = true;
         else
             allowed = false;

     Thus, with `binding_key = 0` attachment is completely disabled, and in general
     one needs to look at the bit-pattern of the keys in base 2, to determine if two keys are compatible.
     For example, `1` [binary: 01] can bind to `3` [binary: 11]  but not to `2` [binary: 10].
     
     <em>default value = all-bits-at-1</em>
     */
    unsigned binding_key;
    
    
    /// detachment rate in the absence of load (also known as `unbinding[0]`)
    /**
     This defines a detachment opportunity that is proportional to time.
     Kramers theory specifies that the detachment rate depends on the force
     in the link:
     
         RATE = unbinding_rate * exp( FORCE / unbinding_force )
     
     where FORCE is the norm of the tension in the link holding the Hand,
     and `unbinding_rate' and `unbinding_force' are two parameters.
     By setting `unbinding_force=inf', unbinding is made independent of load.
     
     Two articles:
         Mechanics of the kinesin step
         Carter, N. & Cross, R. Nature 435, 308–312 (2005).
         http://dx.doi.org/doi:10.1038/nature03528
     and
         Examining kinesin processivity within a general gating framework
         Andreasson et al. eLife 2015;4:e07403
         http://dx.doi.org/10.7554/eLife.07403
     provide similar values for conventional kinesin:

         unbinding_rate = 1 / s
         unbinding_force ~ 2 pN
     
     <em>default value = 0</em>
     (see @ref Stochastic)
     */
    real unbinding_rate;
    
    
    /// characteristic force of unbinding (also known as `unbinding[1]`)
    /**
     @copydetails unbinding_rate
     
     <em>default value = inf</em>
     */
    real unbinding_force;
    
    
    /// if true, the Hand can also bind directly to the tip of fibers
    /**
     The value of `bind_also_end` determines binding from a position that does not
     project inside a fiber backbone, for which the closest point on the fiber
     is one of its tip. Even when attachment is permitted to the tips of the fiber,
     the distance must be shorter than `binding_range` for binding to occur.
     
     By setting 'bind_also_end', you can extend the capture regions of the fibers
     to include one or two hemi-spheres, of radius `binding_range`, at the tips of
     the fibers.
     
     Values: `bind_also_end = {off, minus_end, plus_end, both_ends}`.
     <em>default value = off</em>
     */
    int bind_also_end;
    
    
    /// if true, the Hand can bind only near the tips of the fibers
    /**
     This determines that a Hand can only bind near the ends of the fiber.
     This parameter can be 'none', 'plus_end', 'minus_end' , 'center' or 'both_ends'.
     Binding is allowed on positions located within a distance 'bind_end_range'
     from the specified end ('bind_end_range' is specified as `bind_only_end[1]`).
     
     <em>default value = off</em>
     */
    FiberEnd bind_only_end;
    
    
    /// length cutoff associated with `bind_only_end` where hand may bind (set as `bind_only_end[1]`)
    real bind_end_range;

#if NEW_BIND_ONLY_FREE_END
    /// if true, only bind fiber tip if no other hand is bound already
    bool bind_only_free_end;
#endif
    
    /// detachment parameter, for cases when the hand reaches a growing or a static fiber end
    /**
     A Hand may reach the tip of the fiber on which it is bound, because it has moved,
     and `hold_growing_end` will determine the probability of detachment in this case.
     - if `hold_growing_end = 0`, the hand will detach,
     - if `hold_growing_end = 1`, the hand will be placed exactly at the fiber's end.
     `hold_growing_end` is a probability that must be in [0, 1].
     There are two values: [0] applies to the plus end, and [1] is for the minus end
     
     <em>default = 0, 0</em>
     */
    real hold_growing_end[2];
    
    
    /// detachment parameter, for cases when the Hand is reached by a shrinking fiber end
    /**
     This determines detachment/no-detachment if a Hand is reached by a shrinking fiber end:
     - if `hold_shrinking_end = 0`, the hand will detach,
     - if `hold_shrinking_end = 1`, the hand will be relocated to track the end.
     `hold_shrinking_end` is a probability that must be in [0, 1].
     There are two values: [0] applies to the plus end, and [1] is for the minus end
     To set the hand to track shrinking minus end, use `hold_shrinking_end == 0, 1`
     <em>default = 0, 0</em>
     */
    real hold_shrinking_end[2];
    
    
    /// specialization
    /**
     @copydetails HandGroup
     */
    std::string activity;
    
    
    /// display parameters (see @ref PointDispPar)
    std::string display;
    
    
    /// size of one step
    real step_size;
    
    /// specifies the position of binding within the Lattice cell
    /**
     `site_shift` should be in [0, step_size]:
     - at `0.0`, the attachment position is at the start of the site
     - at `step_size`, the attachment position is at the end of the site
     - at `step_size/2`, the attachment is midway
     [default = step_size/2]
     */
    real site_shift;
    
    /// list of cell's bits covered upon binding to the lattice
    FiberLattice::cell_t footprint;

    /** @} */

public:

    /// derived variable: 1.0 / unbinding_force
    real unbinding_force_inv;
    
    /// derived variable: log(unbinding_rate * time_step)
    real unbinding_rate_log;

    /// derived variable = probability to bind in one time step;
    real binding_prob;
    
    /// derived variable = unbinding_rate * timestep;
    float unbinding_rate_dt;
    
    /// flag to indicate that `display` has a new value
    bool display_fresh;
    
public:

    /// the display parameters for this category of Hand
    PointDisp * disp;
    
    /// constructor
    HandProp(const std::string& n) : Property(n), disp(nullptr) { clear(); }
    
    /// destructor
    ~HandProp() { }
    
    /// return a Hand with this property
    virtual Hand * newHand(HandMonitor*) const;

    /// identifies the property
    std::string category() const { return "hand"; }
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    virtual void read(Glossary&);
    
    /// compute values derived from the parameters
    virtual void complete(Simul const&);
    
    /// perform additional tests for the validity of parameters, given the elasticity
    virtual void checkStiffness(real stiff, real len, real mul, real kT) const;

    /// return 'unload_speed' for the Motor class
    virtual real motorSpeed() const { return 0; }
    
    /// Attachment rate per unit length of fiber
    real bindingSectionRate() const;
    
    /// Attachment probability per unit length of fiber in one time step
    real bindingSectionProb() const;

    /// write all values
    void write_values(std::ostream&) const;
    
    /// return a carbon copy of object
    Property* clone() const { return new HandProp(*this); }
};


#endif

