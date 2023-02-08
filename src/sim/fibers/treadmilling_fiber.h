// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef TREADMILLING_FIBER_H
#define TREADMILLING_FIBER_H

#include "cymdef.h"
#include "vector.h"
#include "fiber.h"
#include "treadmilling_fiber_prop.h"


/// A Fiber with assembly at both ends 
/**
 The growing speed of each end are set independently.
 The basic parameters are:
 
 * `growing_speed`, the base assembly rate in um/s.
 * `growing_force`, the characteristic force of polymerization in pN.
 
 Positive values correspond to assembly, and negative values to disassembly.
 Assembly is exponentially decreased by antagonistic force, and linearly dependent
 on the availability of polymer.  Disassembly always occurs at the specified rate.
 Only the component of the force parallel to the direction of the fiber at the end
 is taken into account:
 
     force = force_vector * fiber_direction;
 
 The projected force is negative ( antagonistic ) if it is directed against fiber assembly.
 
     if ( force < 0 )
         speed = growing_speed * free_polymer * exp(force/growing_force);
     else
         speed = growing_speed * free_polymer;
 
 In this equation, `free_polymer` is a number in [0,1], representing the fraction of free monomers.
 It is defined as:
 
    free_polymer = 1.0 - sum(all_fiber_length) / total_polymer
 
 The length of a fiber will not exceed `fiber:max_length`,
 Fiber shorter than `fiber:min_length` are deleted if `fiber:persistent == 0`.
 
 See the @ref TreadmillingFiberPar.

 @ingroup FiberGroup
 */
class TreadmillingFiber : public Fiber
{   
private:
    
    /// state of plus end
    state_t mStateP;
    
    /// state of minus end
    state_t mStateM;
    
public:
  
    /// constructor
    TreadmillingFiber(TreadmillingFiberProp const*);
    
    /// Property
    TreadmillingFiberProp const* prop() const { return static_cast<TreadmillingFiberProp const*>(Fiber::prop); }
    
    /// destructor
    virtual ~TreadmillingFiber();
        
    //--------------------------------------------------------------------------
    
    /// return assembly/disassembly state of minus end
    state_t endStateM() const { return mStateM; }
    
    /// change state of minus end
    void setEndStateM(state_t s);

    
    /// return assembly/disassembly state of plus end
    state_t endStateP() const { return mStateP; }

    /// change state of plus end
    void setEndStateP(state_t s);
    
    //--------------------------------------------------------------------------
    
    /// Stochastic simulation
    void step();
    
    //--------------------------------------------------------------------------
    
    /// return specification of fiber class
    std::string activity() const { return "treadmill"; }

    /// write to Outputter
    void write(Outputter&) const;
    
    /// read from Inputter
    void readEndStates(Inputter&);

    /// read from Inputter
    void read(Inputter&, Simul&, ObjectTag);
    
};


#endif
