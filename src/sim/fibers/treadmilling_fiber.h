// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef TREADMILLING_FIBER_H
#define TREADMILLING_FIBER_H

#include "cymdef.h"
#include "vector.h"
#include "fiber.h"

class TreadmillingFiberProp;


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
    
    /// state of PLUS_END
    state_t    mStateP;
    
    /// assembly during last time-step
    real       mGrowthP;
    
    /// state of MINUS_END
    state_t    mStateM;
    
    /// assembly during last time-step
    real       mGrowthM;
    
public:
    
    /// the Property of this object
    TreadmillingFiberProp const* prop;
  
    /// constructor
    TreadmillingFiber(TreadmillingFiberProp const*);

    /// destructor
    virtual ~TreadmillingFiber();
        
    //--------------------------------------------------------------------------
    
    /// return assembly/disassembly state of MINUS_END
    state_t     endStateM() const { return mStateM; }
    
    /// change state of MINUS_END
    void        setEndStateM(state_t s);
    
    /// length increment at MINUS_END during last time-step
    real        freshAssemblyM() const { return mGrowthM; }

    
    /// return assembly/disassembly state of PLUS_END
    state_t     endStateP() const { return mStateP; }

    /// change state of PLUS_END
    void        setEndStateP(state_t s);
    
    /// length increment at PLUS_END during last time-step
    real        freshAssemblyP() const { return mGrowthP; }
    
    //--------------------------------------------------------------------------
    
    /// Stochastic simulation
    void        step();
    
    //--------------------------------------------------------------------------
    
    /// return specification of fiber class
    std::string activity() const { return "treadmill"; }

    /// write to Outputter
    void        write(Outputter&) const;
    
    /// read from Inputter
    void        readEndState(Inputter&);

    /// read from Inputter
    void        read(Inputter&, Simul&, ObjectTag);
    
};


#endif
