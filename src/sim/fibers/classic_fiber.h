// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef CLASSIC_FIBER_H
#define CLASSIC_FIBER_H

#include "cymdef.h"
#include "vector.h"
#include "fiber.h"
#include "classic_fiber_prop.h"


/// A Fiber following a standard two-state model of dynamic instability
/**
 This implements the 'classical' two-state model of dynamic instability proposed
 by Mitchison and Kirschner:
 - Growing and Shrinking states persist over multiple monomer addition/removal
 - Transitions between these two states, catastrophes and rescues, are stochastic
 .
 
      Dynamic instability of microtubule growth\n
      T. Mitchison and M. Kirschner\n
      https://www.nature.com/articles/312237a0\n
      https://doi.org/10.1073/pnas.82.2.431
 
 For the analysis, see:
 
      Physical aspects of the growth and regulation of microtubule structures\n
      M. Dogterom and S. Leibler\n
      https://doi.org/10.1103/PhysRevLett.70.1347
 
 In this class, both ends may undergo dynamic instability and do so continuously:
 The length is incremented at each time step by `time_step * speed`.
 
 Moreover,
 - the growth speed is linearly proportional to free tubulin concentration.
 - the growing speed is decreased by antagonistic force,
 - the rate of catastrophes is also affected by force,
 - the shrinking speed is constant.
 
 Growth is reduced under antagonistic force exponentially:
 
     Measurement of the Force-Velocity Relation for Growing Microtubules\n
     M. Dogterom and B. Yurke\n
     http://dx.doi.org/10.1126/science.278.5339.856 \n
     http://www.sciencemag.org/content/278/5339/856.abstract
 
 ...and this increases the catastrophe rate:

     Dynamic instability of MTs is regulated by force\n
     M. Janson, M. de Dood, M. Dogterom.\n
     Figure 2 C\n
     http://dx.doi.org/10.1083/jcb.200301147\n
     http://jcb.rupress.org/content/161/6/1029

 See the @ref ClassicFiberPar.

 @ingroup FiberGroup
 */
class ClassicFiber : public Fiber
{
private:
    
    /// state of minus end
    state_t mStateM;

    /// state of plus end
    state_t mStateP;
     
public:
  
    /// constructor
    ClassicFiber(ClassicFiberProp const*);
    
    /// Property
    ClassicFiberProp const* prop() const { return static_cast<ClassicFiberProp const*>(Fiber::prop); }

    /// destructor
    virtual ~ClassicFiber();
        
    //--------------------------------------------------------------------------
    
    /// return assembly/disassembly state of minus end
    state_t endStateM() const { return mStateM; }

    /// return assembly/disassembly state of plus end
    state_t endStateP() const { return mStateP; }

    
    /// change state of minus end
    void setEndStateM(state_t s);
    
    /// change state of plus end
    void setEndStateP(state_t s);

    /// sub simulation step
    real stepMinusEnd();
    
    /// sub simulation step
    real stepPlusEnd();
    
    /// Stochastic simulation
    void step();
    
    //--------------------------------------------------------------------------
    
    /// return specification of fiber class
    std::string activity() const { return "classic"; }

    /// write to Outputter
    void write(Outputter&) const;
    
    /// read from Inputter
    void readEndStates(Inputter&);

    /// read from Inputter
    void read(Inputter&, Simul&, ObjectTag);
    
};


#endif
