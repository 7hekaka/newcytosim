// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2021

#ifndef DYNAMIC_FIBER_H
#define DYNAMIC_FIBER_H

#include "cymdef.h"
#include "vector.h"
#include "fiber.h"

class DynamicFiberProp;


/// A Fiber with discrete growth and dynamic instability at the PLUS_END
/**
 This implements the microtubule dynamic instability model proposed by
 Brun, Rupp et al. with a 'hard-coded' coupling parameter N=2.
 
 Assembly and disassembly follow discrete steps of size `prop->unit_length`.
 The model keeps track of the state of the two terminal units (ie. tubulin heterodimers).
 This leads to 4 different states, which are mapped to states [GREEN YELLOW ORANGE RED].
 
 
 The growth speed is reduced under antagonistic force by an exponential factor:
 
***Measurement of the Force-Velocity Relation for Growing Microtubules***
Marileen Dogterom and Bernard Yurke
Science Vol 278 pp 856-860; 1997
http://www.sciencemag.org/content/278/5339/856.abstract
 
...and this will increase the catastrophe rate:
 
***Dynamic instability of MTs is regulated by force***
M.Janson, M. de Dood, M. Dogterom.
Journal of Cell Biology Vol 161, Nb 6, 2003
Figure 2 C
http://www.jcb.org/cgi/doi/10.1083/jcb.200301147
 
 
 If you use this model, please cite:
 
    ***A theory of microtubule catastrophes and their regulation***
    Brun L, Rupp B, Ward J, Nedelec F
    PNAS 106 (50) 21173-21178; 2009
    http://www.pnas.org/content/106/50/21173

 The predicted mean time until catastrophe is approximately

    growing_rate = growing_speed / unit_length
    real ctime = growing_rate / ( 3 * hydrolysis_rate * hydrolysis_rate );
 
 The implemented model includes off-rate in the assembly state, as described in:
 
    ***Random Hydrolysis Controls the Dynamic Instability of Microtubules***
    Ranjith Padinhateeri, Anatoly B Kolomeisky, and David Lacoste
    Biophys J 102, 1274–1283 (2012)
    http://dx.doi.org/10.1016/j.bpj.2011.12.059
 
 Moreover:
 
 - Gillespie timers are used for the stochastic model
 - the MINUS_END can be static or shrinking (parameter `shrinking_rate[1]`)
 - a simple rescue mechanism was implented as unhydrolyzed_prob[] (Maud Formanek)
 .
 
 See the @ref DynamicFiberPar.

 @todo DynamicFiber detach_rate could depend on the state of the subunit
 @todo DynamicFiber could keep the entire state vector of the subunits

 Note:
 This class is not fully tested
 @ingroup FiberGroup
 */
class DynamicFiber : public Fiber
{
private:
    
    /// assembly during last time-step
    real mGrowthP;
    real mGrowthM;
    
    /// Gillespie countdown timers for PLUS_END:
    real nextGrowthP;
    real nextHydrolP;
    real nextShrinkP;
    
    /// Gillespie countdown timers for MINUS_END:
    real nextGrowthM;
    real nextHydrolM;
    real nextShrinkM;
    
    /// state of units near the PLUS_END: [0] is terminal, [1] is penultimate unit
    unsigned unitP[2];
    
    /// state of units near the MINUS_END
    unsigned unitM[2];
    
    /// dynamic state of PLUS_END
    state_t mStateP;

    /// dynamic state of MINUS_END
    state_t mStateM;

    /// calculate dynamic state from unit states near PLUS_END
    state_t calculateStateP() const;
    
    /// calculate dynamic state from unit states near MINUS_END
    state_t calculateStateM() const;
   
public:
    
    /// the Property of this object
    DynamicFiberProp const* prop;
  
    /// constructor
    DynamicFiber(DynamicFiberProp const*);

    /// destructor
    virtual ~DynamicFiber();
        
    //--------------------------------------------------------------------------
    
    /// initialize minus end
    void initM();
    
    /// initialize plus end
    void initP();

    /// return assembly/disassembly state of MINUS_END
    state_t endStateM() const;
    
    /// return assembly/disassembly state of PLUS_END
    state_t endStateP() const;
    
    /// change state of MINUS_END
    void setEndStateM(state_t s);
    
    /// change state of PLUS_END
    void setEndStateP(state_t s);
    
    /// length increment at MINUS_END during last time-step
    real freshAssemblyM() const { return mGrowthM; }
    
    /// length increment at PLUS_END during last time-step
    real freshAssemblyP() const { return mGrowthP; }

    //--------------------------------------------------------------------------
    
    /// simulate dynamic instability of PLUS_END
    int stepPlusEnd();
    
    /// simulate dynamic instability of MINUS_END
    int stepMinusEnd();
    
    /// Stochastic simulation
    void step();
    
    //--------------------------------------------------------------------------
    
    /// return specification of fiber class
    std::string activity() const { return "dynamic"; }

    /// write to Outputter
    void write(Outputter&) const;
    
    /// read from Inputter
    void readEndState(Inputter&);

    /// read from Inputter
    void read(Inputter&, Simul&, ObjectTag);
    
};


#endif
