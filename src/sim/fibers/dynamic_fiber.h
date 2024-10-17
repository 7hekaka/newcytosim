// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2021

#ifndef DYNAMIC_FIBER_H
#define DYNAMIC_FIBER_H

#include "cymdef.h"
#include "vector.h"
#include "fiber.h"
#include "dynamic_fiber_prop.h"


/// A Fiber with discrete growth and dynamic instability
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
 - the minus end can be static or shrinking (parameter `shrinking_rate[1]`)
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
    
    /// Gillespie countdown timers for plus end:
    float nextGrowthP;
    float nextHydrolP;
    float nextShrinkP;
    
    /// Gillespie countdown timers for minus end:
    float nextGrowthM;
    float nextHydrolM;
    float nextShrinkM;
    
    /// state of units near the plus end: [0] is terminal, [1] is penultimate unit
    short unitP[2];
    
    /// state of units near the minus end
    short unitM[2];
    
    /// dynamic state of plus end
    state_t mStateP;

    /// dynamic state of minus end
    state_t mStateM;
    
    /// calculate dynamic state from unit states near plus end
    state_t calculateStateP() const;
    
    /// calculate dynamic state from unit states near minus end
    state_t calculateStateM() const;
    
    // convert chewing speed to a unit off rate
    real chewingUnits(int end);
    
    // add unit
    void addUnitM();

    // remove last unit at minus end
    void removeUnitM();
    
    // add unit
    void addUnitP();

    // remove last unit at plus end
    void removeUnitP();

public:
  
    /// constructor
    DynamicFiber(DynamicFiberProp const*);
    
    /// Property
    DynamicFiberProp const* prop() const { return static_cast<DynamicFiberProp const*>(Fiber::prop); }

    /// destructor
    virtual ~DynamicFiber();
        
    //--------------------------------------------------------------------------
    
    /// initialize minus end
    void initM();
    
    /// initialize plus end
    void initP();

    /// return assembly/disassembly state of minus end
    state_t endStateM() const;
    
    /// return assembly/disassembly state of plus end
    state_t endStateP() const;
    
    /// change state of minus end
    void setEndStateM(state_t s);
    
    /// change state of plus end
    void setEndStateP(state_t s);

    //--------------------------------------------------------------------------
    
    /// simulate dynamic instability of plus end
    int stepPlusEnd();

    /// simulate dynamic instability of minus end
    int stepMinusEnd();
    
    /// Stochastic simulation
    void step();
    
    /// calculate the edges for a cut of width `w` around `a` (arguments used for input/output)
    void findSeverEdges(real& a, real& w);

    //--------------------------------------------------------------------------
    
    /// return specification of fiber class
    std::string activity() const { return "dynamic"; }

    /// write to Outputter
    void write(Outputter&) const;
    
    /// read from Inputter
    void readEndStates(Inputter&);

    /// read from Inputter
    void read(Inputter&, Simul&, ObjectTag);
    
};


#endif
