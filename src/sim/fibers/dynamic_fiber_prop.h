// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef DYNAMIC_FIBER_PROP
#define DYNAMIC_FIBER_PROP

#include "fiber_prop.h"

/// Enables code to reduce growth of plus ends that are outside
#define NEW_STALL_OUTSIDE 0

/// additional Property for DynamicFiber
/**
 @ingroup Properties
 */
class DynamicFiberProp : public FiberProp
{
    friend class DynamicFiber;
    
public:
    
    /**
     @defgroup DynamicFiberPar Parameters of DynamicFiber
     @ingroup Parameters
     Inherits @ref FiberPar.
     @{
     */
    
    /// see @ref DynamicFiber
    
    /// Length of discrete units of assembly/disassembly
    real unit_length;
    
    /// Speed of assembly
    real growing_speed[2];
    
    /// Spontaneous disassembly in the growing state, independent of free tubulin concentration
    real growing_off_speed[2];

    /// Characteristic force for polymer assembly (default=+inf)
    /**
     Filament      |   Max. possible    | Measured stall force |
     --------------|--------------------|-----------------------
     Microtubule   | ~ 1GTP per 8nm/13  | ~ 1.33 pN
     Actin         | ~ 1ATP per 4nm/2   | ~ 1 pN

     With 1 ATP bringing ~ 80--100 pN.nm of energy
     
     <em>
     <b>Direct measurement of force generation by actin filament polymerization using an optical trap.</b>\n
     Footer et al.\n
     PNAS vol. 104 no. 7; 2007\n
     http://doi.org/10.1073/pnas.0607052104
     
     <b>Measurement of the Force-Velocity Relation for Growing Microtubules</b>\n
     Marileen Dogterom and Bernard Yurke\n
     Science Vol 278 pp 856-860; 1997\n
     http://www.sciencemag.org/content/278/5339/856.abstract
     </em>
     */
    real growing_force[2];

    /// Hydrolysis rate of G-units, which defines the catastrophe rate
    /**
     Without spontaneous off rate (`growing_off_rate==0`),
     the catastrophe rate is set by
     
         catastrophe_rate = 3 * hydrolysis_rate ^2 / growing_rate;
     
     with `growing_rate = growing_speed / unit_length`
     */
    real hydrolysis_rate[2];

    /// Speed of disassembly
    real shrinking_speed[2];
    
    /// switching rate to the growing state for a fiber shorter than `min_length` (default=0)
    real rebirth_rate[2];

#if NEW_STALL_OUTSIDE
    /// catastrophe rate scaling factor applied if the plus end is outside
    /**
     A value < 1 inhibits catastrophe at the edge; A value > 1 accelerates catastrophes
     */
    real stall_outside;

    /// space used for `stall_outside'
    std::string stall_label;
#endif

    /// The probability of encountering an unhydrolysed tubulin unit while shrinking (default=0)
    real unhydrolyzed_prob[2];
    /// @}
    
private:
    
    /// last message from splash()
    mutable std::string splashed;

    real growing_rate_dt[2];
    real growing_off_rate_dt[2];
    real growing_force_inv[2];
    real hydrolysis_rate_2dt[2];
    real shrinking_rate_dt[2];
    real rebirth_prob[2];
#if NEW_STALL_OUTSIDE
    Space const* stall_space;
#endif

public:
    
    /// constructor
    DynamicFiberProp(const std::string& n) : FiberProp(n) { clear(); }

    /// destructor
    ~DynamicFiberProp() { }
    
    /// return a Fiber with this property
    Fiber* newFiber() const;
    
    /// return a length that is a multiple of the unit_length
    real newFiberLength(Glossary& opt) const;

    /// set default values
    void clear();
       
    /// set using a Glossary
    void read(Glossary&);
    
    /// print some info
    void splash(std::ostream& os) const;

    /// check and derive parameter values
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new DynamicFiberProp(*this); }

    /// write
    void write_values(std::ostream&) const;

};

#endif

