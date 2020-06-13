// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef FIBER_PROP
#define FIBER_PROP

#include "real.h"
#include "vector.h"
#include "property.h"
#include "common.h"
#include "sim.h"

class Field;
class Fiber;
class FiberDisp;
class SingleProp;
class SingleSet;
class Space;


/// compile switches to enable advanced features:
#define OLD_SQUEEZE_FORCE       0
#define NEW_COLINEAR_FORCE      0
#define NEW_FIBER_CHEW          0
#define NEW_FIBER_CONFINE2      1
#define NEW_FIBER_CONFINE_RANGE 0
#define NEW_FIBER_LOOP          0
#define NEW_END_FORCE           1

/// Property for a Fiber
/**
 @ingroup Properties
 */
class FiberProp : public Property
{
    friend class Fiber;
    friend class FiberSet;

public:
    
    /**
     @defgroup FiberPar Parameters of Fiber
     @ingroup Parameters
     These are the parameters for Fiber
     @{
     */
    
    
    /// elastic modulus for bending elasticity
    /**
     The bending elasticity modulus `rigidity` has units of pN * um^2,
     and is related to the persitence length `Lp` via the Boltzman constant and
     absolute temperature (kT = k_B * T):
    
         rigidity = Lp * kT
     
     Many measurments have been made and they agree somewhat:\n
     
     Filament                      |    Lp        |    Rigidity
     ------------------------------|--------------|------------------
     Microtubule                   |   ~ 7200 um  | ~ 30   pN.um^2
     Stabilized Microtubule        |   ~ 2200 um  | ~ 10   pN.um^2
     F-actin                       | ~  9--10 um  | ~ 0.04 pN.um^2
     Phalloidin-stabilized F-actin | ~ 17--18 um  | ~ 0.08 pN.um^2
     
     <em>
     Flexural rigidity of microtubules and actin filaments measured from thermal fluctuations in shape.\n
     Gittes et al.\n JCB vol. 120 no. 4 923-934 (1993)\n
     http://dx.doi.org/10.1083/jcb.120.4.923 \n
     http://jcb.rupress.org/content/120/4/923
     </em>
     
     <em>
     Flexibility of actin filaments derived from thermal fluctuations.\n
     Isambert, H. et al.\n  J Biol Chem 270, 11437–11444 (1995)\n
     http://www.jbc.org/content/270/19/11437
     </em>
     */
    real         rigidity;
    
    
    /// desired distance between vertices
    /**
     `segmentation` is a distance, which affects the precision by which the
     shape of a filament is simulated. Specificially, the number of segments 
     used for a filament of length `L` is the integer `N` that minimizes:

         abs_real( L / N - segmentation )
     
     As a rule of thumb, segmentation should scale with rigidity, depending on 
     the expected magnitude of the forces experienced by the filament:

         segmentation = sqrt(rigidity/force)
         force ~ rigidity / segmentation^2
     
     Generally, a simulation should not be trusted if any filament contains kinks
     (i.e. if the angle between consecutive segments is greater then 45 degrees).
     In that case, the simulation should be redone with a segmentation divided by 2,
     and the segmentation should be reduced until kinks do not appear.
     */
    real         segmentation;
    
    /// Minimum length (this limits the length in some cases)
    real         min_length;
    
    /// Maximum length (this limits the length in some cases)
    real         max_length;

    /// amount of monomer available to make this type of fiber
    /**
     If set, `total_polymer` will limit the total length of all the Fibers of this
     class, by making the assembly rate of the fibers dependent on the amount of
     unused material (ie. 'monomers'):

         assembly_speed = ( 1 - sum_of_all_fiber_length / total_polymer ) * (...)

     This links the assembly for all the fibers within one class, as their assembly
     speed decreases linearly with the total amount of polymer in the cell,
     effectively proportional to the normalized concentration of remaining 'monomers'.
     
     By default `total_polymer = infinite`, and this model is disabled.
     */
    real         total_polymer;
    
    /// if `false`, the fiber will be destroyed if it is shorter than `min_length` (default=`false`)
    bool         persistent;

    /// effective viscosity (if unspecified, simul:viscosity is used)
    /**
     Set the effective `viscosity` to lower or increase the drag coefficient of a particular class of fibers. This makes it possible for example to reduce the total drag coefficient of an aster.
     If unspecified, the global `simul:viscosity` is used.
     */
    real         viscosity;
    
    /// radius used to calculate mobility, corresponding to the radius of the fiber
    real         drag_radius;
    
    /// cut-off on the length of the fiber, above which drag is proportional to length
    real         drag_length;

    /// if true, calculate mobility for a cylinder moving near a immobile planar surface
    /**
     You can select between two possible formulas to calculate viscous drag coefficient:

         if ( fiber:drag_model )
             drag = dragCoefficientSurface();
         else
             drag = dragCoefficientCylinder();

     <hr>
     @copydetails Fiber::dragCoefficientCylinder
     <hr>
     @copydetails Fiber::dragCoefficientSurface
     */
    int          drag_model;
    
    /// distance of fluid between immobile surface and cylinder (set as `drag_model[1]`)
    real         drag_gap;

    
    /// can be set to control which Hands may bind
    /**
     This option limits the binding of Hands to this class of Fiber. To decide
     if a Hand may bind or not, the `binding_keys` of the Hand is compared to
     the `binding_key` of the Fiber, using a BITWISE AND:

         allowed = ( fiber:binding_key BITWISE_AND hand:binding_key )

     Hence attachement if the two `binding_key` do not share any common digit
     in base 2. For most usage, one can thus use powers of 2 to distinguish fibers:
     - microtubule: `binding_key = 1`,
     - actin: `binding_key = 2`,
     - etc.
     .
     However, more complex combinations can be created by using all the bits of binding_key.
     With the example above, a Hand with `binding_key=3` can bind to both type of fibers.
     */
    unsigned int binding_key;
    
    /// if true, a Lattice is associated to this fiber
    int          lattice;
    
    /// unit length associated with Lattice
    real         lattice_unit;
    
    /// if true, associate an analog lattice
    int          mesh;
    
    /// unit length associated with the analog Lattice
    real         mesh_unit;
    
    /// if true, the quantities in the lattice can cut the fiber
    int          mesh_cut_fiber;

    /// flux speed of substance on Lattice (speed<0 is MINUS_END directed)
    real         mesh_flux_speed;
    
    /// loading rate of substance from Field to Lattice
    /**
     This is a binding rate per unit time and per unit length of Fiber.
     Binding is proportional to the concentration of substance in the field.
     */
    real         mesh_binding_rate;

    /// unloading rate of substance from Lattice to Field (unit is 1/second)
    real         mesh_unbinding_rate;
    
    real         mesh_aging_rate;
    
    /// flag controlling the forces exerted by Space on fiber points
    /**
     Possible values:
     - `off` (default)
     - `on` or `surface`
     - `inside`
     - `outside`
     - `plus_end`
     - `minus_end`
     - `both_ends`
     .
     */
    Confinement  confine;
    
    /// stiffness of confinement (also known as `confine[1]`)
    real         confine_stiffness;
    
    /// name of space used for confinement (also known as `confine[2]`)
    std::string  confine_space;
    
#if NEW_FIBER_CONFINE2
    /// flag controlling the forces exerted by Space on fiber points
    /**
     Possible values:
     - `off` (default)
     - `on` or `surface`
     - `inside`
     - `outside`
     .
     */
    Confinement  confine2;
    
    /// stiffness of confinement (also known as `confine[1]`)
    real         confine2_stiffness;
    
    /// name of space used for confinement (also known as `confine[2]`)
    std::string  confine2_space;
#endif
#if NEW_FIBER_CONFINE_RANGE
    real         confine_range[2];
#endif
    
    /// if true, include steric interaction for this object
    /**
     The steric interaction generates a force derived from the potential energy:
     
         E = 1/2 k * ( d - d_0 ) ^ 2
     
     where `d` is the distance between two sections of filament. 
     The force is controlled by two parameters:
     - a stiffness `k`,
     - and equilibrium length `d_0`
     .
     
     This force is repulsive at short range ( d < d_0 ),
     and attractive elsewhere ( d > d_0 ).
     */
    int          steric;
    
    /// radius of repulsive steric interaction (also known as `steric[1]`)
    real         steric_radius;
    
    /// extra radius of attractive steric interaction (also known as `steric[2]`)
    real         steric_range;
    
    /// name of field
    std::string  field;
    
    /// type of glue (interaction between fiber PLUS_END and Space)
    /**
     Parameter fiber:glue is used to create interactions with the boundaries:
     - it creates a Single, everytime a fiber contacts the surface.
     - the Single is deleted if the associated Hand detaches.
     .
    */
    int          glue;
    
    /// name of Single used for glue (set a `glue[1]`)
    std::string  glue_single;
    
#if NEW_COLINEAR_FORCE
    /// a force parallel to the fiber (force per fiber length)
    /**
     This has unit of force per unit length:
     - a positive 'colinear_force' is directed toward the plus end,
     - a negative 'colinear_force' is directed toward the minus end.
     .
     */
    real         colinear_force;
#endif
#if NEW_FIBER_CHEW
    /// maximum speed of disassembly due to chewing (speed)
    real         max_chewing_speed;
#endif
    
    /// specialization
    /**
     @copydetails FiberGroup
     */
    std::string  activity;
    
    /// display string (see @ref FiberDispPar)
    std::string  display;
    
#if OLD_SQUEEZE_FORCE
    /// add a force toward the X-axis
    int  squeeze;
    /// max norm of squeezing force (set as \c squeeze[1])
    real squeeze_force;
    /// range below which squeezing is linear (set as \c squeeze[2])
    real squeeze_range;
#endif
    
#if NEW_END_FORCE
    /// the force vector added to an end of the fiber
    Vector       end_force;
    /// the fiber end to which a force is added (set as end_force[1])
    FiberEnd     end_force_mode;
#endif
    
#if NEW_FIBER_LOOP
    /// if `true`, link MINUS and PLUS ends together to form a loop
    bool         loop;
#endif
    /// @}

    /// derived variable: flag to indicate that `display` has a new value
    bool         display_fresh;
    
    /// derived variable: display
    FiberDisp *  disp;
    
    /// pointer to actual confinement Space, derived from `confine_space`
    Space const* confine_space_ptr;
   
#if NEW_FIBER_CONFINE2
    /// pointer to actual confinement Space, derived from `confine_space`
    Space const* confine2_space_ptr;
#endif
    
    /// derived variable: pointer to associated Field
    Field *      field_ptr;

protected:
    
    /// maximum speed of shrinkage
    real    max_chewing_speed_dt;
    
    /// fraction of unpolymerized monomers in [0, 1]
    real    free_polymer;
    
    /// total length of fiber for this type
    mutable real used_polymer;
    
    /// SingleProp used for glue
    SingleProp * glue_prop;

public:
    
    /// constructor
    FiberProp(const std::string& n) : Property(n), disp(nullptr) { clear(); }
    
    /// destructor
    ~FiberProp() { }
    
    /// return a non-initialized Fiber with this property
    virtual Fiber* newFiber() const;
    
    /// return a Fiber with this property, initialized
    Fiber* newFiber(Glossary& opt) const;
    
    /// identifies the property
    std::string category() const { return "fiber"; }
    
    /// set default values
    virtual void clear();
       
    /// set using a Glossary
    virtual void read(Glossary&);
   
    /// check and derive parameter values
    virtual void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new FiberProp(*this); }

    /// write
    virtual void write_values(std::ostream&) const;

};

#endif

