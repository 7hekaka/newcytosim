// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SIMUL_PROP_H
#define SIMUL_PROP_H


#include "real.h"
#include "vector.h"
#include "property.h"

class SpaceProp;
class Simul;
class Space;


/**
 Enables capacity to simulate constant fluid flow that transports all objects
 Option normally OFF
 */
#define NEW_CYTOPLASMIC_FLOW 0

/**
 @defgroup Parameters All Object parameters
 List of parameters for user-accessible objects.
 */


/// Property for Simul
/** 
 There is normally only one instantiation of this class.
 
 @ingroup Properties
 */
class SimulProp : public Property
{
    friend class Simul;
    
public:
    
    /**
     @defgroup Properties All Object Properties
     List of properties accessible to the user.
     */

    /**
     @defgroup SimulPar Parameters of Simul
     @ingroup Parameters
     @{
     */
    
    
    /// Current time in the simulated world
    double time;

    /// A small interval of time
    /**
     The `time_step` is the amount of time between two consecutive simulation states.
     It controls the precision of the simulation, at the expense of computation.\n
     We expect that the numerical result will converge to the true mathematical solution
     of the equations when `time_step` becomes infinitely small, but we do not necessarily 
     know how fast this convergence will be. Thus it is not usually possible to 
     calculate the precision as a function of the `time_step`.\n
     To check that `time_step` is appropriate for a particular problem, one should 
     always run several simulations where `time_step` is varied systematically.
     
     Useful rules:
     - A smaller time step is always preferable, provided that the time to run the simulation
     remains acceptable,
     - You should always check that the results of two simulations with <em>time_step=h</em>\
     and <em>time_step=h/2</em> are identical. If this is not the case, then `h`\
     is most likely not an appropriate value for `time_step`,
     - If some reaction in the simulation occures with a rate R, then `time_step`\
     should be such that <em>time_step * R << 1</em>. In practice `time_step * R` should not be higher than 0.2.
     .
     */
    double time_step;

    /// Ambient viscosity
    /**
     The viscosity of the medium should be given in units of pN.s/um^2 = Pa.s.
     <em>The default value is 1</em>
     
     Medium                | Viscosity | Reference                             |
     ----------------------|-----------|----------------------------------------
     Water                 |  0.001    | http://en.wikipedia.org/wiki/Viscosity
     C.elegans embryo      |  ~1       | Daniels et al. 10.1529/biophysj.105.080606  http://dx.doi.org/10.1529/biophysj.105.080606
     C.elegans embryo      |  ~0.5     | Garzon-Coral et al. 10.1126/science.aad9745  http://dx.doi.org/10.1126/science.aad9745
     S.pombe               |  ~1       | Tolic et al. PRL 93, 078102 (2004). http://dx.doi.org/10.1103/PhysRevLett.93.078102
     D.melanogaster embryo |  ~0.312   | Mechanical Aspects of Drosophila Gastrulation, Oleg Polyakov PhD Thesis. Princeton U., 9.2013
     Cleared egg cytoplasm |  ~0.02    | Valentine et al. Biophys J 88, 680–689 (2005). http://dx.doi.org/10.1529/biophysj.104.048025
     Cultivated cells      |  ~1       | Kole et al. Mol Bio Cell 15, 3475--84 (2004) http://dx.doi.org/10.1091/mbc.E04-03-0218
     
     Note that non-linear effects are not taken into account in Cytosim. Hydrodynamic effects are also neglected, such that the drag coefficient of a collection of objects is simply the sum of the individual drag coefficients. However, an `effective` viscosity can be set for each object class, and with this option, it is possible to adjust the drag coefficient of the collection to a realistic value.
     */
    real viscosity;
    
#if NEW_CYTOPLASMIC_FLOW
    /// uniform and constant fluid flow
    Vector uniform_flow;
#endif
    
    /// Energy of Brownian motion in the system = Temperature * Boltzman constant
    /**
     <em>kT</em> is the product of the [Boltzmann constant](http://en.wikipedia.org/wiki/Boltzmann_constant) `k`
     by the [Thermodynamic temperature](http://en.wikipedia.org/wiki/Thermodynamic_temperature) in Kelvin:
     - k = 1.38065 x 10^-23 Joule/Kelvin = 13.8065 x 10^-6 pN.um / Kelvin
     - Kelvin = Celsius + 273.15
     .
     
     Celsius   | Kelvin   |  kT            |
     ----------|----------|-----------------
     ~10 C     |  283 K   |  0.0039  pN.um
     ~24 C     |  297 K   |  0.0041  pN.um
     ~31 C     |  303 K   |  0.0042  pN.um
     ~39 C     |  312 K   |  0.0043  pN.um

     <em>default value = 0.0042</em>
     */
    real kT;
    
    
    /// Desired precision in the motion of the objects
    /**
     The motion of the objects is solved with a residual error that is lower than `tolerance * B`, 
     where `B` is the typical Brownian displacement of the objects in one time step.\n
     <em>Thus one should set 0 < tolerance < 0.1</em>\n
     Lowering tolerance increases precision at the expense of CPU time.
     In the special case where 'kT==0', the maximum residual is simply `tolerance`.
     
     <em>default value = 0.05</em>
    */
    real tolerance;
    
    
    /// 32-bits seed for random number generator
    /**
     The simulation uses SFMT, a fast Mersenne Twister to generate pseudo-random numbers
     http://en.wikipedia.org/wiki/Mersenne_twister
     
     The generator is initialized from `random_seed` specified in the config file,
     but if `random_seed == 0`, it is set automatically during initialization.
    
     <em>default value = 0</em>
     */
    unsigned random_seed;

    
    /// Precision threshold for stochastic events
    /** 
     A warning message is issued for a rate K if:
     
         K * time_step > acceptable_prob
     
     In most implementations, a stochastic event (binding/unbinding) may only occur once
     during a time_step, and this becomes inaccurate if ( K * time_step is not small compared to 1 ).
     
     A user may control the `rate overflow' by adjusting `acceptable_prob` and monitoring the
     warning messages.
     
     <em>default value = 0.5</em>
     */
    real acceptable_prob;
    
    
    /// A flag to enable preconditionning when solving the system of equations
    /**
     This parameter can affect the performance greatly, and it is always a good idea
     to try the different possible values of `precondition`:
     - 0 : do not use preconditionning
     - 1 : use an simple banded preconditionner that is quick to calculate
     - 2 : use a symmetric isotropic block, derived from the full matrix of the system
     - 3 : use a isotropic block, derived from the full matrix of the system
     - 4 : use the diagonal block of the full matrix of the system
     - 5 : same as case `4`, converted to another format optimized for SIMD AVX.
     .
     
     Preconditionning is a technique to speed up convergence of iterative methods,
     that can be used to solve the movements of the Mecables in Cytosim.
     It relies on finding a matrix (the preconditionner) that is approximately
     equal to the inverse of the matrix of the system we want to solve.
     Using a good preconditionner should reduce the number of iterations
     needed to converge to a solution, resulting in a potential overall speedup.
     
     However, calculating the preconditionner itself is costly, and performing
     one iteration with preconditionning is also more expensive than without.
     This generally results in a complex tradeoff, and performance will vary.
     In some cases, using a preconditionner can even degrade preformance,
     in particular if some objects have many vertices.
     
     The different preconditionners represent different tradeoff:
     the approximation is improved from 0 to 6, but the CPU cost also increases.
     
     If there is only one filament in the system, `precondition=0` should perform best.
     With many filaments, trying `precondition = [0, 1, 2, 3, 4, 6]' is a good strategy.
     This can be performed automatically by setting 'solve=auto' in the 'run' command.
     <em>default value = 0</em>
     */
    unsigned precondition;

    /// Number of times a preconditionner block can be used
    /**
     Calculating the preconditionner is more costly than applying it.
     
     Thus, for a system that does not evolve too fast, it can be advantageous
     to keep a preconditionner calculated at some time point, for later times.
     Trying different values, starting from `1` up is advised.
     
     This only apply to 'precondition == 7'
     <em>default value = 2</em>
     */
    unsigned precond_span;
    
    /// A flag to control the engine that implement steric interactions between objects
    int steric_mode;
    
    /// Stiffness for repulsive steric interaction (set as steric[1])
    real steric_stiff_push[2];
    
    /// Stiffness for `attractive` steric interaction (set as steric[2])
    real steric_stiff_pull[2];
    
    /// A rectangular volume where steric happens ( infX infY infZ, supX supY supZ )
    Vector steric_region[2];

    /// Grid size used to determine steric interactions
    /**
     Cytosim uses a divide-and-conquer approach to find pairs of objects that are 
     close enough to interact, based on a dividing the Space with a rectangular grid.
     
     `steric_max_range` defines the minimum size of the cells in the grid.
     A finer grid reduces false positives, but increases the amount of memory used
     by the grid, and the number operations needed to establish and clear the grid.
     
     Thus optimal performance is usually obtained for an intermediate value of
     `steric_max_range`. However `steric_max_range` must remain greater than the
     maximum interaction distance, otherwise some interacting pairs will be missed.
     Experimentation is usually necessary to find the best value.
     
     The maximum distance at which an object may interact with a sibling is its
     diameter. Generally, `steric_max_range` should be greater or equal to the
     sum of the radiuses, of any two objects that may interact.
     In the case of fiber, the `interaction-radius` is a combination of the
     segmentation, and the radius: std::sqrt( (4/3*segmentation)^2 + 4*radius^2 )

     If the parameter is not set, cytosim attempts to set it automatically.
     */
    mutable real steric_max_range;
    
    
    /// Grid size used to determine the attachment of Hand to Fiber
    /**
     Cytosim uses a divide-and-conquer approach to detect which Fibers are near
     a given point, without testing every Fiber. This is necessary to determine
     to which Fiber a Hand may bind. The algorithm is based on partitionning Space
     with a grid of rectangular cells of size `binding_grid_step` (see FiberGrid).

     `binding_grid_step` affects the performance of the algorithm, but not its result.
     Smaller values reduce the number of false positives, but requires more memory
     and housekeeping calculations. These requirements also increase with the physical
     dimensions of the system, to the power DIM (the dimensionality).
     
     If the parameter is not set, cytosim attempts to determine it automatically.
     */
    real binding_grid_step;
    
    /// level of verbosity in output
    unsigned verbose;
    
    /// a flag for custom purpose
    unsigned flag;

    /// Name of configuration file (<em>default = config.cym</em>)
    std::string config_file;
    
    /// Name of output property file (<em>default = properties.cmp</em>)
    std::string property_file;
    
    /// Name of output trajectory file (also known as `system`, <em>default = objects.cmo</em>)
    std::string system_file;
    
    /// If `true`, any pre-existing system file will be cleared (<em>default = true</em>)
    bool clear_system_file;
    
    /// if true, do not read unbound single from file (<em>default = 0</em>)
    int skip_free_single;
    
    /// if true, do not read unbound couple from file (<em>default = 0</em>)
    int skip_free_couple;

    /// Display parameters (see @ref DisplayPar)
    std::string display;

    /// @}
    
    /// the time at which the current run will halt
    mutable double stop_time;
    
    /// the time at which the simulation is terminated
    mutable double end_time;

    /// last message from splash()
    mutable std::string splashed;

    /// this is set to true when 'display' is modified, and to 'false' when it is read
    bool display_fresh;

public:
    
    /// constructor
    SimulProp(const std::string& n) : Property(n), config_file("config.cym") { clear(); }
    
    /// destructor
    ~SimulProp() { }
    
    /// identifies the property
    std::string category() const { return "simul"; }
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);

    /// print some info
    void splash(std::ostream&) const;

    /// check and derive parameters
    void complete(Simul const&);

    /// return a carbon copy of object
    Property* clone() const { return new SimulProp(*this); }

    /// write all values
    void write_values(std::ostream&) const;

};

#endif

