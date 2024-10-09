// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.
/**
 This file declares the `Simul` object, and global functions, providing access
 to parameter values that are used frequently across Cytosim's code.
 This avoids including the complete definition of class `Simul` in most case,
 reducing the circular dependencies that would arise if we did so.
 */

#include "real.h"

class Simul;

/// time in the simulated world (@returns SimulProp::time)
real time(Simul const&);

/// time step of the simulation engine (@returns SimulProp::time_step)
real time_step(Simul const&);

/// kT = Boltzmann constant * temperature (@returns SimulProp::kT)
real boltzmann(Simul const&);

/// flag indicating readiness to simulate (@returns Simul::primed)
int primed(Simul const&);

