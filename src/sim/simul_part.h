// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.
/**
 This file declares the Simul object, and some global functions, permitting
 access to values that are used frequently across Cytosim's code.
 This avoids including the full definition of Simul in most case, breaking
 the circular dependencies that would arise if we did so.
 */

#include "real.h"

class Simul;

/// returns SimulProp::time_step
real time_step(Simul const&);

/// returns SimulProp::kT
real boltzmann(Simul const&);

/// returns Simul::primed()
int primed(Simul const&);

