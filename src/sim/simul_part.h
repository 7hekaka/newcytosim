// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.
/**
 This provides an incomplete declaration of Simul, and some global functions
 to access values that are used frequently across Cytosim's code
 */

#include "real.h"

class Simul;

/// returns SimulProp::time_step
real time_step(Simul const&);

/// returns SimulProp::kT
real boltzmann(Simul const&);

/// returns Simul::primed()
int primed(Simul const&);

