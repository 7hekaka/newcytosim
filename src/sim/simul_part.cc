// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#include "simul.h"


/// time_step
real time_step(Simul const& sim) { return sim.prop.time_step; }

/// Boltzmann constant * temperature
real boltzmann(Simul const& sim) { return sim.prop.kT; }

/// ready to simulate
int primed(Simul const& sim) { return sim.primed(); }
