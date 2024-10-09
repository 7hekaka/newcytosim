// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#include "simul.h"

/// time in the simulated world
real time(Simul const& sim) { return sim.prop.time; }

/// time step of the simulation engine
real time_step(Simul const& sim) { return sim.prop.time_step; }

/// kT = Boltzmann constant * temperature
real boltzmann(Simul const& sim) { return sim.prop.kT; }

/// flag indicating readiness to simulate
int primed(Simul const& sim) { return sim.primed(); }
