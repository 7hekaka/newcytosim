
#include <x86intrin.h>

/// keeping time using Intel's cycle counters
unsigned long long rdt_ = 0;

/// start timer
inline void tic() { rdt_ = __rdtsc(); }

/// return time since last 'tic()'
inline double toc(double arg) { return double(__rdtsc()-rdt_) / arg; }

