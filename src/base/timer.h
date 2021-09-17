
#include <x86intrin.h>

#if 0
/* rdtsc */
extern __inline unsigned long long
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
__rdtsc (void)
{
  return __builtin_ia32_rdtsc();
}
#endif

/// keeping time using Intel's cycle counters
unsigned long long rdt_ = 0;

/// start timer
inline void tic() { rdt_ = __rdtsc(); }

/// return time since last 'tic()'
inline double toc(double arg) { return double(__rdtsc()-rdt_) / arg; }

