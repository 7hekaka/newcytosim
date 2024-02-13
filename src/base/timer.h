// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#ifndef TIMER_H
#define TIMER_H

#include <time.h>
#include <sys/time.h>

// using computer's clock
struct timespec tic_t;

/// return current time in milliseconds
inline unsigned long timer()
{
    timespec tv;
    clock_gettime(CLOCK_MONOTONIC, &tv);
    return 1000 * (unsigned long)tv.tv_sec + tv.tv_nsec / 1e6;
}

/// start timer
inline void tick()
{
    clock_gettime(CLOCK_MONOTONIC, &tic_t);
}

/// return time since last 'tick()', divided by 'arg'
inline double tock(double arg = 1)
{
    timespec tv;
    clock_gettime(CLOCK_MONOTONIC, &tv);
    return double(1000*(tv.tv_sec-tic_t.tv_sec) + (tv.tv_nsec-tic_t.tv_nsec)/1e6)/arg;
}


#if ( 0 )
// using the CPU time

#include <ctime>
clock_t tic_t;

/// return current time value
inline clock_t timer()
{
    return clock();
}

/// start timer
void tick()
{
    tic_t = clock();
}

/// return time since last 'tick()', divided by 'arg'
double tock(double arg = CLOCKS_PER_SEC)
{
    return (1e3*( clock() - tic_t )) / arg;
}
#endif

#if ( 0 )

/*
 Using Intel's processor Cycle counting functions,
 multiplying by 2^20 to get into the milli-second range
 but exact conversion depends on the processor frequency and may vary
*/

#include <x86intrin.h>

/* declare rdtsc from assembly if necessary
extern __inline unsigned long long
__attribute__((__gnu_inline__, __always_inline__, __artificial__))
__rdtsc (void)
{
  return __builtin_ia32_rdtsc();
}
*/

/// keeping time using Intel's cycle counters
unsigned long long rdt_ = 0;

/// return current timer value (64 bits)
inline unsigned long long timer() { return __rdtsc() >> 20; }

/// start timer
inline void tick() { rdt_ = __rdtsc(); }

/// return time since last 'tick()', divided by 'arg'
inline double tock(double arg = 1) { return double((__rdtsc()-rdt_) >> 20) / arg; }

#endif

#endif
