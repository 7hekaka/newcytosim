// Cytosim was created by Francois Nedelec. Copyright 2024 Cambridge University.

#include "random_pcg.h"

#include <limits.h>
#include <sys/time.h>
#include <iostream>
#include <random>

namespace PCG32
{
    /// this is a global random number generator
    /**
     The Cytosim simulation engine use `RNG`, which exist independently for each thread.
     The pcg32_state is used by the master thread to coordinate multiple instances
     of the RNG, in multiplay and other App which can run mulplitle independent Cytosims
     This must be initialized!
     */
    uint64_t pcg32_state = 7;
    
    
    /**
     FJN's function to distribute bits
     returns a random integer with exactly `b` bits equal to `1`, randomly positionned.
     */
    uint32_t distribute_bits(unsigned b, uint64_t& state)
    {
        if ( b > 16 )
        {
            if ( b > 31 )
                return ~0U;
            return ~distribute_bits(32-b, state);
        }
        uint32_t i = 0;
        while ( b > 0 )
        {
            uint32_t s = pcg32(state);
            uint32_t x = 1 << ( s & 31 );
            if (!( i & x ))
            {
                i |= x;
                --b;
            }
        }
        return i;
    }
    
}
