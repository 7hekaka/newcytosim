// Melissa E. O’Neill Permuted Congruential Generator
// https://en.wikipedia.org/wiki/Permuted_congruential_generator

#ifndef RANDOM_PCG
#define RANDOM_PCG

#include <stdint.h>

static uint64_t const multiplier = 6364136223846793005u;
static uint64_t const increment  = 1442695040888963407u; // or an arbitrary odd constant

static inline uint32_t rotr32(uint32_t x, unsigned r)
{
    return x >> r | x << (-r & 31);
}

inline uint32_t pcg32(uint64_t& state)
{
    uint64_t x = state;
    unsigned count = (unsigned)(x >> 59);        // 59 = 64 - 5

    state = x * multiplier + increment;
    x ^= x >> 18;                                // 18 = (64 - 27)/2
    return rotr32((uint32_t)(x >> 27), count);   // 27 = 32 - 5
}

inline uint64_t pcg32_init(uint64_t seed)
{
    uint64_t state = seed + increment;
    (void)pcg32(state);
    return state;
}


/**
 FJN's function to distribute bits
 returns a random integer with exactly `b` bits equal to `1`, randomly positionned.
 */
inline uint32_t distribute_bits(unsigned b, uint64_t& state)
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

#endif
