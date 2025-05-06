// Melissa E. O’Neill Permuted Congruential Generator
// https://en.wikipedia.org/wiki/Permuted_congruential_generator

#ifndef RANDOM_PCG
#define RANDOM_PCG

#include <cstdint>

/**
 Cytosim's simulation engine uses `RNG`, which exists independently for each thread.
 PCG32 is used by the master thread to differently seed multiple instances of RNG,
 in multiplay and other Apps which can run many independent Cytosim simulations
 */
namespace PCG32
{
    /// this is a global random number generator
    extern uint64_t pcg32_state;
    
    /// Find random seed from a random device, or from the system's clock
    extern uint32_t get_random_seed();

    /// Melissa E. O’Neill Permuted Congruential Generator
    static uint64_t const multiplier = 6364136223846793005u;
    
    /// Melissa E. O’Neill Permuted Congruential Generator
    static uint64_t const increment  = 1442695040888963407u; // or an arbitrary odd constant
    
    /// Melissa E. O’Neill Permuted Congruential Generator
    static inline uint32_t rotr32(uint32_t x, unsigned r)
    {
        return x >> r | x << (-r & 31);
    }
    
    /// Melissa E. O’Neill Permuted Congruential Generator
    inline uint32_t pcg32(uint64_t& state)
    {
        uint64_t x = state;
        unsigned count = (unsigned)(x >> 59);        // 59 = 64 - 5
        
        state = x * multiplier + increment;
        x ^= x >> 18;                                // 18 = (64 - 27)/2
        return rotr32((uint32_t)(x >> 27), count);   // 27 = 32 - 5
    }
    
    /// Melissa E. O’Neill Permuted Congruential Generator
    inline uint64_t pcg32_init(uint64_t seed)
    {
        uint64_t state = seed + increment;
        (void)pcg32(state);
        return state;
    }
    
    
    /// intialize global `pcg32_state`
    inline void seed_pcg32()
    {
        pcg32_state = get_random_seed();
    }

    /// return next random integer using global `pcg32_state`
    inline uint32_t pcg32()
    {
        return pcg32(pcg32_state);
    }

    ///returns a random integer with exactly `b` bits equal to `1`, randomly positionned.
    uint32_t distribute_bits(unsigned b, uint64_t& state);
    
}
#endif
