// Melissa E. O’Neill Permuted Congruential Generator
// https://en.wikipedia.org/wiki/Permuted_congruential_generator

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
