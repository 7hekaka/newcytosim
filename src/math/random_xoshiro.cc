#include <stdint.h>


typedef uint64_t splitmix64_state;

uint64_t splitmix64(splitmix64_state& state)
{
    uint64_t result = (state += 0x9E3779B97f4A7C15);
    result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9;
    result = (result ^ (result >> 27)) * 0x94D049BB133111EB;
    return result ^ (result >> 31);
}


typedef uint64_t xoshiro256ss_state[4];

// one could do the same for any of the other generators
void xoshiro256_init(xoshiro256ss_state& state, uint64_t seed)
{
    splitmix64_state s(seed);

    state[0] = splitmix64(s);
    state[1] = splitmix64(s);
    state[2] = splitmix64(s);
    state[3] = splitmix64(s);
}


inline uint64_t rol64(uint64_t x, int k)
{
    return (x << k) | (x >> (64 - k));
}

uint64_t xoshiro256ss(xoshiro256ss_state& state)
{
    uint64_t const result = rol64(state[1] * 5, 7) * 9;
    uint64_t const t = state[1] << 17;

    state[2] ^= state[0];
    state[3] ^= state[1];
    state[1] ^= state[2];
    state[0] ^= state[3];

    state[2] ^= t;
    state[3] = rol64(state[3], 45);

    return result;
}
