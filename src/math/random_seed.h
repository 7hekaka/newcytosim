// Cytosim was created by Francois Nedelec. Copyright 2024 Cambridge University.

#include <cstdint>

namespace PCG32
{
    /**
     Find random seed from a random device, or from the system's clock
     */
    uint32_t get_random_seed();

}
