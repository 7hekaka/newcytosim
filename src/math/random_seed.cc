// Cytosim was created by Francois Nedelec. Copyright 2024 Cambridge University.

#include "random_seed.h"

#include <limits.h>
#include <sys/time.h>
#include <iostream>
#include <random>

namespace PCG32
{
    /**
     Combine t and c to form a new integer
     Better than uint32_t(x) in case x is floating point in [0,1]
     Based on code by Lawrence Kirby (fred@genesis.demon.co.uk)
     */
    static uint32_t hash(long t, long c)
    {
        uint32_t h1 = 0;
        unsigned char* p = (unsigned char*) &t;
        for ( size_t i = 0; i < sizeof(t); ++i )
        {
            h1 *= UCHAR_MAX + 2U;
            h1 += p[i];
        }
        uint32_t h2 = 0;
        p = (unsigned char*) &c;
        for ( size_t j = 0; j < sizeof(c); ++j )
        {
            h2 *= UCHAR_MAX + 2U;
            h2 += p[j];
        }
        return h1 ^ h2;
    }
    
    
    /**
     Find random seed from a random device, or from the system's clock
     */
    uint32_t get_random_seed()
    {
        uint32_t s = 0;
        try {
            // read system source if available:
            std::random_device rd;
            s = rd();
            //std::cerr << "random_device ---> seed " << s << '\n';
        }
        catch (const std::exception &e) {
            std::cerr << e.what() << '\n';
            // use clock otherwise
            timespec now;
            clock_gettime(CLOCK_MONOTONIC, &now);
            s = hash(now.tv_sec, now.tv_nsec);
            //std::cerr << "system clock ---> seed " << s << '\n';
        }
        if ( s == 0 )
            return 7;
        return s;
    }
    
}
