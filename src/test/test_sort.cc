// Cytosim was created by Francois Nedelec. Copyright 2024 Cambridge University.

#include <cstdio>
#include <stdlib.h>
#include <cstdint>
#include <algorithm>

#include "timer.h"
#include "random_pcg.h"
#include "random_seed.cc"


uint64_t pcg32_state;


int compare16(const void * A, const void * B)
{
    uint16_t a = *static_cast<uint16_t const*>(A);
    uint16_t b = *static_cast<uint16_t const*>(B);
    return ( a > b ) - ( a < b );
}

int compare32(const void * A, const void * B)
{
    uint32_t a = *static_cast<uint32_t const*>(A);
    uint32_t b = *static_cast<uint32_t const*>(B);
    return ( a > b ) - ( a < b );
}

int compare64(const void * A, const void * B)
{
    uint64_t a = *static_cast<uint64_t const*>(A);
    uint64_t b = *static_cast<uint64_t const*>(B);
    return ( a > b ) - ( a < b );
}

struct stuff
{
    uint32_t key;
    uint32_t load[511];
    stuff& operator = (uint32_t const i) { key = i; return *this; }
    bool operator < (const stuff s) const { return key < s.key; }
};

int compare(const void * A, const void * B)
{
    uint32_t a = static_cast<stuff const*>(A)->key;
    uint32_t b = static_cast<stuff const*>(B)->key;
    return ( a > b ) - ( a < b );
}


template < typename INTEGER, int COMPARE(void const*, void const*) >
void speed_test(size_t cnt, size_t rep)
{
    size_t S = cnt * sizeof(INTEGER);
    INTEGER * val = (INTEGER*)malloc(S);

    // check random number generator
    tick();
    for ( size_t r = 0; r < rep; ++r )
    {
        for ( size_t i = 0; i < cnt; ++i )
            val[i] = pcg32(pcg32_state);
    }
    double t0 = 1000*tock(cnt);

    // test C quicksort
    tick();
    for ( size_t r = 0; r < rep; ++r )
    {
        for ( size_t i = 0; i < cnt; ++i )
            val[i] = pcg32(pcg32_state);
        qsort(val, cnt, sizeof(INTEGER), COMPARE);
    }
    double t1 = 1000*tock(cnt);
    
    // test C++ introsort with a Lambda function
    tick();
    for ( size_t r = 0; r < rep; ++r )
    {
        for ( size_t i = 0; i < cnt; ++i )
            val[i] = pcg32(pcg32_state);
        std::sort(val, val+cnt, [](INTEGER const& a, INTEGER const& b) { return a < b; });
    }
    double t2 = 1000*tock(cnt);
    printf("int%-6lu  size %8lu kB  set %8.2f  qsort %8.2f  std::sort %8.2f\n", 8*sizeof(INTEGER), S>>10, t0, t1, t2);
    free(val);
}


int main(int argc, char* argv[])
{
    size_t cnt = 1<<16, rep = 256;
    if ( argc > 1 )
        cnt = std::max(2, atoi(argv[1]));
    pcg32_state = get_random_seed();
    speed_test<uint16_t, compare16>(cnt, rep);
    speed_test<uint32_t, compare32>(cnt, rep);
    speed_test<uint64_t, compare64>(cnt, rep);
    speed_test<stuff, compare>(cnt, rep);
    printf("done\n");
}
