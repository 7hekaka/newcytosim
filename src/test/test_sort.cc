// Cytosim was created by Francois Nedelec. Copyright 2024 Cambridge University.

#include <cstdio>
#include <stdlib.h>
#include <cstdint>
#include <algorithm>
#include <vector>

#include "timer.h"
#include "random_pcg.h"
#include "random_seed.cc"
using namespace PCG32;

constexpr size_t LOAD = 128;

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

// plain data type
struct stuff
{
    double X, Y, Z;
    uint32_t load[LOAD];
    stuff& operator = (uint32_t const i) { Z = i; return *this; }
    bool operator < (const stuff& s) const { return Z < s.Z; }
};

int compare_stuff(const void * A, const void * B)
{
    uint32_t a = static_cast<stuff const*>(A)->Z;
    uint32_t b = static_cast<stuff const*>(B)->Z;
    return ( a > b ) - ( a < b );
}


// indirect data type with std::vector
class vertex
{
public:
    std::vector<double> vec; // X, Y, Z
    uint32_t load[LOAD];
    
    vertex() : vec(3) { vec.resize(3,0.0); }
    vertex& operator = (uint32_t const i) { vec[2] = i; return *this; }
    bool operator < (const vertex& s) const { return vec[2] < s.vec[2]; }
};

int compare_vertex(const void * A, const void * B)
{
    double a = static_cast<vertex const*>(A)->vec[2];
    double b = static_cast<vertex const*>(B)->vec[2];
    return ( a > b ) - ( a < b );
}

// even more indirect data type with std::vector
class triangle
{
public:
    std::vector<std::vector<double> > pts; // A, B, C
    uint32_t load[LOAD];
    
    triangle() : pts(3) { pts[0].resize(3, 0.0); pts[1].resize(3, 0.0); pts[2].resize(3, 0.0);  }
    triangle& operator = (uint32_t const i) { pts[0][2] = i; return *this; }
    bool operator < (const triangle& s) const { return pts[0][2] < s.pts[0][2]; }
};

int compare_triangle(const void * A, const void * B)
{
    double a = static_cast<triangle const*>(A)->pts[0][2];
    double b = static_cast<triangle const*>(B)->pts[0][2];
    return ( a > b ) - ( a < b );
}




template < typename TYPE, int FUNC(void const*, void const*) >
void speed_test(const char str[], size_t cnt, size_t rep)
{
    size_t S = cnt * sizeof(TYPE);
    TYPE * val = new TYPE[cnt];

    // check random number generator
    tick();
    for ( size_t r = 0; r < rep; ++r )
    {
        for ( size_t i = 0; i < cnt; ++i )
            val[i] = pcg32();
    }
    double t0 = 1000*tock(cnt);

    // test C quicksort
    tick();
    for ( size_t r = 0; r < rep; ++r )
    {
        for ( size_t i = 0; i < cnt; ++i )
            val[i] = pcg32();
        qsort(val, cnt, sizeof(TYPE), FUNC);
    }
    double t1 = 1000*tock(cnt);
    
    // test C++ introsort with a Lambda function
    tick();
    for ( size_t r = 0; r < rep; ++r )
    {
        for ( size_t i = 0; i < cnt; ++i )
            val[i] = pcg32();
        std::sort(val, val+cnt, [](TYPE const& a, TYPE const& b) { return a < b; });
    }
    double t2 = 1000*tock(cnt);
    printf("%-10s   %8lu kB  set %8.2f  qsort %8.2f  std::sort %8.2f\n", str, S>>10, t0, t1, t2);
    delete[] val;
}


int main(int argc, char* argv[])
{
    size_t cnt = 1<<16, rep = 32;
    if ( argc > 1 )
        cnt = std::max(2, atoi(argv[1]));
    seed_pcg32();
    speed_test<uint16_t, compare16>("int16", cnt, rep);
    speed_test<uint32_t, compare32>("int32", cnt, rep);
    speed_test<uint64_t, compare64>("int64", cnt, rep);
    speed_test<stuff, compare_stuff>("stuff", cnt, rep);
    speed_test<vertex, compare_vertex>("vertex", cnt, rep);
    speed_test<triangle, compare_triangle>("triangle", cnt, rep);
    printf("done\n");
}
