// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
// test 3 methods to generate a random event time, when the rate varies in time
// FJN, Oct 2005 - 5 Feb. 2021

#include "random.h"
#include <cstdio>
#include "timer.h"

const size_t maxTime = 512;
const size_t SAMPLES = 1<<18;

real rates[maxTime] = { 0 };


/// this is the standard method: 64s CPU
size_t method1()
{
    for ( size_t i=0; i<maxTime; ++i )
    {
        if (RNG.test(rates[i])) return i;
    }
    return maxTime;
}

/// this is more exact and very slow: 370s CPU (an exponential at each step!)
size_t method2()
{
    for ( size_t i=0; i<maxTime; ++i )
    {
        if ( RNG.preal() < -std::expm1(-rates[i]) )
            return i;
    }
    return maxTime;
}

/// this is exact, and the fastest method: 10s CPU!
size_t method3()
{
    real T = -std::log( RNG.preal() );
    for ( size_t i=0; i<maxTime; ++i )
    {
        T -= rates[i];
        if ( T < 0 ) return i;
    }
    return maxTime;
}


template< size_t (*FUNC)(void) >
void testMethod(const char str[])
{
    unsigned bins[maxTime+1] = { 0 };
    
    tick();
    for ( size_t i = 0; i < SAMPLES; ++i )
        ++bins[FUNC()];
    printf("%s   %5.2f\n", str, tock(SAMPLES));

    if ( *str )
    {
        FILE* file = fopen(str, "w");
        for ( size_t i = 0; i <= maxTime; ++i )
            fprintf(file, "%4lu %6u\n", i, bins[i]);
        fclose(file);
    }
}


int main(int argc, char* argv[])
{
    RNG.seed();

    real rate = 0.005;
    if ( argc > 1 )
        rate = strtod(argv[1], 0);

    for ( size_t i = 0; i < maxTime; ++i )
        rates[i] = ( i % 10 ) * rate;

    testMethod<method1>("test1.txt");
    testMethod<method2>("test2.txt");
    testMethod<method3>("test3.txt");
}

