// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "real.h"
#include "random.h"
#include <cstdio>
#include <bitset>
#include <cstring>
#include <iostream>
#include <climits>
#include "timer.h"
#include "vector.h"
#include "random_vector.h"

#include <random>

template < typename T >
void print_bits(FILE* f, const T& val, char spc)
{
    unsigned char * ptr = (unsigned char*) & val;
    for ( int i = sizeof(T)-1; i >= 0; --i)
    {
        unsigned char byte = ptr[i];
        for ( int j = 0; j < CHAR_BIT; ++j )
        {
            putc('0' + (1 & (byte>>(CHAR_BIT-1))), f);
            byte <<= 1;
        }
        if ( spc ) putc(spc, f);
    }
    putc('\n', f);
}


void test_int()
{
    for (int j=0; j<8; ++j)
    {
        for (int k=0; k<8; ++k)
            printf(" %12u", RNG.pint32());
        printf("\n");
    }
    printf("\n");
    
    for (int j=0; j<8; ++j)
    {
        for (int k=0; k<8; ++k)
            printf(" %+12i", RNG.sint32());
        printf("\n");
    }
    printf("\n");

    for (int j=0; j<8; ++j)
    {
        for (int k=0; k<32; ++k)
            printf(" %2u", RNG.pint32(100));
        printf("\n");
    }
    printf("\n");
    
    for (int j=0; j<8; ++j)
    {
        for (int k=0; k<32; ++k)
            printf(" %2u", RNG.pint32_fair(100));
        printf("\n");
    }
    printf("\n");

    for (int j=0; j<8; ++j)
    {
        for (int k=0; k<32; ++k)
            printf(" %2u", RNG.pint32_slow(99));
        printf("\n");
    }
    printf("\n");
}


void silly_test()
{
    const uint32_t up = 1 << 30;
    
    const uint32_t cnt = 1 << 24;
    uint32_t hit = 0;
    
    for (uint32_t j=0; j<cnt; ++j)
        hit += ( RNG.pint32() < up );

    printf(" prob( pint32() < 1^30 ) = %f\n", hit/(float)cnt);
}


/**
 This assumes IEEE Standard 754 Floating point numbers
32 bits: 1 for sign, 8 for exponents, 23 for fraction
 */
float convertFix(uint32_t x)
{
    constexpr uint32_t FRAC  = 0x7FFFFFU;
    constexpr uint32_t EXPON = 127 << 23;
    uint32_t res = EXPON | ( x & FRAC );
    return *((float*)&res) - 1.0;
}


void testbits()
{
    const int SCALE = 2;
    for ( int i=0; i <= SCALE; ++i )
    {
        float x = i / float(SCALE);
        printf(" %f :", x);
        print_bits(stdout, x, 0);
        // x = -ii / float(SCALE);
        // printf("%f :", x);
        // print_bits(stdout, x, 0);
    }
    
    for ( int i=0; i < 16; ++i )
    {
        float y = convertFix(RNG.pint32());
        printf(" %f :", y);
        print_bits(stdout, y, ' ');
    }
}


#define TEST test
void check_test( const real prob, const size_t MAX )
{
    int cnt = 0, a, b, c;
    for ( size_t jj=0; jj < MAX; ++jj )
    {
        a = RNG.TEST(prob) + RNG.TEST(prob) + RNG.TEST(prob) + RNG.TEST(prob);
        b = RNG.TEST(prob) + RNG.TEST(prob) + RNG.TEST(prob) + RNG.TEST(prob);
        c = RNG.TEST(prob) + RNG.TEST(prob) + RNG.TEST(prob) + RNG.TEST(prob);
        cnt += a + b + c;
    }
    printf("prob = %f measured = %f cnt = %i\n", prob, cnt / double(12*MAX), cnt);
}

void test_RNG(const size_t MAX)
{
    for ( size_t jj=0; jj < MAX; ++jj )
    {
        RNG.preal();RNG.preal();RNG.preal();RNG.preal();RNG.preal();
        RNG.preal();RNG.preal();RNG.preal();RNG.preal();RNG.preal();
    }
}


void test_real()
{
    for (int j=0; j<8; ++j)
    {
        for (int k=0; k<8; ++k)
            printf(" %10f", RNG.sreal());
        printf("\n");
    }

    printf("\n");
    for (int j=0; j<8; ++j)
    {
        for (int k=0; k<8; ++k)
            printf(" %10f", RNG.preal());
        printf("\n");
    }
    
    printf("\n");
    for (int j=0; j<8; ++j)
    {
        for (int k=0; k<8; ++k)
            printf(" %10f", RNG.shalf());
        printf("\n");
    }

    printf("\npfloat:     ");
    float x;
    for ( int kk=0; kk < 10; ++kk )
    {
        x = RNG.pfloat();
        printf(" %+f", x);
    }

    printf("\nsfloat:     ");
    for ( int kk=0; kk < 10; ++kk )
    {
        x = RNG.sfloat();
        printf(" %+f", x);
    }
    
    double d;
    printf("\npdouble:    ");
    for ( int kk=0; kk < 10; ++kk )
    {
        d = RNG.pdouble();
        printf(" %+f", d);
    }

    printf("\nsdouble:    ");
    for ( int kk=0; kk < 10; ++kk )
    {
        d = RNG.sdouble();
        printf(" %+f", d);
    }

    printf("\nsflip:      ");
    for ( int kk=0; kk < 10; ++kk )
    {
        d = RNG.sflip();
        printf(" %+f", d);
    }
    printf("\n");
}

//==============================================================================

void check_uniform(size_t cnt)
{
    const double off = 0.5;
    double avg = 0, var = 0;
    for ( size_t i = 0; i < cnt; ++i )
    {
        real x = RNG.sreal() - off;
        real y = RNG.sreal() - off;
        real z = RNG.sreal() - off;
        real t = RNG.sreal() - off;
        avg += x + y + z + t;
        var += x*x + y*y + z*z + t*t;
    }
    cnt *= 4;
    if ( cnt > 0 )
    {
        avg /= cnt;
        var -= square(avg)*cnt;
    }
    if ( cnt > 1 )
        var /= real(cnt-1);
    printf("UNIFORM      avg = %.12e   var = %.12e\n", avg+off, var);
}


void check_gauss(size_t CNT)
{
    size_t cnt = 0;
    double avg = 0, var = 0;
    const size_t n_max = 1<<6;
    real vec[n_max] = { 0 };
    for ( size_t i = 0; i < CNT; ++i )
    {
        size_t n = RNG.pint32(n_max);
        RNG.gauss_set(vec, n);
        cnt += n;
        for ( size_t u = 0; u < n; ++u )
        {
            avg += vec[u];
            var += vec[u] * vec[u];
        }
    }
    if ( cnt > 0 )
    {
        avg /= cnt;
        var -= square(avg)*cnt;
    }
    if ( cnt > 1 )
        var /= real(cnt-1);
    printf("GAUSSIAN     avg = %.12e   var = %.12e\n", avg, var);
}


void check_prob()
{
    size_t avg = 0;
    size_t cnt = 1 << 28;
    for ( size_t i = 0; i < cnt; ++i )
        avg += RNG.flip_8th();

    printf("8th      prob = %.6f\n", avg/(double)cnt);
}


void check_exponential(size_t cnt)
{
    double ix = INFINITY, iy = INFINITY, iz = INFINITY, it = INFINITY;
    const double off = 1;
    double avg = 0, var = 0;
    for ( size_t i = 0; i < cnt; ++i )
    {
        real x = RNG.exponential() - off;
        real y = RNG.exponential() - off;
        real z = RNG.exponential() - off;
        real t = RNG.exponential() - off;
        ix = std::min(ix, x);
        iy = std::min(iy, y);
        iz = std::min(iz, z);
        it = std::min(it, t);
        avg += x + y + z + t;
        var += x*x + y*y + z*z + t*t;
    }
    ix = std::min(std::min(ix, iz), std::min(ix, iz));
    cnt *= 4;
    if ( cnt > 0 )
    {
        avg /= cnt;
        var -= square(avg)*cnt;
    }
    if ( cnt > 1 )
        var /= real(cnt-1);
    printf("EXPONENTIAL min: %.12e   avg = %.12e   var = %.12e\n", ix+off, avg+off, var);
}


void check_poisson(size_t sup)
{
    for ( size_t n = 0; n < sup; ++n )
    {
        int x = (int)(RNG.gauss() * std::sqrt(n) + n);
        printf("%10lu %9i %9i %9i\n", n, RNG.poisson_knuth(n), RNG.poisson(n), x);
    }
}

//==============================================================================

void speed_test(size_t cnt)
{
    tick();
    size_t c = 0;
    uint32_t u = 10;
    for ( size_t j = 0; j < cnt; ++j )
    {
        u = RNG.pint32(1024);
        c += RNG.pint32(u);
        c += RNG.pint32(u);
        
        c += RNG.pint32(u);
        c += RNG.pint32(u);
        c += RNG.pint32(u);
        
        c += RNG.pint32(u);
        c += RNG.pint32(u);
        c += RNG.pint32(u);
        
        c += RNG.pint32(u);
        c += RNG.pint32(u);
        c += RNG.pint32(u);
    }
    printf("3x integers %5.2f  ( %lu )\n", tock(cnt>>17), c);
}


void speed_test_vector(size_t cnt)
{
    tick();
    Vector3 sum(0,0,0);
    for ( size_t j = 0; j < cnt; ++j )
    {
        sum += Vector3::randU();
        sum += Vector3::randU();
        sum += Vector3::randU();
        sum += Vector3::randU();
    }
    printf("Unit vector %5.2f  ", tock(cnt>>17));
    printf("( %5.2f %5.2f %5.2f )\n", sum.XX, sum.YY, sum.ZZ);
}


int main(int argc, char* argv[])
{
    int mode = 1;
    RNG.seed();

    if ( argc > 1 )
        mode = atoi(argv[1]);
    real rate = 1;
    if ( argc > 2 )
        rate = strtod(argv[2], nullptr);

    switch ( mode )
    {
        case 0:
            check_poisson(1024);
            check_prob();
            break;
            
        case 1:
            check_exponential(1<<26);
            //check_uniform(1<<20);
            //check_gauss(1<<20);
            break;

        case 2:
            for ( int kk=0; kk < 11; ++kk )
                check_test(rate*kk, 5000000);
            break;
            
        case 3:
            printf("sizeof(uint32_t) = %lu\n", sizeof(uint32_t));
            test_int();
            test_real();
            break;
            
        case 4:
            speed_test(1048576);
            speed_test_vector(1048576);
            break;
            
        case 5:
            silly_test();
            break;
    }
}

