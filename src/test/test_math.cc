// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/*
 A test for the Floating Point Exceptions (Signal)
*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <unistd.h>
#include <csignal>
#include <cmath>

/*
 icpc --help
 
 -fp-trap=<arg>[,<arg>,...]
 control floating point traps at program start.  <arg> can be of the
 following values
 [no]divzero   - [Do not] trap on division by zero
 [no]inexact   - [Do not] trap on inexact result
 [no]invalid   - [Do not] trap on invalid operation
 [no]overflow  - [Do not] trap on overflow
 [no]underflow - [Do not] trap on underflow
 [no]denormal  - [Do not] trap on denormal
 all           - enable trap on all of the above
 none          - trap on none of the above
 common        - trap on most commonly used IEEE traps
 (invalid, division by zero, overflow)
 -fp-trap-all=<arg>[,<arg>,...]
 control floating point traps in every routine.  <arg> can be of the
 values specified in -fp-trap
 */

std::ostream& out = std::cout;

#include "real.h"

void signal_handler(int sig)
{
    psignal(sig, "test");
    _exit(sig);
}

void modulo()
{
    out << "   x    fmod remainder";
    for ( real x = -4; x <= 4; x += 0.5 )
    {
        out << "\n" << std::setw(5) << x;
        out << "  " << std::setw(5) << fmod(x, 2.0);
        out << "  " << std::setw(5) << remainder(x, 2.0);
    }
    out << '\n';
}

void divide()
{
    std::div_t dv;
    
    out << "   x    div.quot div.rem";
    for ( int x = -6; x <= 6; x += 1 )
    {
        dv = std::div(x, 3);
        out << "\n" << std::setw(5) << x;
        out << "  " << std::setw(5) << dv.quot;
        out << "  " << std::setw(5) << dv.rem;
    }
    out << '\n';
}

void std_copysign()
{
    out << "copysign(1, +1) = " << std::copysign(1.0, +1.0) << '\n';
    out << "copysign(1, -1) = " << std::copysign(1.0, -1.0) << '\n';
    out << "copysign(1,  0) = " << std::copysign(1.0,  0.0) << '\n';
    out << "copysign(1, -0) = " << std::copysign(1.0, -0.0) << '\n';
}

void infinities()
{
    real i = INFINITY;
    real z = 0;
    real y = 0 / z;
    real x = 1 / z;
    out << " -inf = " << -i << '\n';
    out << " 1.0/0.0 = " << x << '\n';
    out << " 0.0/0.0 = " << y << '\n';
    out << " 0   < inf = " << ( 0 < i ) << '\n';
    out << " inf < inf = " << ( i < i ) << '\n';
    out << " std::isinf(inf) = " << std::isinf(i) << '\n';
    out << " std::min(inf, 1) = " << std::min(i, 1.) << '\n';
    out << " std::max(inf, 1) = " << std::max(i, 1.) << '\n';
    out << " std::min(-inf, 1) = " << std::min(-i, 1.) << '\n';
    out << " std::max(-inf, 1) = " << std::max(-i, 1.) << '\n';
}

void not_numbers()
{
    real x = nan("");
    out << " nan + 1 = " << ( x + 1. ) << '\n';
    out << " nan < 1 = " << ( x < 1. ) << '\n';
    out << " nan > 0 = " << ( x > 0. ) << '\n';
    out << " std::min(nan, 1) = " << std::min(x, 1.) << '\n';
    out << " std::max(nan, 1) = " << std::max(x, 1.) << '\n';
    out << " abs_real(nan) = " << abs_real(x) << '\n';
    out << " std::ceil(nan) = " << std::ceil(x) << '\n';
    out << " std::floor(nan) = " << std::floor(x) << '\n';
}

void numbers()
{
    out << " 1.0 / 0 = " <<  1.0 / 0 << '\n';
    out << "-1.0 / 0 = " << -1.0 / 0 << '\n';
    out << " 0.0 / 0 = " <<  0.0 / 0 << '\n';
    out << "-log(0)  = " << -std::log(0.0) << '\n';
#if ( 1 )
    out << "absf(-2) = " << std::fabs(-2.0) << '\n';
    out << "absf(-1) = " << std::fabs(-1.) << '\n';
    out << "absf(+1) = " << std::fabs(+1.) << '\n';
    out << "absf(+2) = " << std::fabs(+2.) << '\n';
#endif
}


void read_numbers(std::string const& str)
{
    std::stringstream is(str);
    double x, y, z;
    if ( is >> x )
        out << "read x = " << x << '\n';
    if ( is >> y )
        out << "read y = " << y << '\n';
    if ( is >> z )
        out << "read z = " << z << '\n';
    
    char tmp[128] = { 0 };
    is.clear();
    is.readsome(tmp, sizeof(tmp));
    out << "leftover `" << tmp << "'\n";
}


int main(int argc, char* argv[])
{
    if ( argc > 1 )
    {
        read_numbers(argv[1]);
        return 0;
    }
    if ( signal(SIGFPE, signal_handler) == SIG_ERR )
    {
        out << "Could not register SIGFPE handler\n";
        return EXIT_FAILURE;
    }
    divide();
    modulo();
    std_copysign();
    infinities();
    not_numbers();
    numbers();
    out << "test completed" << '\n';
}
