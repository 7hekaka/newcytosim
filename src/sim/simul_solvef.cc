// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "simul.h"

/**
 The motion occurs along the X-axis for all fibers,
 and is directed parallel to the Fiber axis if ( flux_speed > 0 ).
 The horizontal speed is proportional to the X component of the fiber's direction.
 
 If ( flux_speed < 0 ), right pointing fibers move left, while Left pointing fiber move right.
 */
void Simul::solveF(real shift)
{
    for ( Fiber * fib = fibers.first(); fib ; fib=fib->next() )
    {
        real const* pos = fib->data();
        if ( pos[DIM] > pos[0] )
            fib->translate( Vector(-shift, 0, 0) );
        else
            fib->translate( Vector( shift, 0, 0) );
    }
}


