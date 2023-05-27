// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University.

#include "cutter.h"
#include "cutter_prop.h"
#include "glossary.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul_part.h"
#include "hand_monitor.h"
#include "space.h"

//------------------------------------------------------------------------------

Cutter::Cutter(CutterProp const* p, HandMonitor* h)
: Hand(p,h)
{
    nextAct = RNG.exponential();
}


void Cutter::cut()
{
    assert_true( attached() );

    if ( prop()->selective == CUT_TOP_FIBER )
    {
        Hand const* h = otherHand();
        if ( !h || !h->attached() )
            return;
        Vector P = pos();
        Vector Q = h->pos();
#if ( 0 )
        // using Z works only if the boundary is aligned in XY with cell on top
        real PZ = P.ZZ;
        real QZ = Q.ZZ;
#else
        Vector PQ = 0.5 * ( P + Q );
        Space const* spc = h->fiber()->prop->confine_space;
        Vector prj = spc->project(PQ);  // on the edge
        Vector dir = spc->normalToEdge(PQ);
        // calculate distances to the edge:
        real PZ = dot(prj-P, dir);
        real QZ = dot(prj-Q, dir);
#endif
        // do not sever the fiber that is closest to the edge:
        if ( PZ < QZ )
            return;
        std::clog << "selectively cutting crossing fiber at Z " << PZ << " over " << QZ << "\n";
    }
    /**
     Cutting the fiber can invalidate the FiberGrid used for attachment,
     and this becomes a problem if the Cutter is part of a Couple,
     because calls for attachments and actions are intermingled.
     
     This is why severSoon() below will register the position of the cut,
     but the cut will only be performed later in Fiber::step()
     */
    Fiber * fib = modifiableFiber();
    fib->severSoon(abscissa(), prop()->cut_width, prop()->new_end_state[0], prop()->new_end_state[1]);
    //std::clog << "cut " << fiber()->reference() << " @ " << abscissa() << "\n";
    
    // simplest is to detach since the location will be at the tip of new fiber
    detach();
}

//------------------------------------------------------------------------------

void Cutter::stepUnloaded()
{
    assert_true( attached() );
    
    if ( prop()->selective == CUT_ANY_FIBER )
    {
        nextAct -= prop()->cutting_rate_dt;
        
        if ( nextAct < 0 )
        {
            nextAct = RNG.exponential();
            if ( abscissaFromM() < prop()->cutting_range )
                return cut();
        }
    }
    
    real a = hAbs + prop()->line_diffusion_dt * RNG.sreal();
    
    const real M = hFiber->abscissaM();
    const real P = hFiber->abscissaP();
    
    if ( a <= M )
    {
        if ( RNG.test_not(prop()->hold_growing_end) )
            return detach();
        a = M;
    }
    
    if ( a >= P )
    {
        if ( RNG.test_not(prop()->hold_growing_end) )
            return detach();
        a = P;
    }
    
    moveTo(a);
}


void Cutter::stepLoaded(Vector const& force)
{
    assert_true( attached() );
    
    nextAct -= prop()->cutting_rate_dt;
    
    if ( nextAct < 0 )
    {
        nextAct = RNG.exponential();
        if ( abscissaFromM() < prop()->cutting_range )
            return cut();
    }
    
    real load = dot(force, dirFiber());
    real diff = prop()->line_diffusion_dt * RNG.sreal();
    real a = hAbs + diff + prop()->movability_dt * load;
    
    const real M = hFiber->abscissaM();
    const real P = hFiber->abscissaP();
    
    if ( a <= M )
    {
        if ( RNG.test_not(prop()->hold_growing_end) )
            return detach();
        a = M;
    }
    
    if ( a >= P )
    {
        if ( RNG.test_not(prop()->hold_growing_end) )
            return detach();
        a = P;
    }

    moveTo(a);
}

