// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University.

#include "cutter.h"
#include "cutter_prop.h"
#include "glossary.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul_part.h"
#include "hand_monitor.h"
#include "space.h"
#include "smath.h"

//------------------------------------------------------------------------------

Cutter::Cutter(CutterProp const* p, HandMonitor* h)
: Hand(p,h)
{
    nextAct = RNG.exponential();
}


void Cutter::cut()
{
    assert_true( attached() );
    Fiber * fib = modifiableFiber();

    if ( prop()->selective == CUT_TOP_FIBER )
    {
        Hand const* h = otherHand();
        if ( h && h->attached() )
        {
            Vector P = pos();
            Vector Q = h->pos();
            Fiber const* fox = h->fiber();
            FiberProp const* FP = fib->prop;
#if ( 0 )
            // This is old code: using Z works only for certain geometries
            real PZ = P.ZZ;
            real QZ = Q.ZZ;
            std::clog << "Z ";
#else
            Space const* spc = FP->confine_space;
            Vector mid = 0.5 * ( P + Q );
            Vector prj = spc->project(mid);  // on the edge
            Vector dir = spc->normalToEdge(mid); // directed outward
            // calculate distances to the edge:
            real PZ = dot(prj-P, dir);
            real QZ = dot(prj-Q, dir);
#endif
            const real cutoff = FP->steric_radius;
            // do not sever the fiber that is closest to the edge:
            if ( PZ < QZ + cutoff )
                return;
            std::clog << "selective cut " << fib->reference() << " crossing at ";
            std::clog << Cymath::repr(PZ,8,4) << " over " << fox->reference();
            std::clog << " at " << Cymath::repr(QZ,8,4) << "\n";
        }
    }
    /**
     Cutting the fiber can invalidate the FiberGrid used for attachment,
     and this becomes a problem if the Cutter is part of a Couple,
     because calls for attachments and actions are intermingled.
     
     This is why severSoon() below will register the position of the cut,
     but the cut will only be performed later in Fiber::step()
     */
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
        if ( RNG.test_not(prop()->hold_growing_end[1]) )
            return detach();
        a = M;
    }
    
    if ( a >= P )
    {
        if ( RNG.test_not(prop()->hold_growing_end[0]) )
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
        if ( RNG.test_not(prop()->hold_growing_end[1]) )
            return detach();
        a = M;
    }
    
    if ( a >= P )
    {
        if ( RNG.test_not(prop()->hold_growing_end[0]) )
            return detach();
        a = P;
    }

    moveTo(a);
}

