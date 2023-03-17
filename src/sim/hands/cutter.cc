// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "cutter.h"
#include "cutter_prop.h"
#include "glossary.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul_part.h"
#include "hand_monitor.h"

//------------------------------------------------------------------------------

Cutter::Cutter(CutterProp const* p, HandMonitor* h)
: Hand(p,h)
{
    nextAct = RNG.exponential();
}


void Cutter::cut()
{
    assert_true( attached() );

    /**
     Cutting the fiber can invalidate the FiberGrid used for attachment,
     and this becomes a problem if the Cutter is part of a Couple,
     because calls for attachments and actions are intermingled.
     
     This is why sever() below will register the position of the cut,
     but his cut will be performed later in Fiber::step()
     */
    fiber()->sever(abscissa(), prop()->new_end_state[0], prop()->new_end_state[1]);
    //std::clog << "cut " << fiber()->reference() << " @ " << abscissaFromM() << "\n";
    
    // simplest is to detach since the location will be at the tip of new fiber
    detach();
}

//------------------------------------------------------------------------------

void Cutter::stepUnloaded()
{
    assert_true( attached() );

    nextAct -= prop()->cutting_rate_dt;
    
    if ( nextAct < 0 )
    {
        nextAct = RNG.exponential();
        if ( abscissaFromM() < prop()->cutting_range )
            return cut();
    }
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
}

