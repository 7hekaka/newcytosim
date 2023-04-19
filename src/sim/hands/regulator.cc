// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "regulator.h"
#include "regulator_prop.h"
#include "glossary.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul_part.h"
#include "hand_monitor.h"


Regulator::Regulator(RegulatorProp const* p, HandMonitor* h)
: Hand(p,h)
{
    //throw InvalidParameter("the Regulator class is unfinished");
}


    
void Regulator::attach(FiberSite const& s)
{
    Hand::attach(s);
    // freeze the plus end:
    Fiber * fib = modifiableFiber();
    fib->setEndStateP(STATE_WHITE);
}


void Regulator::stepUnloaded()
{
    assert_true( attached() );
}


void Regulator::stepLoaded(Vector const& force)
{
    assert_true( attached() );
}

