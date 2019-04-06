// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "event.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"
#include "parser.h"

void Event::initialize(real time)
{
    stochastic = ( rate > 0 );
    nextEvent = time + RNG.exponential() / rate;
}


void Event::initialize(real time, Glossary& opt)
{
    if (!opt.set(rate, "rate")) rate = 0;
    if (!opt.set(code, "code")) code = "";
    initialize(time);
}


Event::~Event()
{
    //Cytosim::log("destroying Event %p\n", this);
}


/// stochastic firing at specified rate
void Event::step(Simul& sim)
{
    if ( sim.time() > nextEvent )
    {
        sim.relax();
        do {
            nextEvent += RNG.exponential() / rate;
            Parser(sim, 1, 1, 1, 1, 1).evaluate(code, ", in event:code");
        } while ( sim.time() > nextEvent );
        sim.prepare();
    }
}


void Event::write(Outputter& out) const
{
}


void Event::read(Inputter & in, Simul& sim, ObjectTag tag)
{
}
