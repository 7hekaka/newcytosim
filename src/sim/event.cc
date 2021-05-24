// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include "event.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"


void Event::clear()
{
    activity = "";
    rate = 0;
    delay = 0;
    recurrent = true;
    nextTime = 0;
}


void Event::fire_once_at(double time)
{
    nextTime = time;
    recurrent = false;
}


void Event::reload(double now)
{
    if ( recurrent )
    {
        if ( rate > 0 )
            nextTime = now + RNG.exponential() / rate;
        else
            nextTime = now + delay;
    }
    else
    {
        nextTime = INFINITY;
    }
}


Event::Event(double now, Glossary& opt)
{
    double t = now;
    clear();
    opt.set(activity, "activity", "code");
    
    if ( opt.set(t, "time") )
    {
        fire_once_at(t);
    }
    else if ( opt.set(rate, "rate") )
    {
        if ( rate <= 0 )
            throw InvalidParameter("event:rate must be > 0");
        reload(now);
    }
    else if ( opt.set(delay, "interval") || opt.set(delay, "delay") )
    {
        if ( delay <= 0 )
            throw InvalidParameter("event:delay must be > 0");
        reload(now);
    }
    else
        throw InvalidParameter("event:time, rate or delay must be specified");
}


Event::~Event()
{
    //Cytosim::log("destroying Event %p\n", this);
}


/**
 This is called once per time step
 */
void Event::step(Simul& sim)
{
    if ( sim.time() >= nextTime )
    {
        sim.relax();
        // the event can fire multiple time at each time step
        do {
            reload(nextTime);
            sim.evaluate(activity);
        } while ( sim.time() >= nextTime );
        sim.unrelax();
    }
}


void Event::write(Outputter& out) const
{
}


void Event::read(Inputter& in, Simul& sim, ObjectTag tag)
{
}
