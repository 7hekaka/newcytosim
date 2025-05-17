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
    delay = INFINITY;
    recurrent = true;
    multiplexed = true;
    nextTime = 0;
    stop = INFINITY;
}


/** This will make the event fire once, after 'time' */
void Event::fire_once_at(double time)
{
    nextTime = time;
    recurrent = false;
    multiplexed = false;
    delay = INFINITY;
}


/** This will make the event fire once at every time step, after 'time' */
void Event::fire_always_after(double time)
{
    nextTime = time;
    recurrent = true;
    multiplexed = false;
    delay = 0;
}


void Event::reload(double t)
{
    if ( recurrent )
    {
        if ( rate > 0 )
            nextTime = t + RNG.exponential() / rate;
        else
            nextTime = t + delay;
    }
    else
    {
        nextTime = INFINITY;
    }
}


Event::Event(double now, Glossary& opt)
{
    clear();
    opt.set(activity, "activity", "code");
    double t = now;

    if ( opt.set(t, "time") )
    {
        fire_once_at(t);
    }
    else if ( opt.set(rate, "rate") )
    {
        if ( rate < 0 )
            throw InvalidParameter("event:rate must be > 0");
        opt.set(multiplexed, "multiplexed");
        opt.set(t, "start");
        opt.set(stop, "stop");
        reload(t);
    }
    else if ( opt.set(delay, "interval") || opt.set(delay, "delay") )
    {
        if ( delay <= 0 )
            throw InvalidParameter("event:delay must be > 0");
        opt.set(multiplexed, "multiplexed");
        opt.set(t, "start");
        opt.set(stop, "stop");
        reload(t);
    }
    else
    {
        opt.set(t, "start");
        opt.set(stop, "stop");
        fire_always_after(t);
    }
}


Event::~Event()
{
    //Cytosim::log("destroying Event (", activity, ")\n");
}


/**
 This is called once per time step
 */
void Event::step(Simul& sim)
{
    if ( sim.time() >= nextTime )
    {
        if ( sim.time() > stop )
        {
            nextTime = INFINITY;
            return;
        }
        sim.relax();
        // if 'multiplexed', the event can fire multiple times within a time step
        do {
            try {
                //std::clog << "firing Event (" << activity << ")\n";
                sim.perform(activity);
            }
            catch( Exception & e ) {
                Cytosim::warn("Buggy event:code `", activity, "' : ", e.message(), "\n");
            }
            reload(nextTime);
        } while ( multiplexed && sim.time() >= nextTime );
        sim.unrelax();
    }
}


void Event::write(Outputter& out) const
{
}


void Event::read(Inputter& in, Simul& sim, ObjectTag tag)
{
}
