// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "assert_macro.h"
#include "treadmilling_fiber.h"
#include "treadmilling_fiber_prop.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "space.h"


//------------------------------------------------------------------------------

TreadmillingFiber::TreadmillingFiber(TreadmillingFiberProp const* p) : Fiber(p)
{
    mStateM = STATE_WHITE;
    mStateP = STATE_WHITE;
}


TreadmillingFiber::~TreadmillingFiber()
{
}


//------------------------------------------------------------------------------
#pragma mark -

void TreadmillingFiber::setEndStateM(state_t s)
{
    if ( s == STATE_WHITE || s == STATE_GREEN || s == STATE_RED )
        mStateM = s;
    else
        throw InvalidParameter("invalid AssemblyState ("+std::to_string(s)+") for TreadmillingFiber minus end");
}


void TreadmillingFiber::setEndStateP(state_t s)
{
    if ( s == STATE_WHITE || s == STATE_GREEN || s == STATE_RED )
        mStateP = s;
    else
        throw InvalidParameter("invalid AssemblyState ("+std::to_string(s)+") for TreadmillingFiber plus end");
}

//------------------------------------------------------------------------------

void TreadmillingFiber::step()
{
    constexpr unsigned P = 0, M = 1;
    real addP = 0, addM = 0;

    if ( mStateP == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real forceP = projectedForceEndP();
        
        // growth is reduced if free monomers are scarce:
        addP = prop()->growing_speed_dt[P] * prop()->free_polymer;
        
        // antagonistic force (< 0) decreases assembly rate exponentially
        if (( forceP < 0 ) & ( addP > 0 ))
            addP *= std::exp(forceP*prop()->growing_force_inv[P]);
    }
    else if ( mStateP == STATE_RED )
    {
        addP = prop()->shrinking_speed_dt[P];
    }
    
    // minus end dynamics
    if ( mStateM == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real forceM = projectedForceEndM();
        
        // growth is reduced if free monomers are scarce:
        addM = prop()->growing_speed_dt[M] * prop()->free_polymer;
        
        // antagonistic force (< 0) decreases assembly rate exponentially
        if (( forceM < 0 ) & ( addM > 0 ))
            addM *= std::exp(forceM*prop()->growing_force_inv[M]);
    }
    else if ( mStateM == STATE_RED )
    {
        addM = prop()->shrinking_speed_dt[M];
    }

    if ( Fiber::updateLength(addM, addP) )
        Fiber::step();
}


//------------------------------------------------------------------------------
#pragma mark -


void TreadmillingFiber::write(Outputter& out) const
{
    Fiber::write(out);

    // write variables describing the dynamic state of the ends:
    writeMarker(out, DYNAMIC_TAG);
    out.writeUInt16(mStateM);
    out.writeUInt16(0);
    out.writeUInt16(mStateP);
    out.writeUInt16(0);
}


void TreadmillingFiber::readEndStates(Inputter& in)
{
#if BACKWARD_COMPATIBILITY < 54
    if ( in.formatID() < 54 )
    {
        mStateM = in.readUInt16();
        mStateP = in.readUInt16();
    }
    else
#endif
    {
        mStateM = in.readUInt16();
        in.readUInt16();
        mStateP = in.readUInt16();
        in.readUInt16();
    }
}


void TreadmillingFiber::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    //std::clog << " TreadmillingFiber::read(" << tag << ")\n";
    if ( tag == DYNAMIC_TAG )
        readEndStates(in);
    else
    {
#if BACKWARD_COMPATIBILITY < 44
        if ( tag == TAG && in.formatID() < 44 )
            readEndStates(in);
#endif
        Fiber::read(in, sim, tag);
    }
}

