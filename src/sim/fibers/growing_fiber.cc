// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "assert_macro.h"
#include "growing_fiber.h"
#include "growing_fiber_prop.h"
#include "object_set.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul_part.h"
#include "couple.h"
#include "space.h"
#include "simul.h"


//------------------------------------------------------------------------------

GrowingFiber::GrowingFiber(GrowingFiberProp const* p) : Fiber(p)
{
    mStateM = STATE_GREEN;
    mStateP = STATE_GREEN;
}


GrowingFiber::~GrowingFiber()
{
}


//------------------------------------------------------------------------------
#pragma mark -


void GrowingFiber::setEndStateM(state_t s)
{
    if ( s == STATE_WHITE || s == STATE_GREEN )
        mStateM = s;
    else
        throw InvalidParameter("invalid AssemblyState ("+std::to_string(s)+") for GrowingFiber MINUS end");
}


void GrowingFiber::setEndStateP(state_t s)
{
    if ( s == STATE_WHITE || s == STATE_GREEN )
        mStateP = s;
    else
        throw InvalidParameter("invalid AssemblyState ("+std::to_string(s)+") for GrowingFiber PLUS end");
}

//------------------------------------------------------------------------------

void GrowingFiber::step()
{
    constexpr size_t P = 0, M = 1;
    real addP = 0, addM = 0;
    
    // plus end
    if ( prop()->shrink_outside[P] && prop()->confine_space->outside(posEndP()) )
    {
        addP = prop()->shrinking_speed_dt[P];
    }
    else if ( mStateP == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real forceP = projectedForceEndP();
        
        // growth is reduced if free monomers are scarce:
        addP = prop()->growing_speed_dt[P] * prop()->free_polymer;
        
        // antagonistic force (< 0) decreases assembly rate exponentially
        if (( forceP < 0 ) & ( addP > 0 ))
            addP *= std::exp(forceP*prop()->growing_force_inv[P]);
        
        addP += prop()->growing_off_speed_dt[P];
    }
    
    // minus end
    if ( prop()->shrink_outside[M] && prop()->confine_space->outside(posEndM()) )
    {
        addM = prop()->shrinking_speed_dt[M];
    }
    else if ( mStateM == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real forceM = projectedForceEndM();
        
        // growth is reduced if free monomers are scarce:
        addM = prop()->growing_speed_dt[M] * prop()->free_polymer;
        
        // antagonistic force (< 0) decreases assembly rate exponentially
        if (( forceM < 0 ) & ( addM > 0 ))
            addM *= std::exp(forceM*prop()->growing_force_inv[M]);

        addM += prop()->growing_off_speed_dt[M];
    }

    if ( Fiber::updateLength(addM, addP) )
    {
        Fiber::step();
        
        if ( length() > prop()->divide )
        {
            real L = 0.5 * length();
            Fiber * frag = severSegment(L, L);
            objset()->add(frag);
            if ( prop()->divide_couple )
            {
                Couple * C = prop()->divide_couple->newCouple();
                C->attachEnd1(this, PLUS_END);
                C->attachEnd2(frag, MINUS_END);
                simul().add(C);
            }
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark -


void GrowingFiber::write(Outputter& out) const
{
    Fiber::write(out);

    // since states are constant, we write growth rates:
    writeMarker(out, DYNAMIC_TAG);
    if ( mStateM == STATE_GREEN )
        out.writeFloat(cDeltaM);
    else
        out.writeFloat(-0.0);
    if ( mStateP == STATE_GREEN )
        out.writeFloat(cDeltaP);
    else
        out.writeFloat(-0.0);
}


void GrowingFiber::readEndStates(Inputter& in)
{
#if BACKWARD_COMPATIBILITY < 54
    if ( in.formatID() < 54 )
    {
        cDeltaM = in.readFloat();
        if ( in.formatID() > 45 )
            cDeltaP = in.readFloat();
    }
    else if ( in.formatID() < 56 )
    {
        cDeltaM = in.readFloat();
        in.readFloat();
        cDeltaP = in.readFloat();
        in.readFloat();
    }
    else
#endif
    {
        cDeltaM = in.readFloat();
        cDeltaP = in.readFloat();
        // the state are derived from the sign bit:
        mStateM = std::signbit(cDeltaM) ? STATE_WHITE : STATE_GREEN;
        mStateP = std::signbit(cDeltaP) ? STATE_WHITE : STATE_GREEN;
    }
}


void GrowingFiber::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    if ( tag == DYNAMIC_TAG )
        readEndStates(in);
    else
    {
#if BACKWARD_COMPATIBILITY < 44
        if ( tag == TAG && in.formatID() < 44 )
            readEndStates(in);
#endif
#if BACKWARD_COMPATIBILITY < 46
        const real len = length();
#endif
        
        Fiber::read(in, sim, tag);
                
#if BACKWARD_COMPATIBILITY < 46
        if ( tag == TAG && in.formatID() < 46 )
        {
            // adjust growing variable
            cDeltaP = length() - len;
            cDeltaM = 0;
        }
#endif
    }
}

