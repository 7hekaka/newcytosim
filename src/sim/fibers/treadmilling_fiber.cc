// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "assert_macro.h"
#include "treadmilling_fiber.h"
#include "treadmilling_fiber_prop.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul.h"
#include "space.h"


//------------------------------------------------------------------------------

TreadmillingFiber::TreadmillingFiber(TreadmillingFiberProp const* p) : Fiber(p), prop(p)
{
    mStateM = STATE_WHITE;
    mStateP = STATE_WHITE;
    mGrowthM = 0;
    mGrowthP = 0;
}


TreadmillingFiber::~TreadmillingFiber()
{
    prop = nullptr;
}


//------------------------------------------------------------------------------
#pragma mark -

void TreadmillingFiber::setEndStateM(state_t s)
{
    if ( s == STATE_WHITE || s == STATE_GREEN || s == STATE_RED )
        mStateM = s;
    else
        throw InvalidParameter("invalid AssemblyState for TreadmillingFiber MINUS_END");
}


void TreadmillingFiber::setEndStateP(state_t s)
{
    if ( s == STATE_WHITE || s == STATE_GREEN || s == STATE_RED )
        mStateP = s;
    else
        throw InvalidParameter("invalid AssemblyState for TreadmillingFiber PLUS_END");
}

//------------------------------------------------------------------------------

void TreadmillingFiber::step()
{
    constexpr size_t P = 0, M = 1;

    if ( mStateP == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real forceP = projectedForceEndP();
        
        // growth is reduced if free monomers are scarce:
        mGrowthP = prop->growing_speed_dt[P] * prop->free_polymer;
        
        // antagonistic force (< 0) decreases assembly rate exponentially
        if (( forceP < 0 ) & ( mGrowthP > 0 ))
            mGrowthP *= std::exp(forceP*prop->growing_force_inv[P]);
    }
    else if ( mStateP == STATE_RED )
    {
        mGrowthP = prop->shrinking_speed_dt[P];
    }
    else
    {
        mGrowthP = 0;
    }
    
    // MINUS_END dynamics
    if ( mStateM == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real forceM = projectedForceEndM();
        
        // growth is reduced if free monomers are scarce:
        mGrowthM = prop->growing_speed_dt[M] * prop->free_polymer;
        
        // antagonistic force (< 0) decreases assembly rate exponentially
        if (( forceM < 0 ) & ( mGrowthM > 0 ))
            mGrowthM *= std::exp(forceM*prop->growing_force_inv[M]);
    }
    else if ( mStateM == STATE_RED )
    {
        mGrowthM = prop->shrinking_speed_dt[M];
    }
    else
    {
        mGrowthM = 0;
    }
    
    real len = length();
    real inc = mGrowthP + mGrowthM;
    if ( len + inc < prop->min_length )
    {
        if ( !prop->persistent )
        {
            delete(this);
            return;
        }
    }
    else if ( len + inc < prop->max_length )
    {
        if ( mGrowthM != 0 ) growM(mGrowthM);
        if ( mGrowthP != 0 ) growP(mGrowthP);
    }
    else if ( len < prop->max_length )
    {
        // the remaining possible growth is distributed to the two ends:
        inc = ( prop->max_length - len ) / inc;
        mGrowthM *= inc;
        mGrowthP *= inc;
        if ( mGrowthM != 0 ) growM(mGrowthM);
        if ( mGrowthP != 0 ) growP(mGrowthP);
    }
    else // len > prop->max_length
    {
        mGrowthM = 0;
        mGrowthP = 0;
    }

    Fiber::step();
    //std::clog << reference() << " P " << mGrowthP << " M " << mGrowthM << " len " << length() << "\n";
}


//------------------------------------------------------------------------------
#pragma mark -


void TreadmillingFiber::write(Outputter& out) const
{
    Fiber::write(out);

    // write variables describing the dynamic state of the ends:
    writeHeader(out, TAG_DYNAMIC);
    out.writeUInt16(mStateM);
    out.writeUInt16(0);
    out.writeUInt16(mStateP);
    out.writeUInt16(0);
}


void TreadmillingFiber::readEndState(Inputter& in)
{
#ifdef BACKWARD_COMPATIBILITY
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
    if ( tag == TAG_DYNAMIC )
        readEndState(in);
    else
    {
#ifdef BACKWARD_COMPATIBILITY
        if ( in.formatID() < 44 )
            readEndState(in);
#endif
        Fiber::read(in, sim, tag);
    }
}

