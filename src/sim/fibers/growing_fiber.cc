// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "assert_macro.h"
#include "growing_fiber.h"
#include "growing_fiber_prop.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul.h"
#include "space.h"


//------------------------------------------------------------------------------

GrowingFiber::GrowingFiber(GrowingFiberProp const* p) : Fiber(p), prop(p)
{
    mStateM = STATE_GREEN;
    mStateP = STATE_GREEN;
    mGrowthM = 0;
    mGrowthP = 0;
}


GrowingFiber::~GrowingFiber()
{
    prop = nullptr;
}


//------------------------------------------------------------------------------
#pragma mark -


void GrowingFiber::setEndStateM(state_t s)
{
    if ( s == STATE_WHITE || s == STATE_GREEN )
        mStateM = s;
    else
        throw InvalidParameter("invalid AssemblyState for a GrowingFiber");
}


void GrowingFiber::setEndStateP(state_t s)
{
    if ( s == STATE_WHITE || s == STATE_GREEN )
        mStateP = s;
    else
        throw InvalidParameter("invalid AssemblyState for a GrowingFiber");
}

//------------------------------------------------------------------------------

void GrowingFiber::step()
{
    constexpr size_t P = 0, M = 1;

    // PLUS_END
    if ( prop->shrink_outside[P] && prop->confine_space_ptr->outside(posEndP()) )
    {
        mGrowthP = prop->shrinking_speed_dt[P];
    }
    else if ( mStateP == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real forceP = projectedForceEndP();
        
        // growth is reduced if free monomers are scarce:
        mGrowthP = prop->growing_speed_dt[P] * prop->free_polymer;
        
        // antagonistic force (< 0) decreases assembly rate exponentially
        if ( forceP < 0 )
            mGrowthP *= std::exp(forceP*prop->growing_force_inv[P]);
        
        mGrowthP += prop->growing_off_speed_dt[P];
    }
    else
    {
        mGrowthP = 0;
    }
    
    // MINUS_END
    if ( prop->shrink_outside[M] && prop->confine_space_ptr->outside(posEndM()) )
    {
        mGrowthM = prop->shrinking_speed_dt[M];
    }
    else if ( mStateM == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real forceM = projectedForceEndM();
        
        // growth is reduced if free monomers are scarce:
        mGrowthM = prop->growing_speed_dt[M] * prop->free_polymer;
        
        // antagonistic force (< 0) decreases assembly rate exponentially
        if ( forceM < 0 )
            mGrowthM *= std::exp(forceM*prop->growing_force_inv[M]);

        mGrowthM += prop->growing_off_speed_dt[M];
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
            // the fiber is too short, we delete it:
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
        if ( mGrowthM != 0 ) growM(inc*mGrowthM);
        if ( mGrowthP != 0 ) growP(inc*mGrowthP);
    }
    else
    {
        mGrowthM = 0;
        mGrowthP = 0;
    }

    Fiber::step();
}


//------------------------------------------------------------------------------
#pragma mark -


void GrowingFiber::write(Outputter& out) const
{
    Fiber::write(out);

    /// write variables describing the dynamic state of the ends:
    writeHeader(out, TAG_DYNAMIC);
    out.writeFloat(mGrowthM);
    out.writeFloat(0.0);
    out.writeFloat(mGrowthP);
    out.writeFloat(0.0);
}


void GrowingFiber::readEndState(Inputter& in)
{
#ifdef BACKWARD_COMPATIBILITY
    if ( in.formatID() < 54 )
    {
        mGrowthM = in.readFloat();
        if ( in.formatID() > 45 )
            mGrowthP = in.readFloat();
    }
    else
#endif
    {
        mGrowthM = in.readFloat();
        in.readFloat();
        mGrowthP = in.readFloat();
        in.readFloat();
    }
}


void GrowingFiber::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    if ( tag == TAG_DYNAMIC )
        readEndState(in);
    else
    {
#ifdef BACKWARD_COMPATIBILITY
        if ( in.formatID() < 44 )
            readEndState(in);
        const real len = length();
#endif
        
        Fiber::read(in, sim, tag);
                
#ifdef BACKWARD_COMPATIBILITY
        if ( tag == TAG && in.formatID() < 46 )
        {
            // adjust growing variable
            mGrowthP = length() - len;
            mGrowthM = 0;
        }
#endif
    }
}

