// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "assert_macro.h"
#include "growing_fiber.h"
#include "growing_fiber_prop.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul_part.h"
#include "space.h"


//------------------------------------------------------------------------------

GrowingFiber::GrowingFiber(GrowingFiberProp const* p) : Fiber(p), prop(p)
{
    mStateM = STATE_GREEN;
    mStateP = STATE_GREEN;
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
    real addP = 0, addM = 0;
    
    // PLUS_END
    if ( prop->shrink_outside[P] && prop->confine_space_ptr->outside(posEndP()) )
    {
        addP = prop->shrinking_speed_dt[P];
    }
    else if ( mStateP == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real forceP = projectedForceEndP();
        
        // growth is reduced if free monomers are scarce:
        addP = prop->growing_speed_dt[P] * prop->free_polymer;
        
        // antagonistic force (< 0) decreases assembly rate exponentially
        if (( forceP < 0 ) & ( addP > 0 ))
            addP *= std::exp(forceP*prop->growing_force_inv[P]);
        
        addP += prop->growing_off_speed_dt[P];
    }
    
    // MINUS_END
    if ( prop->shrink_outside[M] && prop->confine_space_ptr->outside(posEndM()) )
    {
        addM = prop->shrinking_speed_dt[M];
    }
    else if ( mStateM == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real forceM = projectedForceEndM();
        
        // growth is reduced if free monomers are scarce:
        addM = prop->growing_speed_dt[M] * prop->free_polymer;
        
        // antagonistic force (< 0) decreases assembly rate exponentially
        if (( forceM < 0 ) & ( addM > 0 ))
            addM *= std::exp(forceM*prop->growing_force_inv[M]);

        addM += prop->growing_off_speed_dt[M];
    }

    Fiber::step(addM, addP);
}


//------------------------------------------------------------------------------
#pragma mark -


void GrowingFiber::write(Outputter& out) const
{
    Fiber::write(out);

    // since states are constant, we write growth rates:
    writeHeader(out, TAG_DYNAMIC);
    out.writeFloat(cDeltaM);
    out.writeFloat(cDeltaP);
}


void GrowingFiber::readEndState(Inputter& in)
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
    }
}


void GrowingFiber::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    if ( tag == TAG_DYNAMIC )
        readEndState(in);
    else
    {
#if BACKWARD_COMPATIBILITY < 44
        if ( tag == TAG && in.formatID() < 44 )
            readEndState(in);
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

