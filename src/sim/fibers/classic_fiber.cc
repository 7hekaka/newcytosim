// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "assert_macro.h"
#include "classic_fiber.h"
#include "classic_fiber_prop.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "messages.h"
#include "space.h"


//------------------------------------------------------------------------------

ClassicFiber::ClassicFiber(ClassicFiberProp const* p) : Fiber(p)
{
    mStateM = STATE_WHITE;
    mStateP = STATE_WHITE;
}


ClassicFiber::~ClassicFiber()
{
}


void ClassicFiber::setEndStateM(state_t s)
{
    if (( s==STATE_WHITE ) | ( s==STATE_GREEN ) | ( s==STATE_RED ))
        mStateM = s;
    else
        throw InvalidParameter("invalid AssemblyState ("+std::to_string(s)+") for ClassicFiber minus end");
}


void ClassicFiber::setEndStateP(state_t s)
{
    if (( s==STATE_WHITE ) | ( s==STATE_GREEN ) | ( s==STATE_RED ))
        mStateP = s;
    else
        throw InvalidParameter("invalid AssemblyState ("+std::to_string(s)+") for ClassicFiber plus end");
}


//------------------------------------------------------------------------------
#pragma mark -

/** 
 The catastrophe rate depends on the growth rate of the corresponding tip,
 which is itself reduced by antagonistic force. 
 The correspondance is : 1/rate = a + b * growthSpeed.
 
 For no force on the growing tip: rate = catastrophe_rate * timestep
 For very large forces          : rate = catastrophe_rate_stalled * timestep
 
 cf. `Dynamic instability of MTs is regulated by force`
 M.Janson, M. de Dood, M. Dogterom. JCB 2003, Figure 2 C.
 */
real ClassicFiber::stepMinusEnd()
{
    constexpr unsigned M = 1;
    real add = 0;

    if ( mStateM == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real force = projectedForceEndM();
        
        // growth is reduced if free monomers are scarce:
        add = prop()->growing_speed_dt[M] * prop()->free_polymer;
        
        // antagonistic force (< 0) decreases assembly rate exponentially
        if (( force < 0 ) & ( add > 0 ))
            add *= std::exp(force*prop()->growing_force_inv[M]);

        add += prop()->growing_off_speed_dt[M];
        
        // catastrophe may be constant, or it may depend on the growth rate
        real cata;
        if ( prop()->catastrophe_coef[M] > 0 )
            cata = prop()->catastrophe_rate_stalled_dt[M] / ( 1.0 + prop()->catastrophe_coef[M] * add );
        else
            cata = prop()->catastrophe_rate_dt[M];

        if ( RNG.test(cata) )
            mStateM = STATE_RED;
    }
    else if ( mStateM == STATE_RED )
    {
        add = prop()->shrinking_speed_dt[M];
        
        real rate = prop()->rescue_prob[M];
        if ( RNG.test(rate) )
            mStateM = STATE_GREEN;
    }
    
    return add;
}


real ClassicFiber::stepPlusEnd()
{
    constexpr unsigned P = 0;
    real add = 0;

    if ( mStateP == STATE_GREEN )
    {
        // calculate the force acting on the point at the end:
        real force = projectedForceEndP();
        
        // growth is reduced if free monomers are scarce:
        add = prop()->growing_speed_dt[P] * prop()->free_polymer;
        
        // antagonistic force (< 0) decreases assembly rate exponentially
        if (( force < 0 ) & ( add > 0 ))
            add *= std::exp(force*prop()->growing_force_inv[P]);

        add += prop()->growing_off_speed_dt[P];
        
        // catastrophe may be constant, or it may depend on the growth rate
        real cata;
        if ( prop()->catastrophe_coef[P] > 0 )
            cata = prop()->catastrophe_rate_stalled_dt[P] / ( 1.0 + prop()->catastrophe_coef[P] * add );
        else
            cata = prop()->catastrophe_rate_dt[P];
        
        //printf("ClassicFiber %5u: force %9.2f growth %9.6f cata %9.6f\n", identity(), force, addP, cata);
        
#if NEW_CATASTROPHE_OUTSIDE
        // Catastrophe rate is multiplied if the plus end is outside
        if ( prop()->catastrophe_space && prop()->catastrophe_space->outside(posEndP()) )
        {
            LOG_ONCE("*** Fiber's plus-end catastrophe changed outside Space\n");
            cata *= prop()->catastrophe_outside;
        }
#endif

#if NEW_LENGTH_DEPENDENT_CATASTROPHE
        /*
         Ad-hoc length dependence, used to simulate S. pombe with catastrophe_length=5
         Foethke et al. Molecular Systems Biology 5:241 - 2009
         */
        if ( prop()->catastrophe_length > 0 )
        {
            LOG_ONCE("Using ad-hoc length-dependent catastrophe rate\n");
            cata *= length() / prop()->catastrophe_length;
        }
#endif
        if ( RNG.test(cata) )
            mStateP = STATE_RED;
    }
    else if ( mStateP == STATE_RED )
    {
        real force = -projectedForceEndP();
        
        add = prop()->shrinking_speed_dt[P] * exp(force*prop()->shrinking_force_inv[P]);

        //printf("ClassicFiber %5u: force %9.2f shrink %9.6f\n", identity(), force, addP);
        
        if ( RNG.test(prop()->rescue_prob[P]) )
            mStateP = STATE_GREEN;
    }
    
    return add;
}


void ClassicFiber::step()
{
    real addM = 0;
    constexpr unsigned M = 1;
    // STATE_WHITE is a dormant state from which you can exit by 'rebirth'
    if ( mStateM == STATE_WHITE )
    {
        if ( RNG.test(prop()->rebirth_prob[M]) )
            setEndStateM(STATE_GREEN);
    }
    else
        addM = stepMinusEnd();
    
    real addP = 0;
    constexpr unsigned P = 0;
    // STATE_WHITE is a dormant state from which you can exit by 'rebirth'
    if ( mStateP == STATE_WHITE )
    {
        if ( RNG.test(prop()->rebirth_prob[P]) )
            setEndStateP(STATE_GREEN);
    }
    else
        addP = stepPlusEnd();

    if ( Fiber::updateLength(addM, addP) )
        Fiber::step();
}


//------------------------------------------------------------------------------
#pragma mark -


void ClassicFiber::write(Outputter& out) const
{
    Fiber::write(out);
    
    // write variables describing the dynamic state of the ends, using 8 bytes:
    writeMarker(out, DYNAMIC_TAG);
    out.writeUInt16(mStateM);
    out.writeUInt16(0);
    out.writeUInt16(mStateP);
    out.writeUInt16(0);
}


void ClassicFiber::readEndStates(Inputter& in)
{
#if BACKWARD_COMPATIBILITY < 54
    if ( in.formatID() < 54 )
    {
        unsigned m = 0, p = 0;
        if ( in.formatID() < 42 )
            p = in.readUInt8();
        else
        {
            m = in.readUInt16();
            p = in.readUInt16();
        }
        if ( in.formatID() < 46 )
            setEndStateP(p);
        else
        {
            setEndStateM(m);
            setEndStateP(p);
        }
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



void ClassicFiber::read(Inputter& in, Simul& sim, ObjectTag tag)
{
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

