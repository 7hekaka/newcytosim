// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "smath.h"
#include "assert_macro.h"
#include "dynamic_fiber.h"
#include "dynamic_fiber_prop.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul.h"
#include "space.h"


/**
 By default, both ends are growing
 */
DynamicFiber::DynamicFiber(DynamicFiberProp const* p) : Fiber(p), prop(p)
{
    // set PLUS_END as growing
    unitP[0] = 1;
    unitP[1] = 1;
    unitP[2] = 1;
    mStateP  = calculateStateP();
    mGrowthP = 0;
    
    nextGrowthP = RNG.exponential();
    nextHydrolP = RNG.exponential();
    nextShrinkP = RNG.exponential();

    // set MINUS_END as growing
    unitM[0] = 1;
    unitM[1] = 1;
    unitM[2] = 1;
    mStateM  = calculateStateM();
    mGrowthM = 0;

    nextGrowthM = RNG.exponential();
    nextHydrolM = RNG.exponential();
    nextShrinkM = RNG.exponential();
}


DynamicFiber::~DynamicFiber()
{
    prop = nullptr;
}


//------------------------------------------------------------------------------
#pragma mark -


state_t DynamicFiber::calculateStateM() const
{ 
    return 4 - unitM[0] - 2 * unitM[1];
}


state_t DynamicFiber::endStateM() const
{
    return STATE_WHITE;
    assert_true( mStateM == calculateStateM() );
    return mStateM;
}


void DynamicFiber::setEndStateM(state_t s)
{
    if ( s < 1 || 4 < s )
        throw InvalidParameter("invalid AssemblyState for DynamicFiber MINUS_END");
    
    if ( s != mStateM )
    {
        mStateM = s;
        unitM[1] = ( 4 - s ) / 2;
        unitM[0] = ( 4 - s - 2*unitM[1] );
        assert_true( 0==unitP[0] || unitP[0]==1 );
        assert_true( 0==unitP[1] || unitP[1]==1 );
        assert_true( mStateM == calculateStateM() );
    }
}


/**
 Minus end can shrink
 */
int DynamicFiber::stepMinusEnd()
{
    int res = 0;
    real chewing = 0;
    
    // add chewing rate to stochastic off rate:
    
#if NEW_FIBER_CHEW
    
    // convert chewing rate to stochastic off rate:
    if ( frChewM > prop->max_chewing_speed_dt )
        chewing = prop->max_chewing_speed_dt / prop->unit_length;
    else
        chewing = frChewM / prop->unit_length;
    
    //std::clog << " chewing rate M " << chewing_rate / prop->time_step << '\n';
    frChewM = 0;
    
#endif
    
    nextShrinkM -= prop->shrinking_rate_dt[1] + chewing;
    while ( nextShrinkM <= 0 )
    {
        --res;
        nextShrinkM += RNG.exponential();
    }
    return res;
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 The microscopic state correspond to:
 - STATE_GREEN for growth,
 - STATE_RED for shrinkage
 .
 */
state_t DynamicFiber::calculateStateP() const
{
    return 4 - unitP[0] - 2 * unitP[1];
}


state_t DynamicFiber::endStateP() const
{
    assert_true( mStateP == calculateStateP() );
    return mStateP;
}


void DynamicFiber::setEndStateP(state_t s)
{
    if ( s < 1 || 4 < s )
        throw InvalidParameter("invalid AssemblyState for DynamicFiber PLUS_END");
    
    if ( s != mStateP )
    {
        mStateP = s;
        unitP[1] = ( 4 - s ) / 2;
        unitP[0] = ( 4 - s - 2*unitP[1] );
        assert_true( 0==unitP[0] || unitP[0]==1 );
        assert_true( 0==unitP[1] || unitP[1]==1 );
        assert_true( mStateP == calculateStateP() );
    }
}


/**
 Using a modified Gillespie scheme with a variable rate.
 
 returns the number of units added (if result > 0) or removed (if < 0)
 */
int DynamicFiber::stepPlusEnd()
{
    constexpr size_t P = 0;
    int res = 0;
    real chewing = 0;
    
#if NEW_FIBER_CHEW
    
    // convert chewing rate to stochastic off rate:
    ///@todo implement smooth saturation using logistic function
    if ( frChewP > prop->max_chewing_speed_dt )
        chewing = prop->max_chewing_speed_dt / prop->unit_length;
    else
        chewing = frChewP / prop->unit_length;
    
    //std::clog << " chewing rate P " << chewing / prop->time_step << '\n';
    frChewP = 0;
#endif
    
    if ( mStateP == STATE_RED )
    {
        nextShrinkP -= prop->shrinking_rate_dt[P] + chewing;
        while ( nextShrinkP <= 0 )
        {
        	// remove last unit, with a finite probability that a GTP-tubulin is encountered along the lattice
			unitP[0] = unitP[1];
            unitP[1] = unitP[2];
			unitP[2] = RNG.test(prop->unhydrolyzed_prob[P]);
			--res;
            nextShrinkP += RNG.exponential();
            mStateP = calculateStateP();
            //std::cout << "mStateP = " << mStateP << '\n';
        }
    }
    else
    {
        // calculate the force acting on the point at the end:
        real forceP = projectedForceEndP();
        
        // growth is reduced if free monomers are scarce:
        real growth = prop->growing_rate_dt[P] * prop->free_polymer;
        
        // antagonistic force (< 0) decreases assembly rate exponentially
        if (( forceP < 0 ) & ( growth > 0 ))
            growth *= std::exp(forceP*prop->growing_force_inv[P]);

        real hydrol = prop->hydrolysis_rate_2dt[P];
        
#if OLD_DYNAMIC_ZONE
        // change Hydrolysis rate if PLUS_END is far from origin:
        if ( posEndP().normSqr() > prop->zone_radius_sqr )
            hydrol = prop->zone_hydrolysis_rate_2dt[P];
        
        if ( prop->zone_space_ptr && !prop->zone_space_ptr->inside(posEndP()) )
            hydrol = prop->zone_hydrolysis_rate_2dt[P];
#endif
        
        // @todo detach_rate should depend on the state of the subunit
        real detach = prop->growing_off_rate_dt[P] + chewing;
        
        nextGrowthP -= growth;
        nextShrinkP -= detach;
        nextHydrolP -= hydrol;
        
        while (( nextGrowthP < 0 ) | ( nextShrinkP < 0 ) | ( nextHydrolP < 0 ))
        {
            // Select the earliest event:
            int ii = sMath::arg_min(nextGrowthP/growth, nextHydrolP/hydrol, nextShrinkP/detach);
            
            switch ( ii )
            {
                case 0:
                    // add fresh unit, shifting old terminal to penultimate position
                    unitP[2] = unitP[1] * RNG.test(prop->unhydrolyzed_prob[P]);
                    unitP[1] = unitP[0];
                    unitP[0] = 1;
                    ++res;
                    nextGrowthP += RNG.exponential();
                    break;
                    
                case 1:
                    // hydrolyze one of the unit with equal chance:
                    unitP[RNG.flip()] = 0;
                    nextHydrolP += RNG.exponential();
                    break;

                case 2:
                    // remove last unit, with a finite probability that a GTP-tubulin is encountered along the lattice
                    unitP[0] = unitP[1];
                    unitP[1] = unitP[2];
                    unitP[2] = RNG.test(prop->unhydrolyzed_prob[P]);
                    --res;
                    nextShrinkP += RNG.exponential();
                    break;

            }
            
            mStateP = calculateStateP();
        }
    }
    return res;
}


//------------------------------------------------------------------------------
#pragma mark -

void DynamicFiber::step()
{
    constexpr size_t P = 0;
    // perform stochastic simulation:
    int incP = stepPlusEnd();
    int incM = stepMinusEnd();

    mGrowthP = incP * prop->unit_length;
    mGrowthM = incM * prop->unit_length;

    if ( incM || incP )
    {
        if ( length() + mGrowthM + mGrowthP < prop->min_length )
        {
            // do something if the fiber is too short:
            if ( !prop->persistent )
            {
                delete(this);
                // exit to avoid doing anything with a dead object:
                return;
            }
            // possibly rescue:
            if ( RNG.test(prop->rebirth_prob[P]) )
                setEndStateP(STATE_GREEN);
        }
        else if ( length() + mGrowthM + mGrowthP < prop->max_length )
        {
            if ( incP ) growP(mGrowthP);
            if ( incM ) growM(mGrowthM);
            //std::clog << reference() << " " << mGrowthM << " " << mGrowthP << " " << length() << '\n';
        }
    }

    Fiber::step();
}


//------------------------------------------------------------------------------
#pragma mark -


void DynamicFiber::write(Outputter& out) const
{
    Fiber::write(out);
    
    // write variables describing the dynamic state of the ends:
    writeHeader(out, TAG_DYNAMIC);
    out.writeUInt8(unitM[0]);
    out.writeUInt8(unitM[1]);
    out.writeUInt16(0);
    out.writeUInt8(unitP[0]);
    out.writeUInt8(unitP[1]);
    out.writeUInt16(0);
}


void DynamicFiber::readEndState(Inputter& in)
{
#ifdef BACKWARD_COMPATIBILITY
    if ( in.formatID() < 54 )
    {
        unitM[0] = in.readUInt8();
        unitM[1] = in.readUInt8();
        unitP[0] = in.readUInt8();
        unitP[1] = in.readUInt8();
    }
    else
#endif
    {
        unitM[0] = in.readUInt8();
        unitM[1] = in.readUInt8();
        in.readUInt16();
        unitP[0] = in.readUInt8();
        unitP[1] = in.readUInt8();
        in.readUInt16();
    }
    mStateM = calculateStateM();
    mStateP = calculateStateP();
}


void DynamicFiber::read(Inputter& in, Simul& sim, ObjectTag tag)
{
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
