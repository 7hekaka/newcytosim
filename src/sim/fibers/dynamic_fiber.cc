// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2021

#include "dim.h"
#include "smath.h"
#include "assert_macro.h"
#include "dynamic_fiber.h"
#include "dynamic_fiber_prop.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "space.h"


/**
 By default, the plus end is growing and the minus end is shrinking
 */
DynamicFiber::DynamicFiber(DynamicFiberProp const* p) : Fiber(p), prop(p)
{
    initM();
    initP();
}


DynamicFiber::~DynamicFiber()
{
    prop = nullptr;
}

//------------------------------------------------------------------------------
#pragma mark - MINUS END

/** set MINUS_END as shrinking */
void DynamicFiber::initM()
{
    unitM[0] = 0;
    unitM[1] = 0;
    mStateM  = calculateStateM();
    mGrowthM = 0;

    nextGrowthM = RNG.exponential();
    nextHydrolM = RNG.exponential();
    nextShrinkM = RNG.exponential();
}


state_t DynamicFiber::calculateStateM() const
{ 
    return 4 - unitM[0] - 2 * unitM[1];
}


state_t DynamicFiber::endStateM() const
{
    assert_true( mStateM == calculateStateM() );
    return mStateM;
}


void DynamicFiber::setEndStateM(state_t s)
{
    if ( s < 0 || 4 < s )
        throw InvalidParameter("Invalid AssemblyState for DynamicFiber MINUS_END");
    
    if ( s != mStateM )
    {
        mStateM = s;
        if ( STATE_WHITE < s )
        {
            unitM[1] = ( 4 - s ) / 2;
            unitM[0] = ( 4 - s ) & 1;
            assert_true( 0==unitM[0] || 1==unitM[0] );
            assert_true( 0==unitM[1] || 1==unitM[1] );
            assert_true( mStateM == calculateStateM() );
        }
    }
}


/**
 Minus end can only shrink, and so far the state vector unitM[] is ignored
 */
int DynamicFiber::stepMinusEnd()
{
	constexpr size_t M = 1;
	int res = 0;
    real chewing = 0;
    
    // add chewing rate to stochastic off rate:
#if NEW_FIBER_CHEW
    
    // convert chewing rate to stochastic off rate:
    if ( fChewM > prop->max_chewing_speed_dt )
        chewing = prop->max_chewing_speed_dt / prop->unit_length;
    else
        chewing = fChewM / prop->unit_length;
    
    //std::clog << " chewing rate M " << chewing_rate / prop->time_step << '\n';
    fChewM = 0;
    
#endif
    
    if ( mStateM == STATE_RED )
	{
		nextShrinkM -= prop->shrinking_rate_dt[M] + chewing;
		while ( nextShrinkM <= 0 )
		{
			// remove last unit, with a finite probability that a GTP-tubulin is encountered along the lattice
			unitM[0] = unitM[1];
			unitM[1] = RNG.test(prop->unhydrolyzed_prob[M]);
			--res;
			nextShrinkM += RNG.exponential();
			mStateM = calculateStateM();
		}
	}
    else if ( mStateM > STATE_WHITE )
	{
		// calculate the force acting on the point at the end:
		real forceM = projectedForceEndM();

		// growth is reduced if free monomers are scarce:
		real growth = prop->growing_rate_dt[M] * prop->free_polymer;

		// antagonistic force (< 0) decreases assembly rate exponentially
		if (( forceM < 0 ) & ( growth > 0 ))
			growth *= std::exp(forceM*prop->growing_force_inv[M]);

		real hydrol = prop->hydrolysis_rate_2dt[M];

#if OLD_DYNAMIC_ZONE
        // change Hydrolysis rate if PLUS_END is far from origin:
        if ( posEndM().normSqr() > prop->zone_radius_sqr )
            hydrol = prop->zone_hydrolysis_rate_2dt[M];

        if ( prop->zone_space_ptr && !prop->zone_space_ptr->inside(posEndM()) )
            hydrol = prop->zone_hydrolysis_rate_2dt[M];
#endif

		// @todo detach_rate should depend on the state of the subunit
		real detach = prop->growing_off_rate_dt[M] + chewing;

		nextGrowthM -= growth;
		nextShrinkM -= detach;
		nextHydrolM -= hydrol;

		while (( nextGrowthM < 0 ) | ( nextShrinkM < 0 ) | ( nextHydrolM < 0 ))
		{
			// Select the earliest event:
			int ii = sMath::arg_min(nextGrowthM/growth, nextHydrolM/hydrol, nextShrinkM/detach);

			switch ( ii )
			{
				case 0:
					// add fresh unit, shifting old terminal to penultimate position
					unitM[1] = unitM[0];
					unitM[0] = 1;
					++res;
					nextGrowthM += RNG.exponential();
					break;

				case 1:
					// hydrolyze one of the unit with equal chance:
					unitM[RNG.flip()] = 0;
					nextHydrolM += RNG.exponential();
					break;

				case 2:
					// remove last unit, with a finite probability that a GTP-tubulin is encountered along the lattice
					unitM[0] = unitM[1];
					unitM[1] = RNG.test(prop->unhydrolyzed_prob[M]);
					--res;
					nextShrinkM += RNG.exponential();
					break;

			}
			mStateM = calculateStateM();
		}
	}
    return res;
}


//------------------------------------------------------------------------------
#pragma mark - PLUS END

/** set PLUS_END as growing */
void DynamicFiber::initP()
{
    unitP[0] = 1;
    unitP[1] = 1;
    mStateP  = calculateStateP();
    mGrowthP = 0;
    
    nextGrowthP = RNG.exponential();
    nextHydrolP = RNG.exponential();
    nextShrinkP = RNG.exponential();
}


/**
 The microscopic state correspond to:
 - 1 = STATE_GREEN for growth,
 - 4 = STATE_RED for shrinkage
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
    if ( s < 0 || 4 < s )
        throw InvalidParameter("invalid AssemblyState for DynamicFiber PLUS_END");
    
    if ( s != mStateP )
    {
        mStateP = s;
        unitP[1] = ( 4 - s ) / 2;
        unitP[0] = ( 4 - s ) & 1;
        assert_true( 0==unitP[0] || 1==unitP[0] );
        assert_true( 0==unitP[1] || 1==unitP[1] );
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
    if ( fChewP > prop->max_chewing_speed_dt )
        chewing = prop->max_chewing_speed_dt / prop->unit_length;
    else
        chewing = fChewP / prop->unit_length;
    
    //std::clog << " chewing rate P " << chewing / prop->time_step << '\n';
    fChewP = 0;
#endif
    
    if ( mStateP == STATE_RED )
    {
        nextShrinkP -= prop->shrinking_rate_dt[P] + chewing;
        while ( nextShrinkP <= 0 )
        {
        	// remove last unit, with a finite probability that a GTP-tubulin is encountered along the lattice
			unitP[0] = unitP[1];
            unitP[1] = RNG.test(prop->unhydrolyzed_prob[P]);
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
                    unitP[1] = RNG.test(prop->unhydrolyzed_prob[P]);
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
    constexpr size_t M = 1;

    // perform stochastic simulation:
    int incP = stepPlusEnd();
    int incM = mStateM ? stepMinusEnd() : 0;

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
            if ( mStateM == STATE_RED && RNG.test(prop->rebirth_prob[M]) )
            	setEndStateM(STATE_GREEN);

            if ( mStateP == STATE_RED && RNG.test(prop->rebirth_prob[P]) )
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
#if BACKWARD_COMPATIBILITY < 54
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
#if BACKWARD_COMPATIBILITY < 44
        if ( tag == TAG && in.formatID() < 44 )
            readEndState(in);
#endif
        Fiber::read(in, sim, tag);
    }
}
