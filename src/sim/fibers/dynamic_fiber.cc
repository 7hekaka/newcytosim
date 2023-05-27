// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2021

#include "dim.h"
#include "smath.h"
#include "assert_macro.h"
#include "dynamic_fiber.h"
#include "dynamic_fiber_prop.h"
#include "exceptions.h"
#include "iowrapper.h"
//#include "object_set.h"
//#include "simul_part.h"
#include "space.h"


/**
 By default, the plus end is growing and the minus end is shrinking
 */
DynamicFiber::DynamicFiber(DynamicFiberProp const* p) : Fiber(p)
{
    initM();
    initP();
    stabilized_[0] = 0;
    stabilized_[1] = 0;
}


DynamicFiber::~DynamicFiber()
{
    //std::clog << "chewed " << abscissaM() << "\n";
}


real DynamicFiber::chewingUnits(int end)
{
#if NEW_FIBER_END_CHEW
    const real sup = prop()->max_chewing_speed_dt;
    real res = std::tanh(fChew[end]/sup) * (sup/prop()->unit_length);
    fChew[end] = 0;
    return res;
#else
    return 0;
#endif
}

//------------------------------------------------------------------------------
#pragma mark - MINUS END

/** set minus end as STATE_WHITE (inactive) by default */
void DynamicFiber::initM()
{
    unitM[0] = 0;
    unitM[1] = 2;
    mStateM  = calculateStateM();

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
    return mStateM;
}


void DynamicFiber::setEndStateM(state_t s)
{
    if ( s < 0 || 4 < s )
        throw InvalidParameter("Invalid AssemblyState for DynamicFiber minus end");
    
    if ( s != mStateM )
    {
        mStateM = s;
        unitM[1] = ( 4 - s ) / 2;
        unitM[0] = ( 4 - s ) & 1;
        if ( s != STATE_WHITE )
        {
            assert_true( 0==unitM[0] || 1==unitM[0] );
            assert_true( 0==unitM[1] || 1==unitM[1] );
        }
        assert_true( mStateM == calculateStateM() );
    }
}


void DynamicFiber::addUnitM()
{
    unitM[1] = unitM[0];
    unitM[0] = 1;
    nextGrowthM += RNG.exponential();
}


// remove last unit, while penultimate position can be a GTP survivor
void DynamicFiber::removeUnitM()
{
    unitM[0] = unitM[1];
    unitM[1] = RNG.test(prop()->unhydrolyzed_prob[1]);
    nextShrinkM += RNG.exponential();
    mStateM = calculateStateM();
    //std::cout << reference() << "mStateM = " << mStateM << '\n';
}


/**
 Minus end can only shrink, and so far the state vector unitM[] is ignored
 */
int DynamicFiber::stepMinusEnd()
{
    constexpr size_t M = 1;
    int res = 0;
    
#if ( NEW_FIBER_END_CHEW & 2 )
    real chewed = 0;
    if ( fChew[M] > 0 )
    {
        chewed = chewingUnits(M);
        //std::clog << reference() << " chew M " << chewed/time_step(simul()) << " unit/s\n";
    }
#else
    constexpr real chewed = 0;
#endif
    
    if ( mStateM == STATE_RED )
    {
        nextShrinkM -= prop()->shrinking_rate_dt[M] + chewed;
        while ( nextShrinkM <= 0 )
        {
            removeUnitM();
            --res;
        }
    }
    else if ( mStateM > STATE_WHITE )
    {
        // calculate the force acting on the point at the end:
        real forceM = projectedForceEndM();

        // growth is reduced if free monomers are scarce:
        real growth = prop()->growing_rate_dt[M] * prop()->free_polymer;

        // antagonistic force (< 0) decreases assembly rate exponentially
        if (( forceM < 0 ) & ( growth > 0 ))
            growth *= std::exp(forceM*prop()->growing_force_inv[M]);

        real hydrol = prop()->hydrolysis_rate_2dt[M];

#if OLD_DYNAMIC_ZONE
        // change Hydrolysis rate if plus end is far from origin:
        if ( posEndM().normSqr() > prop()->zone_radius_sqr )
            hydrol = prop()->zone_hydrolysis_rate_2dt[M];

        if ( prop()->zone_space_ptr && !prop()->zone_space_ptr->inside(posEndM()) )
            hydrol = prop()->zone_hydrolysis_rate_2dt[M];
#endif

        // @todo detach_rate should depend on the state of the subunit
        real shrink = prop()->growing_off_rate_dt[M] + chewed;

        nextGrowthM -= growth;
        nextShrinkM -= shrink;
        nextHydrolM -= hydrol;
        
        while (( nextGrowthM < 0 ) | ( nextShrinkM < 0 ) | ( nextHydrolM < 0 ))
        {
            // Select the earliest event (in most cases, only one event will fire up)
            int ii = ( nextHydrolM * growth < nextGrowthM * hydrol );
            if ( nextShrinkM < 0 )
            {
                if (( ii == 0 ) & ( nextShrinkM * growth < nextGrowthM * shrink ))
                    ii = 2;
                else if ( nextShrinkM * hydrol < nextHydrolM * shrink )
                    ii = 2;
            }
            //printf("%s %6.2f %6.2f %6.2f: %i\n", reference().c_str(), nextGrowthM, nextHydrolM, nextShrinkM, ii);
            switch ( ii )
            {
                case 0:
                    assert_true(nextGrowthM < 0);
                    // add fresh unit, shifting old terminal to penultimate position
                    addUnitM();
                    ++res;
                    break;

                case 1:
                    assert_true(nextHydrolM < 0);
                    // hydrolyze one of the unit with equal chance:
                    unitM[RNG.flip()] = 0;
                    nextHydrolM += RNG.exponential();
                    break;

                case 2:
                    assert_true(nextShrinkM < 0);
                    // remove last unit, while penultimate position can be a GTP survivor
                    removeUnitM();
                    --res;
                    break;

            }
            mStateM = calculateStateM();
        }
    }
    return res;
}


//------------------------------------------------------------------------------
#pragma mark - PLUS END

/** set as STATE_GREEN (growing) by default */
void DynamicFiber::initP()
{
    unitP[0] = 1;
    unitP[1] = 1;
    mStateP = calculateStateP();
    
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
    return mStateP;
}


void DynamicFiber::setEndStateP(state_t s)
{
    if ( s < 0 || 4 < s )
        throw InvalidParameter("invalid AssemblyState ("+std::to_string(s)+") for DynamicFiber plus end");
    
    if ( s != mStateP )
    {
        mStateP = s;
        unitP[1] = ( 4 - s ) / 2;
        unitP[0] = ( 4 - s ) & 1;
        if ( s != STATE_WHITE )
        {
            assert_true( 0==unitP[0] || 1==unitP[0] );
            assert_true( 0==unitP[1] || 1==unitP[1] );
        }
        assert_true( mStateP == calculateStateP() );
    }
}


// add unit
void DynamicFiber::addUnitP()
{
    unitP[1] = unitP[0];
    unitP[0] = 1;
    mStateP = calculateStateP();
    nextGrowthP += RNG.exponential();
}


// remove last unit, while penultimate position can be a GTP survivor
void DynamicFiber::removeUnitP()
{
    unitP[0] = unitP[1];
    unitP[1] = RNG.test(prop()->unhydrolyzed_prob[0]);
    nextShrinkP += RNG.exponential();
    mStateP = calculateStateP();
    //std::cout << reference() << "mStateP = " << mStateP << '\n';
}


/**
 Using a modified Gillespie scheme with a variable rate.
 
 returns the number of units added (if result > 0) or removed (if < 0)
 */
int DynamicFiber::stepPlusEnd()
{
    constexpr size_t P = 0;
    int res = 0;
    
#if ( NEW_FIBER_END_CHEW & 1 )
    real chewed = 0;
    if ( fChew[P] > 0 )
    {
        chewed = chewingUnits(P);
        //std::clog << reference() << " chew P " << chewed/time_step(simul()) << " unit/s\n";
    }
#else
    constexpr real chewed = 0;
#endif
    
    if ( mStateP == STATE_RED )
    {
        nextShrinkP -= prop()->shrinking_rate_dt[P] + chewed;
        while ( nextShrinkP <= 0 )
        {
            removeUnitP();
            --res;
        }
    }
    else
    {
        // calculate the force acting on the point at the end:
        real forceP = projectedForceEndP();
        
        // growth is reduced if free monomers are scarce:
        real growth = prop()->growing_rate_dt[P] * prop()->free_polymer;
        
        // antagonistic force (< 0) decreases assembly rate exponentially
        if (( forceP < 0 ) & ( growth > 0 ))
            growth *= std::exp(forceP*prop()->growing_force_inv[P]);

        real hydrol = prop()->hydrolysis_rate_2dt[P];
        
#if OLD_DYNAMIC_ZONE
        // change Hydrolysis rate if plus end is far from origin:
        if ( posEndP().normSqr() > prop()->zone_radius_sqr )
            hydrol = prop()->zone_hydrolysis_rate_2dt[P];
        
        if ( prop()->zone_space_ptr && !prop()->zone_space_ptr->inside(posEndP()) )
            hydrol = prop()->zone_hydrolysis_rate_2dt[P];
#endif
        
        // @todo detach_rate should depend on the state of the subunit
        real shrink = prop()->growing_off_rate_dt[P] + chewed;

        nextGrowthP -= growth;
        nextShrinkP -= shrink;
        nextHydrolP -= hydrol;
        
        while (( nextGrowthP < 0 ) | ( nextShrinkP < 0 ) | ( nextHydrolP < 0 ))
        {
            // Select the earliest event (in most cases, only one event will fire up)
            int ii = ( nextHydrolP * growth < nextGrowthP * hydrol );
            if ( nextShrinkP < 0 )
            {
                if (( ii == 0 ) & ( nextShrinkP * growth < nextGrowthP * shrink ))
                    ii = 2;
                else if ( nextShrinkP * hydrol < nextHydrolP * shrink )
                    ii = 2;
            }
            //printf("%s %6.2f %6.2f %6.2f: %i\n", reference().c_str(), nextGrowthP, nextHydrolP, nextShrinkP, ii);
            switch ( ii )
            {
                case 0:
                    assert_true(nextGrowthP < 0);
                    // add fresh unit, shifting old terminal to penultimate position
                    addUnitP();
                    ++res;
                    break;
                    
                case 1:
                    assert_true(nextHydrolP < 0);
                    // hydrolyze one of the unit with equal chance:
                    unitP[RNG.flip()] = 0;
                    nextHydrolP += RNG.exponential();
                    break;

                case 2:
                    assert_true(nextShrinkP < 0);
                    // remove last unit, while penultimate position can be a GTP survivor
                    removeUnitP();
                    --res;
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
    //std::clog << prop()->name() << " " << this << " " << stabilized_[0] << "  " << stabilized_[1] << "\n";
    real addM = 0;
    constexpr size_t M = 1;
    if ( stabilized_[M] > 0 )
    {
        stabilized_[M] = 0;
    }
    else if ( mStateM == STATE_WHITE )
    {
#if ( NEW_FIBER_END_CHEW & 2 )
        if ( fChew[M] > 0 )
        {
            nextShrinkM -= chewingUnits(M);
            while ( nextShrinkM <= 0 )
            {
                removeUnitM();
                --addM;
            }
        }
#endif
        // STATE_WHITE is a dormant state from which you can exit by 'rebirth'
        if ( RNG.test(prop()->rebirth_prob[M]) )
            setEndStateM(STATE_GREEN);
    }
    else
        addM = stepMinusEnd();
    
    real addP = 0;
    constexpr size_t P = 0;
    if ( stabilized_[P] > 0 )
    {
        // Special case for Kinetochores: grow slowly with no catastrophe...
        real growth = prop()->growing_rate_dt[P] * prop()->free_polymer;
        nextGrowthP -= growth / stabilized_[P];
        while ( nextGrowthP < 0 )
        {
            addUnitP();
            ++addP;
        }
        stabilized_[P] = 0;
    }
    else if ( mStateP == STATE_WHITE )
    {
#if ( NEW_FIBER_END_CHEW & 1 )
        if ( fChew[P] > 0 )
        {
            nextShrinkP -= chewingUnits(P);
            while ( nextShrinkP <= 0 )
            {
                removeUnitP();
                --addP;
            }
        }
#endif
        // STATE_WHITE is a dormant state from which you can exit by 'rebirth'
        if ( RNG.test(prop()->rebirth_prob[P]) )
            setEndStateP(STATE_GREEN);
    }
    else
        addP = stepPlusEnd();

    addM *= prop()->unit_length;
    addP *= prop()->unit_length;
    
    if ( Fiber::updateLength(addM, addP, false) )
        Fiber::step();
}


//------------------------------------------------------------------------------
#pragma mark -


void DynamicFiber::write(Outputter& out) const
{
    Fiber::write(out);
    
    long addM = std::lround(cDeltaM/prop()->unit_length);
    long addP = std::lround(cDeltaP/prop()->unit_length);
    // write variables describing the dynamic state of the ends:
    writeMarker(out, TAG_DYNAMIC);
    out.writeUInt8(unitM[0]);
    out.writeUInt8(unitM[1]);
    out.writeInt16(addM);
    out.writeUInt8(unitP[0]);
    out.writeUInt8(unitP[1]);
    out.writeInt16(addP);
}


void DynamicFiber::readEndStates(Inputter& in)
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
        in.readUInt16(); //cDeltaM = in.readInt16() * prop()->unit_length;
        unitP[0] = in.readUInt8();
        unitP[1] = in.readUInt8();
        in.readUInt16(); //cDeltaP = in.readInt16() * prop()->unit_length;
    }
    mStateM = calculateStateM();
    mStateP = calculateStateP();
}


void DynamicFiber::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    if ( tag == TAG_DYNAMIC )
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
