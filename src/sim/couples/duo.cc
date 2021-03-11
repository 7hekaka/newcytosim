// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "duo.h"
#include "duo_prop.h"
#include "object_set.h"
#include "random.h"
#include "modulo.h"
#include "space.h"
#include "sim.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------

Duo::Duo(DuoProp const* p, Vector const& w)
: Couple(p, w), mActive(0), prop(p)
{
}

Duo::~Duo()
{
    prop = nullptr;
}

//------------------------------------------------------------------------------

void Duo::activate()
{
    mActive = 1;
    gspTime = RNG.exponential();
}

void Duo::deactivate()
{
    mActive = 0;
}

//------------------------------------------------------------------------------

void Duo::stepFF()
{
    diffuse();
    
    // confinement:
    if ( prop->confine == CONFINE_INSIDE )
    {
        prop->confine_space_ptr->bounce(cPos);
    }
    else if ( prop->confine == CONFINE_ON )
    {
        cPos = prop->confine_space_ptr->project(cPos);
    }    
    
    // check activity
    ///@todo better Duo::activation criteria
    if ( prop->activation_space_ptr->inside(cPos) )
        activate();

    
    // activity
    if ( mActive )
    {
        // spontaneous de-activation:
        gspTime -= prop->deactivation_rate_dt;
        if ( gspTime <= 0 )
        {
            deactivate();
            // test fraction of time when it is inactive:
            if ( RNG.test(-gspTime/prop->deactivation_rate_dt) )
                return;
        }
    
        // hands may bind:
        cHand1->stepUnattached(simul(), cPos);
        if ( !prop->trans_activated )
            cHand2->stepUnattached(simul(), cPos);
    }
}


/**
 test for spontaneous de-activation
 */
void Duo::deactivation()
{
    gspTime -= prop->deactivation_rate_dt;
    if ( gspTime <= 0 )
        deactivate();
}


/**
 Simulates:
 - attachment of cHand2
 - attached activity of cHand1
 .
 */
void Duo::stepAF()
{
    if ( mActive && prop->vulnerable )
        deactivation();
    
    //we use cHand1->pos() first, because stepUnloaded() may detach cHand1
    cHand2->stepUnattached(simul(), cHand1->outerPos());

    if ( cHand1->testDetachment() )
        cHand1->stepUnloaded();
    else
        cHand1->detach();
}


/**
 Simulates:
 - attachment of cHand1
 - attached activity of cHand2
 .
 */
void Duo::stepFA()
{
    if ( mActive && prop->vulnerable )
        deactivation();
    
    //we use cHand2->pos() first, because stepUnloaded() may detach cHand2
    cHand1->stepUnattached(simul(), cHand2->outerPos());

    if ( cHand2->testDetachment() )
        cHand2->stepUnloaded();
    else
        cHand2->detach();
}


/**
 Simulates:
 - attached activity of cHand1
 - attached activity of cHand2
 .
 */
void Duo::stepAA()
{
    if ( mActive && prop->vulnerable )
        deactivation();

    Vector f = force();
    real fn = f.norm();
    
    if ( cHand1->testKramersDetachment(fn) )
        cHand1->stepLoaded( f);
    else
        cHand1->detach();
    
    if ( cHand2->testKramersDetachment(fn) )
        cHand2->stepLoaded(-f);
    else
        cHand2->detach();
}


//------------------------------------------------------------------------------

void Duo::write(Outputter& out) const
{
    out.writeUInt8(mActive);
    Couple::write(out);
}


void Duo::read(Inputter& in, Simul& sim, ObjectTag tag)
{
#ifdef BACKWARD_COMPATIBILITY
    if ( in.formatID() > 36 )
#endif
    mActive = in.readUInt8();
    Couple::read(in, sim, tag);
}


