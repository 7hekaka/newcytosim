// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "duo.h"
#include "duo_prop.h"
#include "object_set.h"
#include "random.h"
#include "modulo.h"
#include "space.h"
#include "cymdef.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------

Duo::Duo(DuoProp const* p, Vector const& w)
: Couple(p, w), active_(0), prop(p)
{
}

Duo::~Duo()
{
    prop = nullptr;
}

//------------------------------------------------------------------------------

void Duo::activate()
{
    active_ = 1;
    countdown_ = RNG.exponential();
}

void Duo::deactivate()
{
    active_ = 0;
}

//------------------------------------------------------------------------------

void Duo::stepFF()
{
    diffuse();
    
    // check activity
    ///@todo better Duo::activation criteria
    if ( prop->activation_space_ptr->inside(cPos) )
        activate();
    
    // activity
    if ( active_ )
    {
        // spontaneous de-activation:
        countdown_ -= prop->deactivation_rate_dt;
        if ( countdown_ <= 0 )
        {
            deactivate();
            // test fraction of time when it is inactive:
            if ( RNG.test(-countdown_/prop->deactivation_rate_dt) )
                return;
        }
    
        // hands may bind:
        if ( RNG.flip() )
            cHand1->stepUnattached(simul(), cPos);
        else if ( !prop->trans_activated )
            cHand2->stepUnattached(simul(), cPos);
    }
}


/**
 test for spontaneous de-activation
 */
void Duo::tryDeactivate()
{
    countdown_ -= prop->deactivation_rate_dt;
    if ( countdown_ <= 0 )
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
    if ( active_ && prop->vulnerable )
        tryDeactivate();
    
    //we use cHand1->pos() first, because stepUnloaded() may detach cHand1
    cHand2->stepUnattached(simul(), cHand1->outerPos());

    if ( cHand1->checkDetachment() )
        cHand1->detach();
    else
        cHand1->stepUnloaded();
}


/**
 Simulates:
 - attachment of cHand1
 - attached activity of cHand2
 .
 */
void Duo::stepFA()
{
    if ( active_ && prop->vulnerable )
        tryDeactivate();
    
    //we use cHand2->pos() first, because stepUnloaded() may detach cHand2
    cHand1->stepUnattached(simul(), cHand2->outerPos());

    if ( cHand2->checkDetachment() )
        cHand2->detach();
    else
        cHand2->stepUnloaded();
}


/**
 Simulates:
 - attached activity of cHand1
 - attached activity of cHand2
 .
 */
void Duo::stepAA()
{
    if ( active_ && prop->vulnerable )
        tryDeactivate();

    Vector f = force();
    real fn = f.norm();
    
    if ( cHand1->checkKramersDetachment(fn) )
        cHand1->detach();
    else
        cHand1->stepLoaded( f);
    
    if ( cHand2->checkKramersDetachment(fn) )
        cHand2->detach();
    else
        cHand2->stepLoaded(-f);
}


//------------------------------------------------------------------------------

void Duo::write(Outputter& out) const
{
    out.writeUInt8(active_);
    Couple::write(out);
}


void Duo::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    active_ = in.readUInt8();
    Couple::read(in, sim, tag);
}


