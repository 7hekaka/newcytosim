// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "duo.h"
#include "duo_prop.h"
#include "object_set.h"
#include "random.h"
#include "modulo.h"
#include "space.h"
#include "cymdef.h"
#include "meca.h"

//------------------------------------------------------------------------------

Duo::Duo(DuoProp const* p, Vector const& w)
: Couple(p, w), active_(0)
{
    nextDeact = -1;
    if ( p->fast_diffusion )
        throw InvalidParameter("`fast_diffusion` is incompatible with `activity=duo`");
}

Duo::~Duo()
{
}

//------------------------------------------------------------------------------

void Duo::activate()
{
    active_ = 1;
    nextDeact = RNG.exponential();
}

void Duo::deactivate()
{
    active_ = 0;
}

void Duo::recycle()
{
    objset()->remove(this);
    prop()->stocks.push(this);
}

//------------------------------------------------------------------------------

void Duo::stepFF()
{
    diffuse();
    
    // check activation
    if ( prop()->activation_space )
    {
        if ( prop()->activation_space->inside(cPos) )
            activate();
    }
    
    // activity
    if ( active_ )
    {
        assert_true(nextDeact >= 0);
        // spontaneous de-activation:
        nextDeact -= prop()->deactivation_rate_dt;
        if ( nextDeact < 0 )
        {
            deactivate();
            // test fraction of time when it is inactive:
            if ( RNG.test(-nextDeact/prop()->deactivation_rate_dt) )
                return;
        }
        
        // hands may bind:
        if ( prop()->trans_activated || RNG.flip() )
            cHand1->stepUnattached(simul(), cPos);
        else
            cHand2->stepUnattached(simul(), cPos);
    }
    else if ( prop()->deactivation_mode )
        recycle();
}


/**
 test for spontaneous de-activation
 */
void Duo::tryDeactivate()
{
    nextDeact -= prop()->deactivation_rate_dt;
    if ( nextDeact < 0 )
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
    if ( active_ && prop()->vulnerable )
        tryDeactivate();
    
    if ( active_ )
    {
        //we use cHand1->pos() first, because stepUnloaded() may detach cHand1
        cHand2->stepUnattached(simul(), cHand1->outerPos());
        
        if ( cHand1->checkDetachment() )
            cHand1->detach();
        else
            cHand1->stepUnloaded();
    }
    else
    {
        cHand1->detach();
        if ( prop()->deactivation_mode )
            recycle();
    }
}


/**
 Simulates:
 - attachment of cHand1
 - attached activity of cHand2
 .
 */
void Duo::stepFA()
{
    if ( active_ && prop()->vulnerable )
        tryDeactivate();
    
    if ( active_ )
    {
        //we use cHand2->pos() first, because stepUnloaded() may detach cHand2
        if ( !prop()->trans_activated )
            cHand1->stepUnattached(simul(), cHand2->outerPos());
        
        if ( cHand2->checkDetachment() )
            cHand2->detach();
        else
            cHand2->stepUnloaded();
    }
    else
    {
        cHand2->detach();
        if ( prop()->deactivation_mode )
            recycle();
        return;
    }
}


/**
 Simulates:
 - attached activity of cHand1
 - attached activity of cHand2
 .
 */
void Duo::stepAA()
{
    if ( active_ && prop()->vulnerable )
        tryDeactivate();

    if ( active_ )
    {
        Vector f = Couple::force();
        real mag = f.norm();
        
        if ( cHand1->checkKramersDetachment(mag) )
            cHand1->detach();
        else
            cHand1->stepLoaded( f);
        
        if ( cHand2->checkKramersDetachment(mag) )
            cHand2->detach();
        else
            cHand2->stepLoaded(-f);
    }
    else
    {
        cHand1->detach();
        cHand2->detach();
        if ( prop()->deactivation_mode )
            recycle();
        return;
    }
}


//------------------------------------------------------------------------------

void Duo::write(Outputter& out) const
{
    writeMarker(out, Couple::DUO_TAG);
    out.writeUInt8(active_);
    cHand1->writeHand(out);
    cHand2->writeHand(out);
    if ( !attached1() && !attached2() )
        out.writeFloats(cPos, DIM);
}


void Duo::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    if ( tag == Couple::DUO_TAG )
        active_ = in.readUInt8();
    else
        active_ = 1;
    if ( active_ )
        nextDeact = RNG.exponential();
    Couple::read(in, sim, tag);
}


