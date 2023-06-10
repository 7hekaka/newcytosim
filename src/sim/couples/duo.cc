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
    countdown_ = 0;
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
    
    // check activation
    if ( prop()->activation_space )
    {
        if ( prop()->activation_space->inside(cPos) )
            activate();
    }
    
    // activity
    if ( active_ )
    {
        assert_true(countdown_ > 0);
        // spontaneous de-activation:
        countdown_ -= prop()->deactivation_rate_dt;
        if ( countdown_ <= 0 )
        {
            deactivate();
            // test fraction of time when it is inactive:
            if ( RNG.test(-countdown_/prop()->deactivation_rate_dt) )
                return;
        }
        
        // hands may bind:
        if ( prop()->trans_activated || RNG.flip() )
            cHand1->stepUnattached(simul(), cPos);
        else
            cHand2->stepUnattached(simul(), cPos);
    }
    else if ( prop()->deactivation_mode )
    {
        delete(this);
        return;
    }
}


/**
 test for spontaneous de-activation
 */
void Duo::tryDeactivate()
{
    countdown_ -= prop()->deactivation_rate_dt;
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
            delete(this);
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
            delete(this);
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
            delete(this);
        return;
    }
}


#if NEW_DUO_HAS_TORQUE
void Duo::setInteractions(Meca& meca) const
{
    Interpolation const& pt1 = cHand1->interpolation();
    Interpolation const& pt2 = cHand2->interpolation();
    
    meca.addLink(pt1, pt2, prop()->stiffness);
    
#if ( DIM == 2 )
    if ( prop()->flip )
    {
        Vector2 dir = prop()->rest_dir;
        // flip the angle to match the current configuration of the bond
        sine = std::copysign(dir.YY, cross(pt1.diff(), pt2.diff()));
        dir.YY = sine;
        meca.addTorque(pt1, pt2, dir, prop()->angular_stiffness);
    }
    else
    {
        meca.addTorque(pt1, pt2, prop()->rest_dir, prop()->angular_stiffness);
    }
    //meca.addTorquePoliti(pt1, pt2, dir, prop()->angular_stiffness);
#elif ( DIM == 3 )
    meca.addTorque(pt1, pt2, prop()->rest_dir, prop()->angular_stiffness);
#endif
}
#endif

//------------------------------------------------------------------------------

void Duo::write(Outputter& out) const
{
    writeMarker(out, Couple::TAG_DUO);
    out.writeUInt8(active_);
    cHand1->writeHand(out);
    cHand2->writeHand(out);
    if ( !attached1() && !attached2() )
        out.writeFloats(cPos, DIM);
}


void Duo::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    if ( tag == Couple::TAG_DUO )
        active_ = in.readUInt8();
    else
        active_ = 1;
    if ( active_ )
        countdown_ = RNG.exponential();
    Couple::read(in, sim, tag);
}


