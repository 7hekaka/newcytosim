// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University
#include "single_prop.h"
#include "glossary.h"
#include "messages.h"

#include "simul_prop.h"
#include "hand_prop.h"
#include "single.h"
#include "wrist.h"
#include "wrist_long.h"
#include "picket.h"
#include "picket_long.h"
#include "simul.h"

/**
 @defgroup SingleGroup Single and related
 @ingroup ObjectGroup
 @ingroup NewObject
 @brief A Single contains one Hand, and can thus bind to one Fiber.

 List of classes accessible by specifying single:activity:
 
 `activity`          | Class             | Description                               |
 --------------------|-------------------|--------------------------------------------
 `diffuse` (default) | Single            | a single Hand that is mobile (default)
 `fixed`             | Picket PicketLong | a single Hand anchored at a fixed position

 The Single will actually move only if its diffusion coefficient is set and > 0.
 Another class Wrist is used automatically to anchor a Single to a Mecable.
 
 Example:

     set single grafted
     {
       hand = kinesin
       stiffness = 100
       activity = fixed
     } 

 */
Single * SingleProp::newSingle() const
{
    //std::clog << "SingleProp::newSingle" << '\n';
    if ( activity == "fixed" )
    {
        if ( length > 0 )
            return new PicketLong(this);
        else
            return new Picket(this);
    }
    else if ( activity == "diffuse" )
    {
        return new Single(this);
    }
#if ( 0 )
    throw InvalidParameter("unknown single:activity `"+activity+"'");
#else
    // try to proceed anyhow:
    std::cerr << "WARNING: unknown single:activity `"+activity+"'\n";
#endif
    return new Single(this);
}

/**
 Create Wrist anchored to a Mecable vertex
 */
Wrist * SingleProp::newWrist(Mecable const* mec, const unsigned point) const
{
    //std::clog << "SingleProp::newWrist(length=" << length << ")\n";
    if ( length > 0 )
        return new WristLong(this, mec, point);
    else
        return new Wrist(this, mec, point);
}

/**
 Create Wrist anchored on points interpolated from the Mecable's vertices
 */
Wrist * SingleProp::newWrist(Mecable const* mec, const unsigned ref, Vector const& vec) const
{
    //std::clog << "SingleProp::newWrist(length=" << length << ")\n";
    Wrist * w;
    if ( length > 0 )
        w = new WristLong(this, mec, ref);
    else
        w = new Wrist(this, mec, ref);
    w->rebase(mec, ref, vec);
    return w;
}

//------------------------------------------------------------------------------
#pragma mark -

void SingleProp::clear()
{
    hand           = "";
    hand_prop      = nullptr;
    stiffness      = 0;
    length         = 0;
    diffusion      = 0;
    fast_diffusion = 0;
    fast_reservoir = 0;
    save_unbound  = ~0U;
#if NEW_MOBILE_SINGLE
    speed.reset();
#endif
    activity      = "diffuse";
    confine       = CONFINE_INSIDE;
    //confine_stiff = 0;
    confine_spec = "first";
    confine_space = nullptr;
    uni_counts = 0;
}


void SingleProp::read(Glossary& glos)
{
    glos.set(hand,           "hand");
    glos.set(stiffness,      "stiffness");
    glos.set(length,         "length");
    if ( glos.value_is("diffusion", 0, "fast") )
    {
        diffusion = 100;
        fast_diffusion = 1;
    }
    else
        glos.set(diffusion,  "diffusion");
    glos.set(fast_diffusion, "fast_diffusion");
    glos.set(fast_reservoir, "fast_diffusion", 1);
    glos.set(save_unbound, "save_unbound");
#if NEW_MOBILE_SINGLE
    glos.set(speed,          "speed");
#endif
    if ( glos.set(activity, "activity") )
    {
        // change the default value for confinement
        if ( activity == "fixed" )
        {
            if ( fast_diffusion )
                fast_diffusion = -1;
            confine = CONFINE_OFF;
        }
    }
    
    glos.set(confine, "confine",
        {{"off",    CONFINE_OFF},
        {"on",      CONFINE_ON},
        {"none",    CONFINE_OFF},
        {"surface", CONFINE_ON},
        {"inside",  CONFINE_INSIDE}});
    
    real val;
    if ( glos.set(val, "confine", 1) && val > 0 )
        throw InvalidParameter(name()+":confine[1] is ignored");
    
    glos.set(confine_spec, "confine", 2);
    
#if BACKWARD_COMPATIBILITY < 50
    if ( confine_spec == "current" )
        confine_spec = "last";
#endif
}


void SingleProp::complete(Simul const& sim)
{
    confine_space = sim.findSpace(confine_spec);
    if ( confine != CONFINE_OFF )
    {
        if ( activity=="fixed" )
            throw InvalidParameter(name()+":confine is ignored since activity=fixed");
        if ( confine_space )
        {
            if ( confine_spec.empty() )
                confine_spec = sim.spaces.nameObject(confine_space);
        }
        else
        {
            // this condition may occur when the Property is created before the Space
            if ( primed(sim) )
                throw InvalidParameter(name()+":confine_space `"+confine_spec+"' was not found");
        }
    }

    if ( hand.empty() )
        throw InvalidParameter(name()+":hand must be defined");
    hand_prop = sim.findProperty<HandProp>("hand", hand);
    
    if ( !hand_prop )
        throw InvalidParameter("unknown single:hand '"+hand+"'");

    if ( diffusion < 0 )
        throw InvalidParameter(name()+":diffusion must be >= 0");

    /**
     We want for one degree of freedom to fulfill `var(dx) = 2 D time_step`
     And we use: dx = diffusion_dt * RNG.sreal()
     Since `sreal()` is uniformly distributed, its variance is 1/3,
     and we need `diffusion_dt^2 = 6 D time_step`
     */
    diffusion_dt = std::sqrt(6.0 * diffusion * time_step(sim) * POOL_UNATTACHED);
#if NEW_MOBILE_SINGLE
    speed_dt = speed * time_step(sim);
#endif
    
    if ( stiffness < 0 )
        throw InvalidParameter(name()+":stiffness must be >= 0");

    if ( length < 0 )
        throw InvalidParameter(name()+":length must be >= 0");

    if ( primed(sim) && stiffness > 0 )
    {
        hand_prop->checkStiffness(stiffness, length, 1, boltzmann(sim));
        
        /*
         If the length of a Single (L) is longer than the attachment range of its hands,
         a Couple would place a pair of Fibers at a distance L, thus preventing further
         Singles from linking these two Fibers.
         In most cases, this is not desirable and physically inconsistent.
         */
        if ( length > hand_prop->binding_range )
            throw InvalidParameter(hand_prop->name()+":binding_range must be >= "+name()+":length");
            //Cytosim::warn("Attachment cannot occur because single:length > hand:binding_range\n");
    }
}


//------------------------------------------------------------------------------

void SingleProp::write_values(std::ostream& os) const
{
    write_value(os, "hand",           hand);
    write_value(os, "stiffness",      stiffness);
    write_value(os, "length",         length);
    write_value(os, "diffusion",      diffusion);
    write_value(os, "fast_diffusion", fast_diffusion, fast_reservoir);
#if NEW_MOBILE_SINGLE
    write_value(os, "speed",          speed);
#endif
    write_value(os, "confine",        confine, 0, confine_spec);
    write_value(os, "activity",       activity);
}


//------------------------------------------------------------------------------

real SingleProp::spaceVolume() const
{
    if ( !confine_space )
        throw InvalidParameter("no single:confinement defined for `"+name()+"'");
    
    real res = confine_space->volume();
    
    if ( res <= 0 )
        throw InvalidParameter(name()+":confinement has null volume");
    
    return res;
}
