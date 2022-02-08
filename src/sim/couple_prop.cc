// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "couple_prop.h"
#include "couple.h"
#include "couple_long.h"
#include "hand_prop.h"
#include "messages.h"
#include "glossary.h"
#include "hand_prop.h"
#include "simul_prop.h"
#include "simul.h"
#include <cmath>


/**
 This returns a new Couple if ( prop::length <= 0 ),
 or a CoupleLong if ( prop::length > 0 ).
 */
Couple * CoupleProp::newCouple(Glossary*) const
{
    //std::clog << "CoupleProp::newCouple" << '\n';
    if ( length > 0 )
        return new CoupleLong(this);
    
    return new Couple(this);
}

//------------------------------------------------------------------------------
#pragma mark -

void CoupleProp::clear()
{
    hand1             = "";
    hand2             = "";
    hand1_prop        = nullptr;
    hand2_prop        = nullptr;
    stiffness         = -1;
    length            = 0;
    diffusion         = 0;
    fast_diffusion    = false;
    fast_diffusion_nb = 0;
    trans_activated   = 0;
    min_loop          = 1;
    specificity       = BIND_ALWAYS;
    activity          = "diffuse";
    
    confine           = CONFINE_INSIDE;
    //confine_stiffness = 0;
    confine_space     = "first";
    confine_space_ptr = nullptr;
}


void CoupleProp::read(Glossary& glos)
{
    glos.set(hand1,           "hand1");
    glos.set(hand2,           "hand2");
    glos.set(stiffness,       "stiffness");
    glos.set(length,          "length");
    
    if ( glos.value_is("diffusion", 0, "fast") )
        fast_diffusion = 1;
    else
        glos.set(diffusion,   "diffusion");
    glos.set(fast_diffusion,  "fast_diffusion");
    glos.set(fast_diffusion_nb, "fast_diffusion", 1);
    
    glos.set(trans_activated, "trans_activated");
    // changed 'stiff' to 'min_loop' on 26.04.2020
    glos.set(min_loop, "min_loop", "stiff");

    glos.set(specificity, "specificity", {{"off",          BIND_ALWAYS},
#if BACKWARD_COMPATIBILITY < 50
                                          {"none",         BIND_ALWAYS},
#endif
                                          {"orthogonal",   BIND_ORTHOGONAL},
                                          {"parallel",     BIND_PARALLEL},
                                          {"not_parallel", BIND_NOT_PARALLEL},
                                          {"antiparallel", BIND_ANTIPARALLEL},
                                          {"not_antiparallel", BIND_NOT_ANTIPARALLEL}});
    
    glos.set(activity,        "activity");
    
    glos.set(confine,         "confine", {{"off",     CONFINE_OFF},
                                          {"on",      CONFINE_ON},
                                          {"none",    CONFINE_OFF},
                                          {"surface", CONFINE_ON},
                                          {"inside",  CONFINE_INSIDE}});
    
    real val;
    if ( glos.set(val, "confine", 1) && val > 0 )
        throw InvalidParameter(name()+":confine[1] is ignored");

    //glos.set(confine_stiffness, "confine", 1);
    glos.set(confine_space,   "confine", 2);

    //glos.set(confine_stiffness, "confine_stiffness");
    glos.set(confine_space,  "confine_space");

#if BACKWARD_COMPATIBILITY < 50
    if ( confine_space == "current" )
        confine_space = "last";
#endif
}


void CoupleProp::complete(Simul const& sim)
{
    if ( confine != CONFINE_OFF )
    {
        confine_space_ptr = sim.findSpace(confine_space);
        if ( confine_space_ptr )
            confine_space = confine_space_ptr->name();
        else
        {
            if ( sim.primed() )
                throw InvalidParameter(name()+":confine_space `"+confine_space+"' was not found");
            // this condition occur when the Property is created before the Space
        }
    }
    else
        confine_space_ptr = nullptr;

    if ( length < 0 )
        throw InvalidParameter(name()+":length must be >= 0");
    
    if ( diffusion < 0 )
        throw InvalidParameter(name()+":diffusion must be >= 0");

    /**
     We want for one degree of freedom to fulfill `var(dx) = 2 D time_step`
     And we use: dx = diffusion_dt * RNG.sreal()
     Since `sreal()` is uniformly distributed, its variance is 1/3,
     and we need `diffusion_dt^2 = 6 D time_step`
     */
    diffusion_dt = std::sqrt( 6.0 * diffusion * sim.time_step() * POOL_HAND_ATTACHMENT );

    if ( stiffness < 0 )
        throw InvalidParameter(name()+":stiffness must be specified and >= 0");
    
    if ( hand1.empty() )
        throw InvalidParameter(name()+":hand1 must be defined");
    hand1_prop = sim.findProperty<HandProp>("hand", hand1);
   
    if ( hand2.empty() )
        throw InvalidParameter(name()+":hand2 must be defined");
    hand2_prop = sim.findProperty<HandProp>("hand", hand2);
    
    if ( sim.primed() )
    {
        hand1_prop->checkStiffness(stiffness, length, 2, sim.prop.kT);
        /*
         If the length of a Couple (L) is longer than the attachment range of its hands,
         a Couple would place a pair of Fibers at a distance L, thus preventing further
         Couples from linking these two Fibers.
         In most cases, this is not desirable and physically inconsistent.
         */
        if ( length > hand1_prop->binding_range )
            Cytosim::warn << hand1_prop->name() << ":binding_range should be >= " << name() << ":length\n";

        if ( hand2_prop != hand1_prop )
        {
            hand2_prop->checkStiffness(stiffness, length, 2, sim.prop.kT);
        
            if ( length > hand2_prop->binding_range )
                Cytosim::warn << hand2_prop->name() << ":binding_range should be >= " << name() << ":length\n";
        }
    }
}


void CoupleProp::write_values(std::ostream& os) const
{
    write_value(os, "hand1",           hand1);
    write_value(os, "hand2",           hand2);
    write_value(os, "stiffness",       stiffness);
    write_value(os, "length",          length);
    write_value(os, "diffusion",       diffusion);
    write_value(os, "fast_diffusion",  fast_diffusion, fast_diffusion_nb);
    write_value(os, "trans_activated", trans_activated);
    write_value(os, "min_loop",        min_loop);
    write_value(os, "specificity",     specificity);
    write_value(os, "confine",         confine, 0, confine_space);
    write_value(os, "activity",        activity);
}


real CoupleProp::spaceVolume() const
{
    if ( !confine_space_ptr )
        throw InvalidParameter("no couple:confinement defined for `"+name()+"'");
    
    real res = confine_space_ptr->volume();
    
    if ( res <= 0 )
        throw InvalidParameter(name()+":confinement has null volume");
    
    return res;
}


