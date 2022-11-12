// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "messages.h"
#include "exceptions.h"
#include "glossary.h"
#include "hand_prop.h"
#include "duo_prop.h"
#include "duo.h"
#include "duo_long.h"
#include "simul.h"


/**
 returns a Duo if ( length <= 0 ),
 or a DuoLong if ( length > 0 )
 */
Couple * DuoProp::newCouple(Glossary * opt) const
{
    Duo * res = nullptr;
    //std::clog << "DuoProp::newCouple" << '\n';
    if ( length > 0 )
        res = new DuoLong(this);
    else
        res = new Duo(this);
    
    int a;
    if ( opt  &&  opt->set(a, "active") )
    {
        if ( a )
            res->activate();
    }
    
    return res;
}


void DuoProp::clear()
{
    CoupleProp::clear();
    
    deactivation_rate    = 0;
    activation_space     = "off";
    activation_space_ptr = nullptr;
    vulnerable           = true;
}


void DuoProp::read(Glossary& glos)
{
    CoupleProp::read(glos);
    
    glos.set(deactivation_rate, "deactivation_rate");
    glos.set(activation_space,  "activation_space");
    glos.set(vulnerable, "vulnerable");
}


void DuoProp::complete(Simul const& sim)
{
    CoupleProp::complete(sim);
    
    activation_space_ptr = sim.findSpace(activation_space);
    
    if ( primed(sim)  &&  !activation_space_ptr )
        throw InvalidParameter("duo:activation_space not found!");

    if ( deactivation_rate < 0 )
        throw InvalidParameter("deactivation_rate should be >= 0");
    
    deactivation_rate_dt = deactivation_rate * time_step(sim) * POOL_UNATTACHED;
    
    /// print predicted decay distance in verbose mode:
    if ( primed(sim) && sim.prop.verbose )
    {
        real L = std::sqrt(diffusion / deactivation_rate);
        std::clog << name() << ":deactivation_rate " << deactivation_rate;
        std::clog << "  length " << L << "\n";
    }
}


void DuoProp::write_values(std::ostream& os) const
{
    CoupleProp::write_values(os);
    write_value(os, "activation_space",  activation_space);
    write_value(os, "deactivation_rate", deactivation_rate);
    write_value(os, "vulnerable", vulnerable);
}

