// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

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
    
    deactivation_rate = 0;
    deactivation_type = 0;
    activation = "off";
    vulnerable = true;
    activation_space = nullptr;
    activation_beads = nullptr;
}


void DuoProp::read(Glossary& glos)
{
    CoupleProp::read(glos);
    
    glos.set(deactivation_rate, "deactivation", "deactivation_rate");
    glos.set(deactivation_type, "deactivation", 1, {{"normal", 0}, {"delete", 1}});
    glos.set(activation, "activation", "activation_space");
    glos.set(vulnerable, "vulnerable");
}


void DuoProp::complete(Simul const& sim)
{
    CoupleProp::complete(sim);
    
    activation_space = sim.findSpace(activation);
#if ( 0 )
    if ( !activation_space )
    {
        activation_beads = sim.findSolidProp(activation);
    }
    
    if ( primed(sim) && !activation_space && !activation_beads )
        throw InvalidParameter("duo:activation not found!");
#endif
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
    write_value(os, "activation",  activation_space);
    write_value(os, "deactivation", deactivation_rate, deactivation_type);
    write_value(os, "vulnerable", vulnerable);
}

