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
Couple * DuoProp::newCouple() const
{
    //std::clog << "DuoProp::newCouple" << '\n';
    if ( length > 0 )
        return new DuoLong(this);
    else
        return new Duo(this);
}


void DuoProp::clear()
{
    CoupleProp::clear();
    
    deactivation_rate = 0;
    deactivation_mode = 0;
    activation = "off";
    vulnerable = true;
    activation_space = nullptr;
}

void DuoProp::read(Glossary& glos)
{
    CoupleProp::read(glos);
    
    glos.set(deactivation_rate, "deactivation", "deactivation_rate");
    glos.set(deactivation_mode, "deactivation", 1, {{"normal", 0}, {"delete", 1}});
    glos.set(activation, "activation", "activation_space");
    glos.set(vulnerable, "vulnerable");
}


void DuoProp::splash(std::ostream& os) const
{
    std::ostringstream oss;
    real L = std::sqrt(diffusion/deactivation_rate);
    oss << std::setw(10) << name();
    oss << ": deactivation_rate " << deactivation_rate;
    oss << "  traveled_distance " << L << "\n";
    if ( oss.str() != splashed )
    {
        splashed = oss.str();
        os << splashed;
    }
}


void DuoProp::complete(Simul const& sim)
{
    CoupleProp::complete(sim);
    
    activation_space = sim.findSpace(activation);
    
    if ( deactivation_rate < 0 )
        throw InvalidParameter("deactivation_rate should be >= 0");
    
    deactivation_rate_dt = deactivation_rate * time_step(sim) * POOL_UNATTACHED;

    /// print predicted decay distance in verbose mode:
    if ( primed(sim) && sim.prop.verbose )
        splash(std::clog);
}


void DuoProp::write_values(std::ostream& os) const
{
    CoupleProp::write_values(os);
    write_value(os, "activation",  activation_space);
    write_value(os, "deactivation", deactivation_rate, deactivation_mode);
    write_value(os, "vulnerable", vulnerable);
}

