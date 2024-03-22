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
#if NEW_DUO_HAS_TORQUE
    rest_angle = 0;
    rest_dir.set(1, 0);
    angular_stiffness = 0;
    flip = true;
#endif
}

void DuoProp::read(Glossary& glos)
{
    CoupleProp::read(glos);
    
    glos.set(deactivation_rate, "deactivation", "deactivation_rate");
    glos.set(deactivation_mode, "deactivation", 1, {{"normal", 0}, {"delete", 1}});
    glos.set(activation, "activation", "activation_space");
    glos.set(vulnerable, "vulnerable");
#if NEW_DUO_HAS_TORQUE
    if ( glos.has_key("torque") )
    {
        glos.set(angular_stiffness, "torque");
        glos.set(rest_angle, "torque", 1);
    }
    else
    {
        glos.set(angular_stiffness, "angular_stiffness");
        glos.set(rest_angle, "angle", 0);
    }
    glos.set(flip, "flip");
#endif
}


void DuoProp::splash(std::ostream& os) const
{
    std::ostringstream oss;
    real L = std::sqrt(diffusion/deactivation_rate);
    oss << std::setw(10) << name() << ": deactivation_rate " << deactivation_rate;
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
    
#if NEW_DUO_HAS_TORQUE
    rest_dir.XX = std::cos(rest_angle);
    rest_dir.YY = std::sin(rest_angle);
    if ( DIM == 3 ) rest_dir.YY = abs_real(rest_dir.YY);
    if ( angular_stiffness < 0 )
        throw InvalidParameter("The angular stiffness, torque[0] must be >= 0");
#endif

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
#if NEW_DUO_HAS_TORQUE
    write_value(os, "torque", angular_stiffness, rest_angle);
    write_value(os, "flip", flip);
#endif
}

