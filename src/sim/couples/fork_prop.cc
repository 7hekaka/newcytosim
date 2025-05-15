// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "dim.h"
#include "exceptions.h"
#include "glossary.h"
#include "simul_prop.h"
#include "fork_prop.h"
#include "fork.h"


Couple * ForkProp::newCouple() const
{
    //std::clog << "ForkProp::newCouple" << '\n';
    return new Fork(this);
}


void ForkProp::clear()
{
    CoupleProp::clear();
    rest_angle = 0;
    rest_dir.set(1, 0);
    angular_stiffness = 0;
    flip = true;
}


void ForkProp::read(Glossary& glos)
{
    CoupleProp::read(glos);
    
    if ( glos.has_key("torque") )
    {
        glos.set(angular_stiffness, "torque");
        glos.set(rest_angle, "torque", 1);
    }
    else
    {
        glos.set(angular_stiffness, "angular_stiffness");
        glos.set(rest_angle, "rest_angle", "angle");
    }
    glos.set(flip, "flip");
}


void ForkProp::complete(Simul const& sim)
{
    CoupleProp::complete(sim);
    
    rest_dir.XX = std::cos(rest_angle);
    rest_dir.YY = std::sin(rest_angle);
    if ( DIM == 3 ) rest_dir.YY = abs_real(rest_dir.YY);
#if ( 0 )
    if ( rest_angle < 0 || rest_dir.YY < 0 )
        throw InvalidParameter("The equilibrium angle should be defined in [0, pi]");
#endif

    if ( angular_stiffness < 0 )
        throw InvalidParameter("The angular stiffness, torque[0] must be >= 0");
}


void ForkProp::write_values(std::ostream& os) const
{
    CoupleProp::write_values(os);
    write_value(os, "torque", angular_stiffness, rest_angle);
    write_value(os, "flip", flip);
}

