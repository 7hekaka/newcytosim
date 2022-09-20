// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include <cmath>
#include "space_prop.h"
#include "space_dynamic_prop.h"
#include "glossary.h"
#include "simul_prop.h"
#include "cymdef.h"
#include "simul_part.h"


// include spaces that use SpaceDynamicProp
#include "space_lid.h"
#include "space_dynamic_disc.h"
#include "space_dynamic_sphere.h"
#include "space_dynamic_ellipse.h"



/**
 @defgroup SpaceGroup Space and Geometry
 @ingroup ObjectGroup
 @ingroup NewObject
 @brief A Space defines a confined region
 
 A Space is created by specifying shape and dimensions:
 
     set space NAME
     {
        shape = SHAPE
     }
 
     new NAME
     {
        PARAMETER = DIMENSIONS
     }
 
 PARAMETER is usually 'length' or 'radius', but also 'height' or 'width'
 DIMENSIONS is a single REAL or a comma-separated list of REAL.
 
 List of Space with variable geometry:
 
 SHAPE              | Class                | PARAMETER        |
 -------------------|----------------------|--------------------
 `lid`              | SpaceLid             | width height
 `dynamic_disc`     | SpaceDynamicDisc     | radius
 `dynamic_sphere`   | SpaceDynamicSphere   | radius
 `dynamic_ellipse`  | SpaceDynamicEllipse  | DIM lengths
 
 Example:
 
     set space cell
     {
         shape = dynamic_disc
     }
     new cell
     {
         radius = 5
         viscosity = 1
     }
 
 */

//------------------------------------------------------------------------------

Space * SpaceDynamicProp::newSpace() const
{
    const std::string& s = SpaceProp::shape;
    if ( s=="lid" )                   return new SpaceLid(this);
    if ( s=="dynamic_disc" )          return new SpaceDynamicDisc(this);
    if ( s=="dynamic_sphere" )        return new SpaceDynamicSphere(this);
    if ( s=="dynamic_ellipse" )       return new SpaceDynamicEllipse(this);
#if BACKWARD_COMPATIBILITY < 50
    if ( s=="contractile" )           return new SpaceDynamicEllipse(this);
#endif
    //std::cerr << "Warning: unknown dynamic Space shape `"+s+"'\n";
    return nullptr;
}


void SpaceDynamicProp::clear()
{
    SpaceProp::clear();
    viscosity     = INFINITY;
    viscosity_rot = INFINITY;
    tension       = 0;
    volume        = 0;
    
    mobility_dt     = 0;
    mobility_rot_dt = 0;
}


void SpaceDynamicProp::read(Glossary& glos) 
{
    SpaceProp::read(glos);
    glos.set(viscosity,     "viscosity");
    glos.set(viscosity_rot, "viscosity", 1);
    glos.set(tension,       "tension");
    glos.set(volume,        "volume");
}


void SpaceDynamicProp::complete(Simul const& sim)
{
    SpaceProp::complete(sim);
    const real tau = time_step(sim);

    if ( viscosity > 0 )
        mobility_dt = tau / viscosity;
    else if ( primed(sim) )
        throw InvalidParameter("space:viscosity must be > 0");
    
    if ( viscosity_rot > 0 )
        mobility_rot_dt = tau / viscosity_rot;
    else if ( primed(sim) )
        throw InvalidParameter("space:viscosity[1] (rotational viscosity) must be > 0");

    if ( tension < 0 )
        throw InvalidParameter("tension must be >= 0");
}


void SpaceDynamicProp::write_values(std::ostream& os) const
{
    SpaceProp::write_values(os);
    write_value(os, "viscosity", viscosity, viscosity_rot);
    write_value(os, "tension",   tension);
    write_value(os, "volume",    volume);
}

