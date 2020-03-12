// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <cmath>
#include "space_prop.h"
#include "property_list.h"
#include "glossary.h"
#include "simul_prop.h"
#include "simul.h"
#include "sim.h"

// To create new spaces
#include "space_lid.h"
#include "space_disc.h"
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
 
 List of known `shape` with variable geometry:
 
 GEOMETRY           | Class                | PARAMETER        |
 -------------------|----------------------|--------------------
 `lid`              | SpaceLid             | width height
 `disc`             | SpaceDisc            | radius
 `dynamic_sphere`   | SpaceDynamicSphere   | radius
 
 Example:
 
     set space cell
     {
         shape = disc
     }
     new cell
     {
         radius = 5
         viscosity = 1
     }
 
 */

//------------------------------------------------------------------------------
void DynamicSpaceProp::complete(Simul const& sim) 
{
	SpaceProp::complete(sim);
	
    if ( viscosity > 0 )
        mobility_dt = sim.time_step() / viscosity;
    else if ( sim.ready() )
        throw InvalidParameter("space:viscosity must be > 0");
    
    if ( viscosity_rot > 0 )
        mobility_rot_dt = sim.time_step() / viscosity_rot;
    else if ( sim.ready() )
        throw InvalidParameter("space:viscosity[1] (rotational viscosity) must be > 0");
        
}

Space * DynamicSpaceProp::newSpace() const
{
    const std::string& s = SpaceProp::shape;
   
    if ( s=="lid" )                            return new SpaceLid(this);
    if ( s=="disc" )                           return new SpaceDisc(this);
    if ( s=="dynamic_sphere" )                 return new SpaceDynamicSphere(this);
    //std::cerr << "Warning: unknown Space shape `"+s+"'\n";
    return nullptr;
}

void DynamicSpaceProp::clear()
{
    SpaceProp::clear();
    viscosity     = INFINITY;
    viscosity_rot = INFINITY;
    mobility_dt   = 0;
    mobility_rot_dt = 0;
}

void DynamicSpaceProp::read(Glossary& glos) 
{
    SpaceProp::read(glos);
    glos.set(viscosity,     "viscosity");
    glos.set(viscosity_rot, "viscosity", 1);
}

void DynamicSpaceProp::write_values(std::ostream& os) const
{
    SpaceProp::write_values(os);
    write_value(os, "viscosity",  viscosity, viscosity_rot);
}

