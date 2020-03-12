// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <cmath>
#include "space_dynamic_ellipse_prop.h"
#include "space_dynamic_ellipse.h"
#include "dynamic_space_prop.h"
#include "property_list.h"
#include "glossary.h"


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
 `dynamic_ellipse`  | SpaceDynamicEllipse  | length
 
 Example:
 
     set space cell
     {
         shape = dynamic_ellipse
     }
     new cell
     {
         dimensions = 10,5,1
         tension = 10
         viscosity = 1
     }
 
 */

//------------------------------------------------------------------------------


Space * SpaceDynamicEllipseProp::newSpace() const
{
	return new SpaceDynamicEllipse(this);
}


Space * SpaceDynamicEllipseProp::newSpace(Glossary& opt) const
{
    Space * spc = newSpace();
    
    if ( spc )
    {
        // normal way to set the size:
        spc->resize(opt);
    }
    return spc;
}


void SpaceDynamicEllipseProp::clear()
{
    tension = 0 ;
    volume  = 0  ;
	DynamicSpaceProp::clear();
}

void SpaceDynamicEllipseProp::read(Glossary& glos)
{
    DynamicSpaceProp::read(glos);
	glos.set(tension, "tension");

}



void SpaceDynamicEllipseProp::complete(Simul const& sim) 
{
	DynamicSpaceProp::complete(sim);
	
	if	(tension < 0)
		throw InvalidParameter("tension must be positive");
}


void SpaceDynamicEllipseProp::write_values(std::ostream& os) const
{
    DynamicSpaceProp::write_values(os);
	write_value(os, "tension", tension);
}