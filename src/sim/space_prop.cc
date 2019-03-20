// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_prop.h"
#include "filepath.h"
#include "glossary.h"
#include "property_list.h"
#include "simul_prop.h"
#include "simul.h"
#include "sim.h"

#include "space.h"
#include "space_force.h"
#include "space_square.h"
#include "space_sphere.h"
#include "space_polygon.h"
#include "space_polygonZ.h"
#include "space_capsule.h"
#include "space_banana.h"
#include "space_torus.h"
#include "space_dice.h"
#include "space_strip.h"
#include "space_periodic.h"
#include "space_ellipse.h"
#include "space_cylinder.h"
#include "space_cylinderZ.h"
#include "space_cylinderP.h"
#include "space_ring.h"
#include "space_beads.h"
#include "space_tee.h"

#if NEW_SPACES
#include "space_mesh.h"
#include "space_rotate.h"
#endif

#if NEW_DYNAMIC_SPACES
#include "space_lid.h"
#include "space_disc.h"
#include "space_dynamic_sphere.h"
#include "space_dynamic_ellipse.h"
#endif

/**
 @defgroup SpaceGroup Space and Geometry
 @ingroup ObjectGroup
 @ingroup NewObject
 @brief A Space defines a confined region
 
 A Space is created by specifying a geometry:
 
     set space NAME
     {
        geometry = GEOMETRY DIMENSIONS
     }
 
 
 DIMENSIONS is usually a list of numbers.
 
 List of known `geometry`:
 
 GEOMETRY      | Class                | DIMENSIONS                             |
 --------------|----------------------|-----------------------------------------
 `rectangle`   | SpaceSquare          | sizeX sizeY sizeZ
 `sphere`      | SpaceSphere          | radius
 `polygon`     | SpacePolygon         | file_name height
 `polygonZ`    | SpacePolygonZ        | file_name
 `capsule`     | SpaceCapsule         | half_length radius
 `torus`       | SpaceTorus           | radius thickness
 `banana`      | SpaceBanana          | total_length width radius_of_curvature
 `dice`        | SpaceDice            | sizeX sizeY sizeZ radius
 `strip`       | SpaceStrip           | sizeX sizeY sizeZ
 `periodic`    | SpacePeriodic        | sizeX sizeY sizeZ
 `ellipse`     | SpaceEllipse         | sizeX sizeY sizeZ
 `cylinder`    | SpaceCylinder        | half_length radius
 `cylinderZ`   | SpaceCylinderZ       | radius bottom top
 `cylinderP`   | SpaceCylinderP       | half_length radius
 `ring`        | SpaceRing            | half_length radius
 `tee`         | SpaceTee             | half_length radius arm_position arm_length
 `mesh`        | SpaceMesh            | mesh file

 Dynamic Space with variable geometry:
 
 GEOMETRY           | Class                | DIMENSIONS        |
 -------------------|----------------------|--------------------
 `lid`              | SpaceLid             | width height
 `disc`             | SpaceDisc            | radius
 `dynamic_sphere`   | SpaceDynamicSphere   | radius
 `dynamic_ellipse`  | SpaceDynamicEllipse  | sizeX sizeY sizeZ
 
 Example:
 
     set space cell
     {
          geometry = sphere 5
     }
 
 */
Space * SpaceProp::newSpace(SpaceProp const* sp, Glossary& opt)
{
    const std::string s = sp->shape;
    
    if ( s=="rectangle" || s=="square" )       return new SpaceSquare(sp);
    if ( s=="circle" || s=="sphere" )          return new SpaceSphere(sp);
    if ( s=="polygon" )                        return new SpacePolygon(sp, opt);
    if ( s=="polygonZ" )                       return new SpacePolygonZ(sp, opt);
    if ( s=="capsule" || s=="spherocylinder" ) return new SpaceCapsule(sp);
    if ( s=="banana" )                         return new SpaceBanana(sp);
    if ( s=="torus" )                          return new SpaceTorus(sp);
    if ( s=="dice" )                           return new SpaceDice(sp);
    if ( s=="strip" || s=="half_periodic" )    return new SpaceStrip(sp);
    if ( s=="periodic" )                       return new SpacePeriodic(sp);
    if ( s=="ellipse" || s=="ellipsoid" )      return new SpaceEllipse(sp);
#if ( DIM >= 3 )
    if ( s=="cubic" )                          return new SpaceSquare(sp);
    if ( s=="cylinder" )                       return new SpaceCylinder(sp);
    if ( s=="cylinderZ" )                      return new SpaceCylinderZ(sp);
    if ( s=="cylinderP" )                      return new SpaceCylinderP(sp);
#elif ( DIM == 2 )
    if ( s=="cylinder" )                       return new SpaceSquare(sp);
    if ( s=="cylinderP" )                      return new SpaceStrip(sp);
#else
    if ( s=="cylinder" )                       return new SpaceSquare(sp);
    if ( s=="cylinderP" )                      return new SpacePeriodic(sp);
#endif
    if ( s=="ring" )                           return new SpaceRing(sp);
    if ( s=="tee" )                            return new SpaceTee(sp);
#if NEW_SPACES
    if ( s=="mesh" )                           return new SpaceMesh(sp, opt);
    if ( s=="force" )                          return new SpaceForce(sp, opt);
    if ( s=="beads" )                          return new SpaceBeads(sp);
#endif
#if NEW_DYNAMIC_SPACES
    if ( s=="lid" )                            return new SpaceLid(sp);
    if ( s=="disc" )                           return new SpaceDisc(sp);
    if ( s=="dynamic_sphere" )                 return new SpaceDynamicSphere(sp);
    if ( s=="dynamic_ellipse" )                return new SpaceDynamicEllipse(sp);
    // backward compatibility:
    if ( s=="contractile" )                    return new SpaceDynamicEllipse(sp);
#endif
    
#if ( 1 )
    std::cerr << "INCIDENT: using unbounded Space instead of unknown class `"+s+"'\n";
    return new Space(sp);
#endif
    return nullptr;
}


Space * SpaceProp::newSpace(Glossary& opt) const
{
    Space * spc = newSpace(this, opt);
    
    if ( spc )
    {
        // set dimensions:
        if ( dimensions.length() )
            spc->readLengths(dimensions);
        
        std::string dim;
        if ( opt.set(dim, "dimensions") )
            spc->readLengths(dim);
    }
    return spc;
}


//------------------------------------------------------------------------------

void SpaceProp::clear()
{
    geometry      = "";
    shape         = "undefined";
    dimensions    = "";
    shape_spec    = "";
    display       = "";
    display_fresh = false;
    
#if NEW_DYNAMIC_SPACES
    tension       = 0;
    volume        = 0;
    viscosity     = INFINITY;
    viscosity_rot = INFINITY;
    mobility_dt   = 0;
    mobility_rot_dt = 0;
#endif
    
    dimensions_old = "";
}

void SpaceProp::read(Glossary& glos)
{    
    glos.set(shape,        "shape");
    glos.set(shape_spec,   "shape", 1);
    glos.set(dimensions,   "dimension") || glos.set(dimensions, "dimensions");
#ifdef BACKWARD_COMPATIBILITY
    glos.set(dimensions,   "spec");  // format 36
#endif
    glos.set(geometry,     "geometry");
    
#if NEW_DYNAMIC_SPACES
    glos.set(tension,       "tension");
    glos.set(volume,        "volume");
    glos.set(viscosity,     "viscosity");
    glos.set(viscosity_rot, "viscosity", 1);
#endif
    
    if ( glos.set(display, "display") )
        display_fresh = true;
}

//------------------------------------------------------------------------------

void SpaceProp::complete()
{
    if ( !geometry.empty() )
    {
        std::istringstream iss(geometry);
        iss >> shape;
        
        if ( iss.fail() )
            throw InvalidParameter("invalid geometry `"+geometry+"' for Space");
        
        int c;
        while ( isspace(iss.peek()) )
            iss.get();
        
        c = iss.peek();
        if ( !isdigit(c) && c != '+' && c != '-' )
            iss >> shape_spec;
        
        while ( isspace(iss.peek()) )
            iss.get();

        // get remaining characters as a whole:
        if ( iss.good() )
            dimensions = geometry.substr(iss.tellg());
    }
}


void SpaceProp::complete(Simul const& sim)
{
    complete();
    
#if NEW_DYNAMIC_SPACES
    if ( viscosity > 0 )
        mobility_dt = sim.prop->time_step / viscosity;
    else if ( sim.ready() )
        throw InvalidParameter("space:viscosity must be > 0");
    
    if ( viscosity_rot > 0 )
        mobility_rot_dt = sim.prop->time_step / viscosity_rot;
    else if ( sim.ready() )
        throw InvalidParameter("space:viscosity[1] (rotational viscosity) must be > 0");
#endif

    /*
     If the dimensions have changed, update any Space with this property.
     This is necessary to make 'change space:dimension' work.
     */
    if ( dimensions != dimensions_old )
    {
        for ( Space * spc = sim.spaces.first(); spc; spc=spc->next() )
        {
            if ( spc->prop == this )
            {
                spc->readLengths(dimensions);
                // allow Simul to update:
                if ( spc == sim.space() )
                    const_cast<Simul&>(sim).changeSpace(spc);
            }
        }
        dimensions_old = dimensions;
    }
}

//------------------------------------------------------------------------------

void SpaceProp::write_values(std::ostream& os) const
{
    //write_value(os, "geometry",   geometry);
    if ( shape_spec.empty() )
        write_value(os, "shape",  shape);
    else
        write_value(os, "shape",  shape, shape_spec);
    write_value(os, "dimensions", dimensions);
#if NEW_DYNAMIC_SPACES
    write_value(os, "tension",    tension);
    write_value(os, "volume",     volume);
    write_value(os, "viscosity",  viscosity, viscosity_rot);
#endif
    write_value(os, "display",    "("+display+")");
}



