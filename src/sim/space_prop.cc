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
Space * SpaceProp::newSpace() const
{
    const std::string& s = SpaceProp::shape;
    
    if ( s=="rectangle" || s=="square" )       return new SpaceSquare(this);
    if ( s=="circle" || s=="sphere" )          return new SpaceSphere(this);
    if ( s=="polygon" )                        return new SpacePolygon(this);
    if ( s=="polygonZ" )                       return new SpacePolygonZ(this);
    if ( s=="capsule" || s=="spherocylinder" ) return new SpaceCapsule(this);
    if ( s=="banana" )                         return new SpaceBanana(this);
    if ( s=="torus" )                          return new SpaceTorus(this);
    if ( s=="dice" )                           return new SpaceDice(this);
    if ( s=="strip" || s=="half_periodic" )    return new SpaceStrip(this);
    if ( s=="periodic" )                       return new SpacePeriodic(this);
    if ( s=="ellipse" || s=="ellipsoid" )      return new SpaceEllipse(this);
#if ( DIM >= 3 )
    if ( s=="cubic" )                          return new SpaceSquare(this);
    if ( s=="cylinder" )                       return new SpaceCylinder(this);
    if ( s=="cylinderZ" )                      return new SpaceCylinderZ(this);
    if ( s=="cylinderP" )                      return new SpaceCylinderP(this);
#elif ( DIM == 2 )
    if ( s=="cylinder" )                       return new SpaceSquare(this);
    if ( s=="cylinderP" )                      return new SpaceStrip(this);
#else
    if ( s=="cylinder" )                       return new SpaceSquare(this);
    if ( s=="cylinderP" )                      return new SpacePeriodic(this);
#endif
    if ( s=="ring" )                           return new SpaceRing(this);
    if ( s=="tee" )                            return new SpaceTee(this);
#if NEW_SPACES
    if ( s=="mesh" )                           return new SpaceMesh(this);
    if ( s=="force" )                          return new SpaceForce(this);
    if ( s=="beads" )                          return new SpaceBeads(this);
#endif
#if NEW_DYNAMIC_SPACES
    if ( s=="lid" )                            return new SpaceLid(this);
    if ( s=="disc" )                           return new SpaceDisc(this);
    if ( s=="dynamic_sphere" )                 return new SpaceDynamicSphere(this);
    if ( s=="dynamic_ellipse" )                return new SpaceDynamicEllipse(this);
    // backward compatibility:
    if ( s=="contractile" )                    return new SpaceDynamicEllipse(this);
#endif
    
#if ( 1 )
    std::cerr << "Warning: substituting unbounded Space for unknown `"+s+"'\n";
    return new Space(this);
#endif
    return nullptr;
}


Space * SpaceProp::newSpace(Glossary& opt) const
{
    Space * spc = newSpace();
    
    if ( spc )
    {
#ifdef BACKWARD_COMPATIBILITY
        std::string str = dimensions;
        if ( str.size() || opt.set(str, "dimensions") )
        {
            std::stringstream iss(str);
            real len[8] = { 0 };
            int d = 0;
            while ( d < 8 )
            {
                real x = 0;
                iss >> x;
                if ( iss.fail() )
                    break;
                len[d++] = x;
            }
            if ( d > 0 )
            {
                spc->setLengths(len);
                return spc;
            }
        }
#endif
        // normal way to set the size:
        spc->resize(opt);
    }
    return spc;
}


//------------------------------------------------------------------------------

void SpaceProp::clear()
{
    shape         = "";
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
}


void SpaceProp::read(Glossary& glos)
{    
    if ( glos.set(shape, "shape") )
    {
#ifdef BACKWARD_COMPATIBILITY
        glos.set(dimensions, "dimensions");
    }
    else
    {
        std::string str;
        if ( glos.set(str, "geometry") )
        {
            std::stringstream iss(str);
            iss >> shape >> std::ws;
            std::getline(iss, dimensions);
            if ( dimensions.empty() )
                throw InvalidParameter("space:geometry should contains dimensions");
        }
#endif
    }

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

void SpaceProp::complete(Simul const& sim)
{
    if ( shape.empty() )
        throw InvalidParameter("space:shape must be defined");

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
}

//------------------------------------------------------------------------------

void SpaceProp::write_values(std::ostream& os) const
{
    //write_value(os, "geometry",   geometry);
    write_value(os, "shape",  shape);
#if NEW_DYNAMIC_SPACES
    write_value(os, "tension",    tension);
    write_value(os, "volume",     volume);
    write_value(os, "viscosity",  viscosity, viscosity_rot);
#endif
    write_value(os, "display",    "("+display+")");
}

