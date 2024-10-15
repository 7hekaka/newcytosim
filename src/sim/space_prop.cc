// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University
#include "space_prop.h"
#include "glossary.h"
#include "messages.h"
#include "simul_prop.h"
#include "cymdef.h"

#include "space.h"
#include "space_banana.h"
#include "space_capsule.h"
#include "space_cylinder.h"
#include "space_cylinderZ.h"
#include "space_cylinderP.h"
#include "space_bicylinder.h"
#include "space_dice.h"
#include "space_disc.h"
#include "space_ellipse.h"
#include "space_force.h"
#include "space_periodic.h"
#include "space_polygon.h"
#include "space_polygonZ.h"
#include "space_ring.h"
#include "space_sphere.h"
#include "space_square.h"
#include "space_strip.h"
#include "space_torus.h"

#if NEW_SPACES
#include "space_mesh.h"
#include "space_rotate.h"
#endif


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
 
 PARAMETER are specific to each shape but usually include 'length' or 'radius'.
 The last columun in the  list below gives the parameters associated with each shape.
 It is often possible to specify the radius or the diameter (but not both).
 
 The parameter values (DIMENSIONS) is a single REAL or a comma-separated list of REAL.
 
 List of known `shape`:
 
 SHAPE         | Class                | PARAMETERS
 --------------|----------------------|-------------------------------------
 `rectangle`   | SpaceSquare          | length = REAL, REAL, REAL;
 `sphere`      | SpaceSphere          | radius = REAL or diameter = REAL;
 `polygon`     | SpacePolygon         | file=FILE; height=REAL;
 `polygonZ`    | SpacePolygonZ        | file=FILE;
 `capsule`     | SpaceCapsule         | length=REAL; radius=REAL;
 `torus`       | SpaceTorus           | radius=REAL; width=REAL;
 `banana`      | SpaceBanana          | length=REAL; radius=REAL; curvature=REAL;
 `dice`        | SpaceDice            | length = REAL, REAL, REAL; edge=REAL
 `strip`       | SpaceStrip           | length = REAL, REAL, REAL;
 `periodic`    | SpacePeriodic        | length=REAL, REAL, REAL
 `ellipse`     | SpaceEllipse         | radius=REAL, REAL, REAL;
 `cylinder`    | SpaceCylinder        | length=REAL; radius=REAL;
 `cylinderZ`   | SpaceCylinderZ       | radius=REAL; bottom=REAL; top=REAL; edge=REAL
 `cylinderP`   | SpaceCylinderP       | length=REAL radius=REAL;
 `bicylinder`  | SpaceBicylinder      | radius=REAL;
 `ring`        | SpaceRing            | length=REAL; radius=REAL;
 `disc`        | SpaceDisc            | radius=REAL; bottom=REAL; top=REAL;
 `mesh`        | SpaceMesh            | file=FILE
 
 Example:
 
     set space cell
     {
         shape = sphere
     }
     new cell
     {
         radius = 5
     }
 
 */
Space * SpaceProp::newSpace() const
{
    std::string s = shape;

    if ( s=="circle" ) s = "sphere";
    if ( s=="rectangle" ) s = "square";
    if ( s=="spherocylinder" ) s = "capsule";
    if ( s=="semi_periodic" ) s = "strip";
    if ( s=="cube" ) s = "square";
#if ( DIM == 2 )
    if ( s=="cylinder" ) s = "square";
    if ( s=="cylinderP" ) s = "strip";
#elif ( DIM == 1 )
    if ( s=="cylinder" ) s = "square";
    if ( s=="cylinderP" ) s = "periodic";
    if ( s=="capsule" ) s = "square";
#endif

    if ( s=="square" )     return new SpaceSquare(this);
    if ( s=="sphere" )     return new SpaceSphere(this);
    if ( s=="polygon" )    return new SpacePolygon(this);
    if ( s=="polygonZ" )   return new SpacePolygonZ(this);
    if ( s=="capsule" )    return new SpaceCapsule(this);
    if ( s=="banana" )     return new SpaceBanana(this);
    if ( s=="torus" )      return new SpaceTorus(this);
    if ( s=="dice" )       return new SpaceDice(this);
    if ( s=="strip" )      return new SpaceStrip(this);
    if ( s=="periodic" )   return new SpacePeriodic(this);
    if ( s=="ellipse" )    return new SpaceEllipse(this);
#if ( DIM >= 3 )
    if ( s=="disc" )       return new SpaceDisc(this);
    if ( s=="cylinder" )   return new SpaceCylinder(this);
    if ( s=="cylinderZ" )  return new SpaceCylinderZ(this);
    if ( s=="cylinderP" )  return new SpaceCylinderP(this);
    if ( s=="bicylinder" ) return new SpaceBicylinder(this);
#endif
    if ( s=="ring" )       return new SpaceRing(this);
#if NEW_SPACES
    if ( s=="mesh" )       return new SpaceMesh(this);
    if ( s=="force" )      return new SpaceForce(this);
#endif
    
    //std::cerr << "Warning: unknown Space shape `"+s+"'\n";
    return nullptr;
}


Space * SpaceProp::newSpace(Glossary& opt) const
{
    Space * spc = newSpace();
    
    if ( spc )
    {
        //std::clog << "new Space `"+spc->prop->shape+"'\n";
#if BACKWARD_COMPATIBILITY < 50
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
    dimensions    = "";
    display       = "";
    display_fresh = false;
}


void SpaceProp::read(Glossary& glos)
{    
    if ( glos.set(shape, "shape") )
    {
#if BACKWARD_COMPATIBILITY < 50
        glos.set(dimensions, "dimensions");
    }
    if ( dimensions.empty() )
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

    if ( glos.set(display, "display") )
        display_fresh = true;
}


void SpaceProp::complete(Simul const&)
{
    if ( shape.empty() )
        throw InvalidParameter("space:shape must be defined");

}


void SpaceProp::write_values(std::ostream& os) const
{
    //write_value(os, "geometry",   geometry);
    write_value(os, "shape",      shape);
    write_value(os, "dimensions", dimensions);
    write_value(os, "display", "("+display+")");
}

