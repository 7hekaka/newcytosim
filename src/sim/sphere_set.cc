// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "sphere_set.h"
#include "sphere_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"


void SphereSet::remove(Object * obj)
{
    ObjectSet::remove(obj);
    simul_.singles.deleteWrists(obj);
}


Property* SphereSet::newProperty(const std::string& cat, const std::string& nom, Glossary&) const
{
    if ( cat == "sphere" )
        return new SphereProp(nom);
    return nullptr;
}


Object * SphereSet::newObject(const ObjectTag tag, PropertyID pid)
{
    if ( tag == Sphere::TAG )
    {
        SphereProp * p = simul_.findProperty<SphereProp>("sphere", pid);
        return new Sphere(p);
    }
    throw InvalidIO("Warning: unknown Sphere tag `"+std::to_string(tag)+"'");
    return nullptr;
}

/**
 @copydetails Sphere::build
 */
ObjectList SphereSet::newObjects(Property const* p, Glossary& opt)
{
    SphereProp const* pp = static_cast<SphereProp const*>(p);
    // set radius if provided as argument
    real rad = -1;
    if ( !opt.set(rad, "radius" ) || rad <= 0 )
        throw InvalidParameter("parameter `radius` should be specified and > 0");
    
    // possibly add some variability
    real dev = 0, inf = 0;
    if ( opt.set(dev, "radius", 1) && opt.set(inf, "radius", 2) )
    {
        real r;
        do
            r = rad + dev * RNG.gauss();
        while ( r < inf );
        rad = r;
    }
    
    Sphere * obj = new Sphere(pp, rad);
    return obj->build(opt, simul_);
}


void SphereSet::writeSet(Outputter& out) const
{
    if ( size() > 0 )
    {
        out.write("\n#section "+title());
        writePool(out, pool_);
    }
}


void SphereSet::foldPositions(Modulo const* m) const
{
    for ( Sphere * o=SphereSet::first(); o; o=o->next() )
        o->foldPosition(m);
}

