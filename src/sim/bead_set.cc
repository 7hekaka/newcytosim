// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "bead_set.h"
#include "bead_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"


Property* BeadSet::newProperty(const std::string& cat, const std::string& nom, Glossary&) const
{
    if ( cat == "bead" )
        return new BeadProp(cat, nom);
    return nullptr;
}


Object * BeadSet::newObject(const ObjectTag tag, PropertyID pid)
{
    if ( tag == Bead::TAG )
    {
        BeadProp * p = simul_.findProperty<BeadProp>("bead", pid);
        return new Bead(p, Vector(0,0,0), 0);
    }
    throw InvalidIO("Warning: unknown Bead tag `"+std::to_string(tag)+"'");
    return nullptr;
}

/**
 @ingroup NewObject

 A Bead is defined by its center and only the radius can be specified:

     new bead NAME
     {
       radius = VALUE, DEVIATION, MINIMUM
     }
 
 Variability around the mean radius is added if 'DEVIATION' and 'MINIMUM' are specified.
 All values must be positive.

 <h3> How to add Single </h3>

 Singles can only be attached at the center of the Bead:

     new bead NAME
     {
       radius = REAL
       attach = SINGLE [, SINGLE] ...
     }
 
 Where `SINGLE` is string containing at most 2 words: `[INTEGER] NAME`,
 where `INTEGER` specifies the number of Singles and `NAME` their name.
 
 For example if `grafted` is the name of a Single, one can use:
 
     new bead NAME
     {
       attach = 10 grafted
     }

 */

ObjectList BeadSet::newObjects(const std::string& name, Glossary& opt)
{
    ObjectList res;
    real rad = -1;
    size_t inx = 2;

    std::string var = "point1";
    if ( opt.has_key(var) )
    {
        if ( opt.value(var, 0) != "center" )
            throw InvalidParameter("position of `point1` must be `center'");
        opt.set(rad, var, 1);
    }
    else
    {
        inx = 0;
        var = "attach";
        opt.set(rad, "radius");
        
        // possibly add some variability in the radius:
        real dev = 0, inf = 0;
        if ( opt.set(dev, "radius", 1) && opt.set(inf, "radius", 2) )
        {
            real r;
            do
                r = rad + dev * RNG.gauss();
            while ( r < inf );
            rad = r;
        }
    }
    
    if ( rad <= 0 )
        throw InvalidParameter("bead:radius must be specified and > 0");

    BeadProp * p = simul_.findProperty<BeadProp>("bead", name);
    Bead * obj = new Bead(p, Vector(0,0,0), rad);
    
    res.push_back(obj);
    
    std::string str;
    // attach anchored Singles:
    while ( opt.set(str, var, inx++) )
        res.append(simul_.singles.makeWrists(obj, 0, 1, str));

    return res;
}


void BeadSet::write(Outputter& out) const
{
    if ( size() > 0 )
    {
        out.writeLine("\n#section "+title());
        writeObjects(out, pool_);
    }
}


void BeadSet::remove(Object * obj)
{
    ObjectSet::remove(obj);
    simul_.singles.removeWrists(obj);
}


void BeadSet::foldPositions(Modulo const* m) const
{
    for ( Bead * o=first(); o; o=o->next() )
        o->foldPosition(m);
}
