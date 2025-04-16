// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "primitives.h"
#include "bead_set.h"
#include "bead_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"


void BeadSet::steps()
{
#if ( 0 )
    Bead * obj = first();
    while ( obj )
    {
        Bead * nxt = obj->next();
        obj->step();
        obj = nxt;
    }
#endif
    if ( size() > 1 ) shuffle();
}


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

ObjectList BeadSet::newObjects(Property const* p, Glossary& opt)
{
    BeadProp const* pp = static_cast<BeadProp const*>(p);
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
        std::string str;
        if ( opt.set_block(str, '[', "radius") ) // some code specified
        {
            // get a range of radius
            float a = 0, b = 0;
            if ( 2 != sscanf(str.c_str(), "%f, %f", &a, &b) )
                throw InvalidParameter("expected range ([REAL, REAL]) in Bead's radius");
            rad = RNG.real_uniform(a, b);
        }
        else
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

    Bead * obj = new Bead(pp, Vector(0,0,0), rad);
    
    std::string str;
#if NEW_SOLID_CLAMP
    // clamp position set with 'new'
    if ( opt.set(str, "clamp") )
    {
        if ( str == "position" )
            opt.set(obj->clamp_place, "position");
        else
            obj->clamp_place = Cytosim::readPosition(str, nullptr); //obj->prop->confine_space);
    }
    opt.set(obj->clamp_stiff, "clamp", 1);
    if ( obj->clamp_stiff < 0 )
        throw InvalidParameter("clamp[0] (stiffness) should be >= 0");
#endif
    
    // create list with one object:
    ObjectList res(obj);

    // attach anchored Singles:
    while ( opt.set(str, var, inx++) )
        simul_.singles.makeWrists(res, obj, 0, 1, str);
    return res;
}


void BeadSet::writeSet(Outputter& out) const
{
    if ( size() > 0 )
    {
        out.write("\n#section "+title());
        writePool(out, pool_);
    }
}


void BeadSet::defrostMore()
{
    Object * i;
    while (( i = ice_.front() ))
    {
        ice_.pop_front();
        //std::clog << "delete " << i->reference() << "\n";
        inventory_.unassign(i);
        i->objset(nullptr);
        simul_.singles.deleteWrists(i);
        delete(i);
    }
}


void BeadSet::remove(Object * obj)
{
    ObjectSet::remove(obj);
    simul_.singles.deleteWrists(obj);
}


void BeadSet::foldPositions(Modulo const* m) const
{
    for ( Bead * o=first(); o; o=o->next() )
        o->foldPosition(m);
}
