// Cytosim was created by Francois Nedelec.
// Copyright Cambridge University, 2019

#include "event_set.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"
#include "tubule_prop.h"
#include "tubule.h"


Tubule * TubuleSet::first() const
{
    return static_cast<Tubule*>(pool_.front());
}


Tubule * TubuleSet::findID(ObjectID n) const
{
    return static_cast<Tubule*>(inventory_.get(n));
}


void TubuleSet::step()
{
    for ( Tubule * e=first(); e; e=e->next() )
        e->step();
}


Property* TubuleSet::newProperty(const std::string& cat, const std::string& nom, Glossary&) const
{
    if ( cat == "tubule" )
        return new TubuleProp(nom);
    return nullptr;
}


Object * TubuleSet::newObject(const ObjectTag tag, PropertyID pid)
{
    if ( tag == Tubule::TAG )
    {
        TubuleProp * p = simul_.findProperty<TubuleProp>("tubule", pid);
        return new Tubule(p);
    }
    throw InvalidIO("Warning: unknown Tubule tag `"+std::to_string(tag)+"'");
    return nullptr;
}


/**
 @defgroup NewTubule How to create a Tubule
 @ingroup NewObject

 Specify a new Tubule:
 
     new tubule NAME
     {
     }
 */
ObjectList TubuleSet::newObjects(const Property* p, Glossary& opt)
{
    TubuleProp const* pp = static_cast<TubuleProp const*>(p);
    Tubule * obj = new Tubule(pp);
    return obj->build(pp->radius, opt, simul_);
}


void TubuleSet::writeSet(Outputter& out) const
{
    if ( size() > 0 )
    {
        out.write("\n#section "+title());
        writeObjects(out, pool_);
    }
}
