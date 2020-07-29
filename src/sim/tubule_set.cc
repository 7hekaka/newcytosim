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
    return static_cast<Tubule*>(nodes.front());
}


Tubule * TubuleSet::findID(ObjectID n) const
{
    return static_cast<Tubule*>(inventory.get(n));
}


void TubuleSet::step()
{
    for ( Tubule * e=first(); e; e=e->next() )
        e->step(simul);
}


Property* TubuleSet::newProperty(const std::string& cat, const std::string& nom, Glossary&) const
{
    if ( cat == "tubule" )
        return new TubuleProp(nom);
    return nullptr;
}


Object * TubuleSet::newObject(const ObjectTag tag, size_t num)
{
    if ( tag == Tubule::TAG )
    {
        TubuleProp * p = simul.findProperty<TubuleProp>("tubule", num);
        return new Tubule(p);
    }
    std::cerr << "Warning: unknown Tubule tag `"+std::string(1,tag)+"' requested\n";
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
ObjectList TubuleSet::newObjects(const std::string& name, Glossary& opt)
{
    ObjectList res;
    TubuleProp * p = simul.findProperty<TubuleProp>("tubule", name);
    Tubule * o = new Tubule(p);
    res.push_back(o);
    res.append(o->build(opt, simul));
    return res;
}


void TubuleSet::write(Outputter& out) const
{
    if ( size() > 0 )
    {
        out.put_line("\n#section "+title(), out.binary());
        writeNodes(out, nodes);
    }
}
