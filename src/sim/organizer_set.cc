// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "organizer_set.h"
#include "organizer.h"
#include "mecapoint.h"
#include "glossary.h"
#include "nucleus.h"
#include "bundle.h"
#include "aster.h"
#include "fake.h"
#include "solid.h"
#include "simul.h"

//------------------------------------------------------------------------------

void OrganizerSet::steps()
{
    Organizer * obj = first();
    while ( obj )
    {
        Organizer * nxt = obj->next();
        obj->step();
        obj = nxt;
    }
    if ( size() > 1 ) shuffle();
}

//------------------------------------------------------------------------------

Property* OrganizerSet::newProperty(const std::string& cat, const std::string& nom, Glossary&) const
{
    if ( cat == "aster" )   return new AsterProp(nom);
    if ( cat == "bundle" )  return new BundleProp(nom);
    if ( cat == "nucleus" ) return new NucleusProp(nom);
    if ( cat == "fake" )    return new FakeProp(nom);
    return nullptr;
}


Object * OrganizerSet::newObject(const ObjectTag tag, PropertyID pid)
{
    if ( tag == Organizer::ASTER_TAG )
    {
        AsterProp * p = simul_.findProperty<AsterProp>("aster", pid);
        return new Aster(p);
    }
    
    if ( tag == Organizer::BUNDLE_TAG )
    {
        BundleProp * p = simul_.findProperty<BundleProp>("bundle", pid);
        return new Bundle(p);
    }
    
    if ( tag == Organizer::NUCLEUS_TAG )
    {
        NucleusProp * p = simul_.findProperty<NucleusProp>("nucleus", pid);
        return new Nucleus(p);
    }
    
    if ( tag == Organizer::FAKE_TAG )
    {
        FakeProp * p = simul_.findProperty<FakeProp>("fake", pid);
        return new Fake(p);
    }
    
    throw InvalidIO("Warning: unknown Organizer tag `"+std::to_string(tag)+"'");
    return nullptr;
}


ObjectList OrganizerSet::newObjects(Property const* p, Glossary& opt)
{
    Organizer * obj = nullptr;
    
    if ( p->category() == "aster" )
        obj = new Aster(static_cast<AsterProp const*>(p));
    else if ( p->category() == "bundle" )
        obj = new Bundle(static_cast<BundleProp const*>(p));
    else if ( p->category() == "nucleus" )
        obj = new Nucleus(static_cast<NucleusProp const*>(p));
    else if ( p->category() == "fake" )
        obj = new Fake(static_cast<FakeProp const*>(p));

    ObjectList res;
    if ( obj )
    {
        res = obj->build(opt, simul_);
    }
    return res;
}


void OrganizerSet::writeSet(Outputter& out) const
{
    if ( size() > 0 )
    {
        out.write("\n#section "+title());
        writePool(out, pool_);
    }
}

//------------------------------------------------------------------------------

ObjectID OrganizerSet::findOrganizerID(const Mecable * m) const
{
    ObjectID res = 0;
    for ( Organizer const* o=first(); o; o=o->next() )
        if ( o->check(m) )
            res = std::max(res, o->identity());

    return res;
}


Aster * OrganizerSet::pickAster(std::string s) const
{
    return Aster::toAster(pickObject("aster", s));
}

//------------------------------------------------------------------------------

void OrganizerSet::report(std::ostream& os) const
{
    if ( size() > 0 )
    {
        unsigned total = 0;
        os << '\n' << title();
        for ( Property const* i : simul_.properties.find_all("aster", "bundle") )
        {
            size_t cnt = count(i);
            os << '\n' << std::setw(10) << cnt << " " << i->name();
            ++total;
        }
        for ( Property const* i : simul_.properties.find_all("nucleus", "fake") )
        {
            size_t cnt = count(i);
            os << '\n' << std::setw(10) << cnt << " " << i->name();
            ++total;
        }
        if ( total > 1 )
            os << '\n' << std::setw(10) << size() << " total";
    }
}


