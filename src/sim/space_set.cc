// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University
#include "space_set.h"
#include "space_prop.h"
#include "space_dynamic_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"
#include "space.h"
#include "modulo.h"

//---------------------------- GLOBAL VARIABLES --------------------------------

/**
 set current Space to `spc`. (spc==NULL is a valid argument).
 */
void SpaceSet::setMaster(Space const* spc)
{
    if ( spc != master_ )
    {
        master_ = spc;
        
#if ( 0 )
        if ( spc )
            std::clog << "Space master <--- " << spc->prop->name() << '\n';
        else
            std::clog << "Space master <--- NULL" << '\n';
#endif
    }
    
#if ENABLE_PERIODIC_BOUNDARIES
    modulo = nullptr;

    if ( master_ )
        modulo = master_->getModulo();
#endif
}


real SpaceSet::maxExtension() const
{
    real rad = 0;
    for ( Space const* s = first(); s; s=s->next() )
        rad = std::max(rad, s->maxExtension());
    return rad;
}

//------------------------------------------------------------------------------

Property * SpaceSet::newProperty(const std::string& cat, const std::string& nom, Glossary& opt) const
{
    std::string s;
    if ( cat == "space" )
    {
        if ( opt.peek(s, "shape") )
            ;
#if BACKWARD_COMPATIBILITY < 50
        // "geometry" was used before 2018
        else if ( opt.peek(s, "geometry") )
        {
            std::stringstream iss(s);
            iss >> s;
        }
#endif
        if ( s=="lid" )              return new SpaceDynamicProp(nom);
        if ( s=="dynamic_disc" )     return new SpaceDynamicProp(nom);
        if ( s=="dynamic_sphere" )   return new SpaceDynamicProp(nom);
        if ( s=="dynamic_ellipse")   return new SpaceDynamicProp(nom);
#if BACKWARD_COMPATIBILITY < 50
        if ( s=="contractile" )      return new SpaceDynamicProp(nom);
#endif
        return new SpaceProp(nom);
    }
    return nullptr;
}


void SpaceSet::steps()
{
    Space * obj = first();
    while ( obj )
    {
        Space * nxt = obj->next();
        obj->step();
        // delete object that have been flagged
        if ( ! obj->prop ) eraseObject(obj);
        obj = nxt;
    }
    if ( size() > 1 ) shuffle();
}


void SpaceSet::erase()
{
    ObjectSet::erase();
    
    // simul has lost its current Space:
    setMaster(nullptr);
}

/**
 This will change the Simul current Space if it was not set
*/
void SpaceSet::link(Object * obj)
{
    assert_true(obj->tag() == Space::TAG);
    //std::clog << "SpaceSet::add " << obj << '\n';
    ObjectSet::link(obj);
    
    Space const* m = master();
    if ( !m || obj->identity() < m->identity() )
        setMaster(static_cast<Space*>(obj));
}

/**
 If the Simulation current Space is deleted,
 the 'oldest' remaining Space is chosen to replace it.
 */
void SpaceSet::unlink(Object * obj)
{
    //std::clog << "SpaceSet::remove " << obj << '\n';
    ObjectSet::unlink(obj);

    if ( obj == master() )
    {
        /*
         if the current space was deleted, use the oldest Space available
         */
        Space * spc = first();
        
        for ( Space * s=spc; s; s=s->next() )
            if ( s->identity() < spc->identity() )
                spc = s;
        
        setMaster(spc);
    }
}

//------------------------------------------------------------------------------

Object * SpaceSet::newObject(const ObjectTag tag, PropertyID pid)
{
    if ( tag == Space::TAG )
    {
        SpaceProp * p = simul_.findProperty<SpaceProp>("space", pid);
        Space * s = p->newSpace();
        return s;
    }
    throw InvalidIO("Warning: unknown Space tag `"+std::to_string(tag)+"'");
    return nullptr;
}

/**
 The dimensions of a Space can be specified when it is created
 
     new cell
     {
        length = 3, 4
     }
 
 */
ObjectList SpaceSet::newObjects(Property const* p, Glossary& opt)
{
    SpaceProp const* pp = static_cast<SpaceProp const*>(p);
    Space * obj = pp->newSpace(opt);

    if ( !obj )
    {
        throw InvalidParameter("unknown space:shape `"+pp->shape+"'");
        //std::cerr << "Warning: substituting unbounded Space for unknown `"+p->shape+"'\n";
        //obj = new Space(p);
    }
    
    return ObjectList(obj);
}


void SpaceSet::writeSet(Outputter& out) const
{
    if ( size() > 0 )
    {
        out.write("\n#section "+title());
        writePool(out, pool_);
    }
}


void SpaceSet::report(std::ostream& os) const
{
    if ( size() > 0 )
    {
        os << '\n' << title();
        PropertyList plist = simul_.properties.find_all(title());
        for ( Property const* i : plist )
        {
            SpaceProp const* p = static_cast<SpaceProp const*>(i);
            size_t cnt = count(p);
            os << '\n' << std::setw(10) << cnt << ' ' << p->name();
            os << " ( " << p->shape << " )";
        }
        if ( plist.size() > 1 )
            os << '\n' << std::setw(10) << size() << " total";
    }
}
