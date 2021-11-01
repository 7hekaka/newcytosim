// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
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
 This is a global variable that is initialized in Simul
 It is used to implement periodic boundary conditions
 */
Modulo const* modulo = nullptr;


/// static variable of SpaceSet:
Space const* SpaceSet::master_ = nullptr;

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
            std::clog << "Space master: " << spc->prop->name() << '\n';
        else
            std::clog << "Space master: NULL" << '\n';
#endif
    }
    
    modulo = nullptr;

    if ( master_ )
        modulo = master_->getModulo();
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


void SpaceSet::step()
{
    for ( Space * sp = first(); sp; sp=sp->next() )
        sp->step();
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
ObjectList SpaceSet::newObjects(const std::string& name, Glossary& opt)
{
    SpaceProp * p = simul_.findProperty<SpaceProp>("space", name);
    Space * obj = p->newSpace(opt);

    if ( !obj )
    {
        throw InvalidParameter("unknown space:shape `"+p->shape+"'");
        //std::cerr << "Warning: substituting unbounded Space for unknown `"+p->shape+"'\n";
        //obj = new Space(p);
    }
    
    ObjectList res(1, 4);
    res.push_back(obj);
    return res;
}


void SpaceSet::write(Outputter& out) const
{
    if ( size() > 0 )
    {
        out.writeLine("\n#section "+title());
        writeObjects(out, pool_);
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
            size_t cnt = count(match_property, p);
            os << '\n' << std::setw(10) << cnt << ' ' << p->name();
            os << " ( " << p->shape << " )";
        }
        if ( plist.size() > 1 )
            os << '\n' << std::setw(10) << size() << " total";
    }

}
