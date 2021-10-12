// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "cymdef.h"
#include "simul.h"
#include "messages.h"
#include "glossary.h"
#include "iowrapper.h"
#include "exceptions.h"
#include "hand_prop.h"
#include "simul_prop.h"
#include "backtrace.h"
#include "modulo.h"
#include "tictoc.h"

#include "tubule.h"
#include "fiber.h"
#include "field.h"
#include "event.h"
#include "parser.h"

extern Modulo const* modulo;

#include "simul_step.cc"
#include "simul_file.cc"
#include "simul_custom.cc"
#include "simul_report.cc"
#include "simul_solve.cc"
#include "simul_solvef.cc"

#include <csignal>

//---------------------------  global variables/functions ---------------------

void out_of_memory_handler()
{
    (void) write(STDERR_FILENO, "\n* * * * *\n", 11);
    (void) write(STDERR_FILENO, "Cytosim: memory allocation failed", 33);
    (void) write(STDERR_FILENO, "\n* * * * *\n", 11);
    print_backtrace();
    _exit(1);
}

void termination_handler()
{
    (void) write(STDERR_FILENO, "\n* * * * *\n", 11);
    (void) write(STDERR_FILENO, "Cytosim: uncaught exception", 27);
    (void) write(STDERR_FILENO, "\n* * * * *\n", 11);
    print_backtrace();
    abort();
}

void signal_handler(int sig)
{
    (void) write(STDERR_FILENO, "\n* * * * *\n", 11);
    psignal((unsigned)sig, "Cytosim");
    (void) write(STDERR_FILENO, "* * * * *\n", 10);
    print_backtrace();
    _exit(sig);
}

//------------------------------------------------------------------------------
#pragma mark -

Simul::Simul()
: prop(nullptr), spaces(*this), fields(*this),
fibers(*this), spheres(*this), beads(*this), solids(*this),
singles(*this), couples(*this), organizers(*this), tubules(*this), events(*this)
{
    pMeca1D = nullptr;
    parser_ = nullptr;
    primed_ = 0;
#if POOL_HAND_ATTACHMENT
    dontAttach = 0;
#endif
    autoPrecond = 0;
    autoCounter = 0;
    for ( size_t u = 0; u < 8; ++u )
    {
        autoCPU[u] = 0;
        autoCNT[u] = 0;
    }
    
    prop = new SimulProp("undefined");
}

Simul::~Simul()
{
    erase_all(1);
    delete(pMeca1D);
    delete(prop);
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 This will initialize the simulation by registering callbacks.
 You should still call Simul::prepare() before calling step()
 */
void Simul::initialize(Glossary & glos)
{
    // Register a function to be called if operator new fails:
    std::set_new_handler(out_of_memory_handler);
    
    // Register a function to be called upon abortion:
    std::set_terminate(termination_handler);
    
    // Register a function to be called for Floating point exceptions:
    if ( signal(SIGFPE, signal_handler) == SIG_ERR )
        std::cerr << "Could not register SIGFPE handler\n";
    
    // read parameters
    prop->read(glos);
}

//------------------------------------------------------------------------------

double Simul::time() const
{
    return prop->time;
}

real Simul::time_step() const
{
    return prop->time_step;
}

void Simul::erase_all(bool erase_properties)
{
    //std::cerr << "Simul::erase()\n";
    relax();
    tubules.erase();
    organizers.erase();
    fibers.erase();
    singles.erase();
    couples.erase();
    spheres.erase();
    beads.erase();
    solids.erase();
    fields.erase();
    spaces.erase();
    events.erase();
 
    prop->time = 0;
    modulo = nullptr;
    
    if ( erase_properties )
        properties.erase();
}


size_t Simul::nbObjects() const
{
    return  (  organizers.size()
             + singles.size()
             + couples.size()
             + fibers.size()
             + beads.size()
             + solids.size()
             + spheres.size()
             + spaces.size()
             + fields.size() );
}


void Simul::foldPositions() const
{
    if ( modulo )
    {
        fibers.foldPositions(modulo);
        beads.foldPositions(modulo);
        solids.foldPositions(modulo);
        spheres.foldPositions(modulo);
        singles.foldPositions(modulo);
        couples.foldPositions(modulo);
    }
}


void Simul::evaluate(std::string const& code)
{
    //Using the parser which has been given at the start of the config file
    assert_true(parser_);
    parser_->evaluate(code);
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 Convert Object pointer to Mecable* if possible
 */
Mecable* Simul::toMecable(Object * obj)
{
    if ( obj )
    switch( obj->tag() )
    {
        case  Fiber::TAG:  return static_cast<Mecable*>(obj);
        case   Bead::TAG:  return static_cast<Mecable*>(obj);
        case  Solid::TAG:  return static_cast<Mecable*>(obj);
        case Sphere::TAG:  return static_cast<Mecable*>(obj);
    }
    return nullptr;
}

/**
 Find an object from a Human-friendly representation, such as
 fiber1
 single1
 */
Mecable * Simul::findMecable(const std::string& arg) const
{
    Object  * obj = fibers.findObject(arg, fibers.title());
    if (!obj) obj = solids.findObject(arg, solids.title());
    if (!obj) obj = spheres.findObject(arg, spheres.title());
    if (!obj) obj = beads.findObject(arg, beads.title());
    return static_cast<Mecable*>(obj);
}


void Simul::add(Object * w)
{
    //std::clog << " Simul::add(" << w->reference() << ")" << '\n';
    assert_true(w);
    ObjectSet * set = findSetT(w->tag());
    set->add(w);
}


void Simul::add(ObjectList const& objs)
{
    //std::clog << " Simul::add("<< objs.size() <<" objects):" << '\n';
    for ( Object * obj : objs )
        if ( obj )
            add(obj);
}


void Simul::remove(Object * w)
{
    assert_true( w->objset() );
    w->objset()->remove(w);
}


void Simul::remove(ObjectList const& objs)
{
    //std::clog << " Simul::remove("<< objs.size() <<" objects):" << '\n';
    for ( Object * obj : objs )
        if ( obj )
            remove(obj);
}


void Simul::erase(Object * w)
{
    //std::clog << "Simul::erase " << w->reference() << '\n';
    remove(w);
    delete(w);
}


void Simul::erase(ObjectList const& objs)
{
    //std::clog << " Simul::erase(" << objs.size() << " objects):\n";
    for ( Object * obj : objs )
        if ( obj )
        {
            //std::clog << " Simul::erase(" << obj << ")\n";
            remove(obj);
            delete(obj);
        }
}


void Simul::erase(ObjectList& objs, size_t cnt)
{
    if ( cnt == 1 )
    {
        erase(objs.pick_one());
    }
    else
    {
        if ( cnt < objs.size() )
        {
            // limit the list to a random subset
            objs.shuffle();
            objs.truncate(cnt);
        }
        
        //std::clog << "simul:deleting " << objs.size() << " " << set->title() << '\n';
        erase(objs);
    }
}


void Simul::mark(ObjectList const& objs, ObjectMark mrk)
{
    //std::clog << " Simul::erase("<< objs.size() <<" objects):" << '\n';
    for ( Object * i : objs )
        i->mark(mrk);
}

//------------------------------------------------------------------------------
#pragma mark -

Space const* Simul::findSpace(std::string const& str) const
{
    if ( str == "first" )
        return static_cast<Space*>(spaces.inventory_.first());

    if ( str == "last" )
        return static_cast<Space*>(spaces.inventory_.last());
    
    Property * sp = properties.find("space", str);
    
    if ( sp )
        return spaces.findObject(sp);
    
    return nullptr;
}

/**
 This is used primarily to parse the configuration file,
 argument is the full class name
 */
ObjectSet * Simul::findSet(const std::string& cat)
{
    //std::clog << "findSet("<<kind<<")\n";
    if ( cat == spaces.title() )     return &spaces;
    if ( cat == fields.title() )     return &fields;
    if ( cat == fibers.title() )     return &fibers;
    if ( cat == beads.title() )      return &beads;
    if ( cat == solids.title() )     return &solids;
    if ( cat == spheres.title() )    return &spheres;
    if ( cat == singles.title() )    return &singles;
    if ( cat == couples.title() )    return &couples;
    if ( cat == tubules.title() )    return &tubules;
    if ( cat == organizers.title() ) return &organizers;
    if ( cat == "aster" )            return &organizers;
    if ( cat == "bundle" )           return &organizers;
    if ( cat == "nucleus" )          return &organizers;
    if ( cat == "fake" )             return &organizers;
    if ( cat == events.title() )     return &events;
    return nullptr;
}


/**
 This is used primarily to read the binary trajectory file,
 using a single character to refer to each class in Cytosim
 */
ObjectSet * Simul::findSetT(const ObjectTag tag)
{
    switch( tag )
    {
        case        Couple::TAG: return &couples;
        case        Single::TAG: return &singles;
        case  Single::TAG_WRIST: return &singles;
        case         Fiber::TAG: return &fibers;
        case     Fiber::TAG_ALT: return &fibers;
        case Fiber::TAG_DYNAMIC: return &fibers;
        case Fiber::TAG_LATTICE: return &fibers;
        case Fiber::TAG_FIBMESH: return &fibers;
#if BACKWARD_COMPATIBILITY < 57
        case 'l': return &fibers; // TAG_LATTICE before 23/06/2021
        case 'L': return &fibers; // TAG_FIBMESH before 23/06/2021
#endif
        case          Bead::TAG: return &beads;
        case         Solid::TAG: return &solids;
        case        Sphere::TAG: return &spheres;
        case         Field::TAG: return &fields;
        case         Space::TAG: return &spaces;
        case        Tubule::TAG: return &tubules;
        case         Event::TAG: return &events;
        case Organizer::TAG_NUCLEUS: return &organizers;
        case Organizer::TAG_BUNDLE:  return &organizers;
        case Organizer::TAG_FAKE:    return &organizers;
        case Organizer::TAG_ASTER:   return &organizers;
        case   Object::NULL_TAG: return nullptr;
    }
    return nullptr;
}

//------------------------------------------------------------------------------
#pragma mark -

/* There can only be one SimulProp and it is already created */
void Simul::rename(std::string const& arg)
{
    if ( prop->name() == "undefined" )
    {
        prop->rename(arg);
        //std::clog << "Simul is named `" << arg << "'\n";
    }
    else if ( prop->name() != arg )
        throw InvalidSyntax("only one `simul' can be defined");
}


bool Simul::isCategory(const std::string& name) const
{
    if ( name == "hand" )
        return true;
    if ( name == "simul" )
        return true;

    return const_cast<Simul*>(this)->findSet(name);
}


Property* Simul::findProperty(const std::string& cat, const std::string& nom) const
{
    if ( cat == "simul" && nom == prop->name() )
        return prop;

    if ( cat.empty() || nom.empty() )
        throw InvalidSyntax("unexpected syntax");

    return properties.find(cat, nom);
}


Property* Simul::findProperty(const std::string& nom) const
{
    if ( nom == prop->name() )
        return prop;

    if ( nom.empty() )
        throw InvalidSyntax("unexpected syntax");

    return properties.find(nom);
}


PropertyList Simul::findAllProperties(const std::string& cat) const
{
    if ( cat == "simul" )
    {
        PropertyList list;
        list.push_back(prop);
        return list;
    }
    if ( cat.empty() )
        throw InvalidSyntax("unexpected syntax");

    return properties.find_all(cat);
}


/**
 @defgroup ObjectGroup List of objects
 
 The command `set simul` will define the global parameters.
 The `simul` is automatically created, and you cannot use 'new simul'.

 Objects       | Base class    | Parameters    |
 --------------|---------------|----------------
 `simul`       |  Simul        | @ref SimulPar  
 
 
 These objects cannot move:
 
 Class Name    | Base class    | Parameters       | Specialization   |
 --------------|---------------|------------------|-------------------
 `space`       |  Space        | @ref SpacePar    | @ref SpaceGroup
 `field`       |  Field        | @ref FieldPar    | -                 
 `event`       |  Event        | @ref EventPar    | -
 
 
 `Mecables` can move or deform, and come in 4 basic forms:
 
 Class Name    | Base class    | Parameters       | Specialization   |
 --------------|---------------|------------------|-------------------
 `fiber`       |  Fiber        | @ref FiberPar    | @ref FiberGroup
 `bead`        |  Bead         | @ref SolidPar    | -
 `solid`       |  Solid        | @ref SolidPar    | -
 `sphere`      |  Sphere       | @ref SpherePar   | -

 
 A `Hand` is an object that can bind to fiber, but it can only be used
 as a sub-part of `Single` or `Couple`.

 Class Name    | Base class    | Parameters       | Specialization   |
 --------------|---------------|------------------|-------------------
 `hand`        |  Hand         | @ref HandPar     | @ref HandGroup
 
 
 `Single` and `Couple` contain one or two `Hand` respectively:

 Class Name    | Base class    | Parameters       | Specialization   |
 --------------|---------------|------------------|-------------------
 `single`      |  Single       | @ref SinglePar   | @ref SingleGroup
 `couple`      |  Couple       | @ref CouplePar   | @ref CoupleGroup
 
 
 
 The `Organizers` describe composite objects build from multiple Mecables:
 
 Organizers    | Base class    | Parameters      |
 --------------|---------------|------------------
 `aster`       |  Aster        | @ref AsterPar    
 `bundle`      |  Bundle       | @ref BundlePar   
 `nucleus`     |  Nucleus      | @ref NucleusPar  
 `fake`        |  Fake         | @ref FakePar     
 .
 
 */
Property* Simul::newProperty(const std::string& cat, const std::string& nom, Glossary& glos)
{
    if ( cat.empty() || nom.empty() )
        throw InvalidSyntax("unexpected syntax");

    if ( cat == "simul" )
    {
        assert_true(prop);
        rename(nom);
        return prop;
    }
    
    if ( isCategory(nom) )
        throw InvalidSyntax("`"+nom+"' is a reserved keyword");
    
    Property * p = findProperty(nom);
    
    if ( p )
        throw InvalidSyntax("property `"+nom+"' is already defined");
    
    if ( cat == "hand" )
    {
        p = HandProp::newProperty(nom, glos);
        properties.deposit(p);
    }
    else
    {
        ObjectSet * set = findSet(cat);
        
        if ( !set )
            throw InvalidSyntax("unknown class `"+cat+"'");
        
        p = set->newProperty(cat, nom, glos);
        properties.deposit(p);
    }
    
    return p;
}

