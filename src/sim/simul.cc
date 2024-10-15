// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University

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
#include "random_seed.h"

#include "fiber.h"
#include "field.h"
#include "event.h"
#include "parser.h"

const char Simul::TRAJECTORY[] = "objects.cmo";

#include "simul_step.cc"
#include "simul_file.cc"
#include "simul_custom.cc"
#include "simul_report.cc"
#include "simul_solve.cc"

//------------------------------------------------------------------------------
#pragma mark -

Simul::Simul()
: prop(""), spaces(*this), fields(*this), fibers(*this),
spheres(*this), beads(*this), solids(*this), singles(*this),
couples(*this), organizers(*this), events(*this)
{
    pMeca1D = nullptr;
    parser_ = nullptr;
    primed_ = 0;
#if POOL_UNATTACHED > 1
    doAttachCounter = 0;
#elif POOL_UNATTACHED < 1
    ABORT_NOW(" POOL_UNATTACHED must be >= 1");
#endif
    autoPrecond = 0;
    autoCounter = 0;
    for ( size_t u = 0; u < 8; ++u )
    {
        autoCPU[u] = 0;
        autoCNT[u] = 0;
    }
    fresh_ = 1;
}

Simul::~Simul()
{
    eraseObjects();
    eraseProperties();
    delete(pMeca1D);
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 This will initialize the Random number generator if needed
 */
void Simul::initCytosim()
{
    if ( !RNG.seeded() )
    {
        if ( prop.random_seed )
            RNG.seed(prop.random_seed);
        else {
            prop.random_seed = PCG32::get_random_seed();
            RNG.seed(prop.random_seed);
        }
    }
}

//------------------------------------------------------------------------------

void Simul::eraseObjects()
{
    //fprintf(stderr, "Simul:%p:eraseObjects()\n", this);
    organizers.erase();
    fibers.erase();
    spheres.erase();
    beads.erase();
    solids.erase();
    fields.erase();
    spaces.erase();
    events.erase();
    singles.erase();
    couples.erase();
    
    prop.time = 0;
    prop.end_time = INFINITY;
#if ENABLE_PERIODIC_BOUNDARIES
    modulo = nullptr;
#endif
    primed_ = 0;
    fresh_ = 1;
}

void Simul::eraseProperties()
{
    properties.erase();
    prop.clear();
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


/// Using the parser which is set by Simul::parser()
void Simul::perform(std::string const& code)
{
    if ( parser_ )
        parser_->evaluate(code);
    else
        throw InvalidParameter("no parser specified!");
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
 Find an object from a user-defined specification,
 such as `fiber1` or `single1`
 */
Mecable * Simul::pickMecable(const std::string& arg) const
{
    Object  * obj = fibers.pickObject("fiber", arg);
    if (!obj) obj = solids.pickObject("solid", arg);
    if (!obj) obj = spheres.pickObject("sphere", arg);
    if (!obj) obj = beads.pickObject("bead", arg);
    return static_cast<Mecable*>(obj);
}


Object * Simul::pickMovable(const std::string& arg) const
{
    Object *  obj = fibers.pickObject("fiber", arg);
    if (!obj) obj = solids.pickObject("solid", arg);
    if (!obj) obj = spheres.pickObject("sphere", arg);
    if (!obj) obj = beads.pickObject("bead", arg);
    if (!obj) obj = couples.pickObject("couple", arg);
    if (!obj) obj = singles.pickObject("single", arg);
    return obj;
}


void Simul::add(Object * w)
{
    assert_true(w);
    ObjectSet * set = findSetT(w->tag());
    set->add(w);
    //std::clog << " Simul::add(" << w->reference() << ")" << '\n';
}


void Simul::add(ObjectList const& objs)
{
    //std::clog << " Simul::add("<< objs.size() <<" objects):" << '\n';
    for ( Object * obj : objs )
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
        remove(obj);
}


void Simul::eraseObject(Object * w)
{
    //std::clog << "Simul::erase " << w->reference() << '\n';
    remove(w);
    delete(w);
}


void Simul::eraseObjects(ObjectList const& objs)
{
    //std::clog << " Simul::eraseObjects(" << objs.size() << " objects):\n";
    for ( Object * obj : objs )
    {
        //std::clog << " Simul::erase(" << obj << ")\n";
        remove(obj);
        delete(obj);
    }
}

//------------------------------------------------------------------------------
#pragma mark -


ObjectFlag Simul::setUniqueFlags() const
{
    size_t f = 0;
    for ( Fiber * F = fibers.firstID(); F; F = fibers.nextID(F) )
        F->flag(ObjectFlag(++f));
    for ( Solid * S = solids.firstID(); S; S = solids.nextID(S) )
        S->flag(ObjectFlag(++f));
    for ( Bead  * B = beads.firstID(); B; B = beads.nextID(B) )
        B->flag(ObjectFlag(++f));
    for ( Sphere* O = spheres.firstID(); O; O = spheres.nextID(O) )
        O->flag(ObjectFlag(++f));
    if ( f != (ObjectFlag)f )
        throw InvalidParameter("ObjectFlag overflow in setUniqueFlags()");
    return ObjectFlag(f);
}


void Simul::setFlags(ObjectFlag f) const
{
    for ( Fiber * F=fibers.first(); F; F=F->next() )
        F->flag(f);
    for ( Solid * S=solids.first(); S; S=S->next() )
        S->flag(f);
    for ( Bead  * B=beads.first(); B; B=B->next() )
        B->flag(f);
    for ( Sphere* O=spheres.first(); O; O=O->next() )
        O->flag(f);
}


void Simul::changeFlags(ObjectFlag f, ObjectFlag g) const
{
    for ( Fiber * F=fibers.first(); F; F=F->next() )
        if ( F->flag() == f ) F->flag(g);
    for ( Solid * S=solids.first(); S; S=S->next() )
        if ( S->flag() == f ) S->flag(g);
    for ( Bead  * B=beads.first(); B; B=B->next() )
        if ( B->flag() == f ) B->flag(g);
    for ( Sphere* O=spheres.first(); O; O=O->next() )
        if ( O->flag() == f ) O->flag(g);
}

//------------------------------------------------------------------------------
#pragma mark -

/** This should be equivalent to ObjectSet::findObject() */
Space const* Simul::findSpace(std::string spec) const
{
    if ( spec == "first" )
        return static_cast<Space const*>(spaces.inventory_.first());

    if ( spec == "last" )
        return static_cast<Space const*>(spaces.inventory_.last());
    
    if ( spec == "master" )
        return static_cast<Space const*>(spaces.master());

    // check if a Space name was specified:
    Property * sp = properties.find("space", spec);
    
    if ( sp )
        return pickSpace(sp);

    // check if `name + id` was specified:
    long num = 0;
    if ( Tokenizer::split_polysymbol(spec, num) )
    {
        Object* obj = spaces.findObject(spaces.title(), spec, num);
        return static_cast<Space const*>(obj);
    }
    
    return nullptr;
}

Field * Simul::pickField(const Property * p) const
{
    return static_cast<Field*>(fields.pickObject(p));
}

/**
 This is used primarily to parse the configuration file,
 argument is the full class name
 */
ObjectSet * Simul::findSet(const std::string& cat)
{
    //std::clog << "findSet("<<cat<<")\n";
    if ( cat == spaces.title() )     return &spaces;
    if ( cat == fields.title() )     return &fields;
    if ( cat == fibers.title() )     return &fibers;
    if ( cat == beads.title() )      return &beads;
    if ( cat == solids.title() )     return &solids;
    if ( cat == spheres.title() )    return &spheres;
    if ( cat == singles.title() )    return &singles;
    if ( cat == couples.title() )    return &couples;
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
    assert_true(islower(tag));
    switch( tag )
    {
        case        Couple::TAG: return &couples;
        case    Couple::DUO_TAG: return &couples;
        case        Single::TAG: return &singles;
        case  Single::WRIST_TAG: return &singles;
        case         Fiber::TAG: return &fibers;
        case Fiber::COMPACT_TAG: return &fibers;
#if BACKWARD_COMPATIBILITY < 57
        case 'l': return &fibers; // LATTICE_TAG before 23/06/2021
        case 'L': return &fibers; // FIBMESH_TAG before 23/06/2021
#endif
        case          Bead::TAG: return &beads;
        case         Solid::TAG: return &solids;
        case        Sphere::TAG: return &spheres;
        case         Field::TAG: return &fields;
        case         Space::TAG: return &spaces;
        case         Event::TAG: return &events;
        case Organizer::NUCLEUS_TAG: return &organizers;
        case  Organizer::BUNDLE_TAG: return &organizers;
        case    Organizer::FAKE_TAG: return &organizers;
        case   Organizer::ASTER_TAG: return &organizers;
#if BACKWARD_COMPATIBILITY < 60
        case 'v': return &fibers; // NULL_TAG before 2/12/2023
#endif
        case Object::NULL_TAG: return nullptr;
    }
    return nullptr;
}

//------------------------------------------------------------------------------
#pragma mark -

void Simul::rename(std::string const& arg)
{
    if ( prop.name() != arg )
    {
        if ( arg == "display" || isCategory(arg) )
            throw InvalidSyntax("`"+arg+"' is a reserved keyword");
        if ( prop.name().size() )
            throw InvalidSyntax("Simul `"+prop.name()+"' already defined");
        prop.rename(arg);
        //std::clog << "Simul renamed `" << arg << "'\n";
    }
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
    if ( cat == "simul" && nom == prop.name() )
        return &prop;

    if ( cat.empty() || nom.empty() )
        throw InvalidSyntax("findProperty(void, void)");

    return properties.find(cat, nom);
}


Property* Simul::findProperty(const std::string& nom) const
{
    if ( nom == "simul" || nom == prop.name() )
        return &prop;

    if ( nom.empty() )
        throw InvalidSyntax("findProperty(void)");

    return properties.find(nom);
}


PropertyList Simul::findAllProperties(const std::string& cat) const
{
    if ( cat == "simul" )
    {
        PropertyList list;
        list.push_back(&prop);
        return list;
    }
    if ( cat.empty() )
        throw InvalidSyntax("findAllProperty(void)");

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
Property* Simul::makeProperty(const std::string& cat, const std::string& nom, Glossary& glos)
{
    if ( cat.empty() || nom.empty() )
        throw InvalidSyntax("unexpected syntax");
    
    /* We do not permit using a class name to name a property,
     as this is tempting, but would create confusion in the config file */
    if ( nom == "display" || isCategory(nom) )
        throw InvalidSyntax("`"+nom+"' is a reserved keyword");

    if ( cat == "simul" )
    {
        rename(nom);
        return &prop;
    }
    
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


/**
 read and set parameter for some object in the simulation,
 syntax is `CLASS:PARAMETER=VALUE` or `CLASS.PARAMETER=VALUE`
 Example:
 `simul.viscosity=1`
 `kinesin.stall_force=1` if `kinesin` is a defined class
*/
bool Simul::readParameter(const char* arg) const
{
    char tmp[512] = { 0 };
    strncpy(tmp, arg, sizeof(tmp));
    char * val = tmp;
    char * cls = strsep(&val, ":.");
    char * tok = strsep(&val, "=");
    printf("Simul::readParameter(%s|%s|%s)\n", cls, tok, val);
    Property * P = findProperty(cls);
    if ( P && val )
    {
        Glossary glos(tok, val);
        //P->write_values_diff(std::clog, true);
        P->read(glos);
        //std::cerr << " understood `" << argv[n] << "' :\n";
        //P->write_values_diff(std::clog, true);
        unsigned res = glos.num_reads(tok);
        if ( ! res ) {
            print_yellow(stderr, "Warning: could not read any value from `"+std::string(arg)+"'\n");
        }
        return res;
    }
    return 0;
}
