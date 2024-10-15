// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include <fstream>
#include <cstdio>

#include "cymdef.h"
#include "primitives.h"
#include "interface.h"
#include "stream_func.h"
#include "exceptions.h"
#include "simul_prop.h"
#include "tokenizer.h"
#include "evaluator.h"
#include "messages.h"
#include "glossary.h"
#include "time_date.h"
#include "filepath.h"
#include "simul.h"
#include "event.h"


// Use second definition to trace execution
#define VLOG(ARG) ((void) 0)
//#define VLOG(ARG) std::clog << ARG << '\n';

//------------------------------------------------------------------------------

Interface::Interface(Simul* s)
: sim_(s)
{
}


SimulProp const& Interface::simulProp() const
{
    return sim_->prop;
}

bool Interface::isCategory(std::string const& name) const
{
    return sim_->isCategory(name);
}

void Interface::eraseSimul(bool arg) const
{
    sim_->eraseObjects();
    if ( arg )
        sim_->eraseProperties();
}


ObjectSet* Interface::findClass(std::string const& name, Property*& pp)
{
    ObjectSet * set = nullptr;
    pp = sim_->properties.find(name);
    if ( pp )
        set = sim_->findSet(pp->category());
    else
        set = sim_->findSet(name);
    return set;
}


Object* Interface::findObject(std::string const& name, Property*& pp)
{
    ObjectSet * set = nullptr;
    // search for object name, eg. `microtubule1`:
    long num = 0;
    std::string str = name;
    if ( Tokenizer::split_polysymbol(str, num) )
    {
        pp = sim_->properties.find(str);
        if ( pp )
        {
            set = sim_->findSet(pp->category());
            if ( set )
                return set->findObject(pp, num);
        }
    }
    return nullptr;
}

//------------------------------------------------------------------------------
#pragma mark -


/**
 This creates a new Property
 
 Property::complete() is called after a property is set.
 This ensures that inconsistencies are detected as early as possible.
 
 In addition, we call complete() for all Properties, when the simulation is
 about to start.
 */
Property* Interface::execute_set(std::string const& cat, std::string const& name, Glossary& def)
{
    VLOG("+SET " << cat << " `" << name << "'");
    
    Property* pp = sim_->makeProperty(cat, name, def);
    
    if ( !pp )
        throw InvalidSyntax("failed to create property of class `"+cat+"'");
    
    pp->read(def);
    pp->complete(*sim_);
    
    return pp;
}


void Interface::change_property(Property * pp, Glossary& def)
{
    pp->read(def);
    pp->complete(*sim_);
    
    /*
     Specific code to make 'change space:dimension' work.
     This is needed as dimensions belong to Space, and not to SpaceProp
     */
    if ( pp->category() == "space" )
    {
        // update any Space with this property:
        for ( Space * s = sim_->spaces.first(); s; s=s->next() )
        {
            if ( s->prop == pp )
            {
                s->resize(def);
                // allow Simul to update periodic:
                if ( s == sim_->spaces.master() )
                    sim_->spaces.setMaster(s);
            }
        }
    }
}


// in this form, 'name' designates the property name
Property * Interface::execute_change(std::string const& name, Glossary& def, bool strict)
{
    Property * pp = sim_->findProperty(name);
    
    if ( pp )
    {
        VLOG("-CHANGE " << pp->category() << " `" << name << "'");
        change_property(pp, def);
    }
    else
    {
        if ( strict )
        {
            InvalidParameter e("unknown property `"+name+"'");
            e << sim_->properties.all_names(PREF);
            throw e;
        }
        else
        {
            VLOG("unknown change |" << name << "|");
        }
    }
    return pp;
}


void Interface::execute_change_all(std::string const& cat, Glossary& def)
{
    PropertyList plist = sim_->findAllProperties(cat);
    
    for ( Property * i : plist )
    {
        VLOG("+CHANGE " << i->category() << " `" << i->name() << "'");
        change_property(i, def);
    }
    /*
    if ( plist.size() == 0 )
        throw InvalidSyntax("could not find any property of class `"+cat+"'");
     */
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 Define a placement = ( position, orientation ) from parameters in `opt'
 */
bool Interface::read_placement(Isometry& iso, Glossary& opt)
{
    std::string str;
    Space const* spc = sim_->spaces.master();
    
    // Space specified as second argument to 'position'
    if ( opt.set(str, "position", 1) )
    {
        spc = sim_->findSpace(str);
        if ( ! spc )
            throw InvalidParameter("unknown Space `"+str+"'");
    }
    
    // Position
    Vector vec(0,0,0);
    if ( opt.set_block(str, '[', "position") )
    {
        // can specify another object to copy its position
        Movable * obj = sim_->pickMovable(str);
        if ( obj )
            vec = obj->position();
        else
            throw InvalidParameter("Could not find Object `"+str+"'");
    }
    else if ( opt.set(str, "position") )
        vec = Cytosim::readPosition(str, spc);
    else if ( spc )
        vec = spc->place();
    
    if ( !vec.valid() )
        return false;
    
    iso.mov = vec;
    
    // Rotation applied before the translation
    if ( opt.set_block(str, '[', "direction") )
    {
        // can specify another object to copy its orientation
        Movable * obj = sim_->pickMovable(str);
        if ( obj )
            vec = obj->direction();
        else
            throw InvalidParameter("Could not find Object `"+str+"'");
        iso.rot = Rotation::randomRotationToVector(vec);
    }
    else if ( opt.set_from_least_used_value(str, "direction") )
    {
        vec = Cytosim::readDirection(str, iso.mov, spc);
        iso.rot = Rotation::randomRotationToVector(vec);
    }
    else if ( opt.set(str, "rotation") )
    {
        iso.rot = Cytosim::readRotation(str);
    }
    else if ( opt.set(str, "orientation") )
    {
        size_t sci = 0;
        iso.rot = Cytosim::readOrientation(str, sci, iso.mov, spc);
        size_t t = Cytosim::has_trail(str, sci);
        if ( t )
            throw InvalidSyntax("unexpected `"+str.substr(t)+"' in `orientation = "+str+"'");
    }
    else
        iso.rot = Rotation::randomRotation();
    
    // Another rotation can be specified, to be applied after the translation
    if ( opt.set(str, "rotation", 1) )
    {
        iso.rotate(Cytosim::readRotation(str));
    }
    else if ( opt.set(str, "orientation", 1) )
    {
        size_t sci = 0;
        Rotation rot = Cytosim::readOrientation(str, sci, iso.mov, spc);
        size_t t = Cytosim::has_trail(str, sci);
        if ( t )
            throw InvalidSyntax("unexpected `"+str.substr(t)+"' in `orientation = "+str+"'");
        iso.rotate(rot);
    }
    
    return true;
}


enum PlacementType { PLACE_NOT = 0, PLACE_ANYWHERE, PLACE_INSIDE, PLACE_EDGE,
                     PLACE_OUTSIDE, PLACE_ALL_INSIDE };


/**
 
     new INTEGER CLASS NAME
     {
       position = POSITION
       placement = PLACEMENT, SPACE_NAME, CONDITION
       nb_trials = INTEGER
     }
 
 PLACEMENT can be:
 - `inside` (default), the object is created with a majority of vertices inside the Space
 - `all_inside`: the object is created only if all its vertices are inside
 - `outside`, the object is created only if it is placed outside the Space
 - `surface`, the position is projected on the edge of current Space
 - `anywhere`, the object is placed in a random position with random orientation
 - `off`: translation/rotation are not applied
 .
 
 By default, the specifications are relative to the first Space to be defined,
 but a different space can be specified as second argument of PLACEMENT.
 
 You can set the density of objects with `nb_trials=1`:
 
     new 100 grafted
     {
       position = ( rectangle 10 10 )
       nb_trials = 1
     }
 
 In this way an object will be created only if its randomly chosen position falls
 inside the Space, and the density will thus be exactly what is specified from the
 `position` range (here 100/10*10 = 1 object per squared micrometer).
 */
bool Interface::find_placement(Isometry& iso, Glossary& opt, int placement)
{
    std::string str;
    Space const* spc = sim_->spaces.master();
    if ( opt.set(str, "placement", 1) )
        spc = sim_->findSpace(str);

    // generate a new position:
    bool valid = read_placement(iso, opt);
    
    if ( !valid )
        return 0;
    
    // check any conditions to the position:
    bool has_condition = opt.set(str, "placement", 2);
    if ( has_condition )
    {
        Evaluator evaluator{{"X", iso.mov.x()}, {"Y", iso.mov.y()}, {"Z", iso.mov.z()},
            {"R", iso.mov.norm()}, {"P", RNG.preal()}};
        try {
            if ( 0 == evaluator.eval(str) )
                return 0;
        }
        catch( Exception& e ) {
            e.message(e.message()+" in `"+str+"'");
            throw;
        }
    }
    
    if ( !spc || placement == PLACE_ANYWHERE )
        return 1;
    
    if ( placement == PLACE_EDGE )
    {
        iso.mov = spc->project(iso.mov);
        return 1;
    }
    
    if ( spc->inside(iso.mov) )
    {
        if ( placement == PLACE_INSIDE || placement == PLACE_ALL_INSIDE )
            return 1;
    }
    else
    {
        if ( placement == PLACE_OUTSIDE )
            return 1;
    }
    
    return 0;
}


bool all_points_inside(ObjectList const& objs, Space const* spc)
{
    for ( Object * i : objs )
    {
        Mecable * mec = Simul::toMecable(i);
        if ( mec && ! mec->allPointsInside(spc) )
            return false;
    }
    return true;
}

/**
 This would usually create ONE object of type 'pp', placed according to `opt`
 */
ObjectList Interface::new_object(ObjectSet* set, Property const* prp, Glossary& opt)
{
    ObjectList objs;
    long max_trials = 1024;
    opt.set(max_trials, "nb_trials");
    long nb_trials = max_trials;
    Glossary::dict_type<PlacementType> keys{
        {"off",       PLACE_NOT},
#if BACKWARD_COMPATIBILITY < 50
        {"none",       PLACE_NOT},
#endif
        {"anywhere",   PLACE_ANYWHERE},
        {"inside",     PLACE_INSIDE},
        {"all_inside", PLACE_ALL_INSIDE},
        {"outside",    PLACE_OUTSIDE},
        {"surface",    PLACE_EDGE}};
    
    while ( --nb_trials >= 0 )
    {
        objs = set->newObjects(prp, opt);
        
#ifndef NDEBUG
        // check for `nullptr` in list, which should not happen:
        if ( objs.count(nullptr) )
        {
            std::clog << "Cytosim void slots in newObjects(" << prp->name() << ")\n";
            objs.remove_pack(nullptr);
        }
#endif
        
        // early bailout for immobile objects:
        if ( objs.size()==1 && !objs[0]->mobile() )
            break;
        
        PlacementType placement = PLACE_INSIDE;
        opt.set(placement, "placement", keys);
        if ( placement == PLACE_NOT )
            break;
        
        // find possible position & rotation:
        Isometry iso;
        if ( find_placement(iso, opt, placement) )
        {
            // place object at this position:
            for ( Object * obj : objs )
                obj->move(iso);
            // special case for which we check all vertices:
            bool okay = true;
            if ( placement == PLACE_ALL_INSIDE )
            {
                std::string str;
                Space const* spc = sim_->spaces.master();
                if ( opt.set(str, "placement", 1) )
                    spc = sim_->findSpace(str);
                okay = all_points_inside(objs, spc);
            }
            if ( okay )
                break;
        }
        else
        {
            // no suitable placement found, delete new objects:
            for ( Object* i : objs )
                if ( ! i->linked() )
                    delete(i);
            objs.clear();
            continue;
        }
        /*
         objects that were just created by newObjects() are not yet linked and
         will be deleted. Older objects will be moved back to their original position
         */
        iso.inverse();
        for ( Object* i : objs )
        {
            if ( ! i->linked() )
                delete(i);
            else
                i->move(iso);
        }
        objs.clear();
    }
    
    if ( objs.empty() )
    {
        std::string name = prp ? prp->name() : "object";
        if ( max_trials > 1 )
            Cytosim::log("could not place `", name, "' after ", max_trials, " trials\n");
        return objs;
    }

    // optionally mark the objects:
    ObjectMark mk = 0;
    if ( opt.value_is("mark", 0, "random") )
        mk = RNG.pint32();
    if ( mk || opt.set(mk, "mark") )
    {
        for ( Object * i : objs )
            i->mark(mk);
    }
    
    // optionally enlist buddies:
    std::string str;
    if ( opt.set(str, "buddy") )
    {
        Mecable * bud = sim_->pickMecable(str);
        if ( !bud )
            throw InvalidParameter("could not find buddy `"+str+"'");
        for ( Object * i : objs )
        {
            Mecable * mec = Simul::toMecable(i);
            if ( mec ) mec->enlist(bud);
        }
    }
    
    // set identity if specified
    ObjectID id = 0;
    if ( opt.set(id, "identity") )
    {
        if ( set->identifyObject(id) )
            throw InvalidParameter("identity "+std::to_string(id)+" is already assigned");
        objs.front()->setIdentity(id);
    }

    // translation after placement
    Vector vec;
    if ( opt.set(vec, "translation") )
        ObjectSet::translateObjects(objs, vec);
    
    //std::clog << "new_object " << objs.size() << " " << prp->name() << "\n";
    return objs;
}


/**
 Create `cnt` objects of type 'name', according to specifications.
 It is possible to make an object without an associated Property
 */
ObjectList Interface::execute_new(std::string const& cat, std::string const& name, Glossary& opt, size_t cnt)
{
    ObjectList res;
    ObjectSet * set = nullptr;
    Property const* prp = sim_->properties.find(name);
    
    if ( cat.empty() && prp )
    {
        set = sim_->findSet(prp->category());
        if ( !set )
            throw InvalidSyntax("could not determine the class of `"+name+"'");
    }
    else if ( cat.empty() )
        throw InvalidSyntax("could not determine the class of `"+name+"'");
    else
    {
        set = sim_->findSet(cat);
        if ( !set )
            throw InvalidSyntax("undefined class `"+cat+"'");
    }
    
    size_t amount = set->size();
    /// allow to set a desired number of objects:
    size_t target = 0;
    if ( opt.set(target, "nb_objects") )
    {
        if ( target < amount )
        {
            ObjectList objs = set->collect(amount-target);
            set->eraseObjects(objs);
            return res;
        }
        // create enough objects to reach target:
        cnt = target - amount;
    }

    // syntax sugar: distribute objects regularly between two points
    std::string var = "position_range";
#if BACKWARD_COMPATIBILITY < 60
    if ( !opt.has_key(var) )
        var = "range";
#endif
    if ( opt.has_key(var) )
    {
        Vector A, B;
        if ( !opt.set(A, var) || !opt.set(B, var, 1) )
            throw InvalidParameter("two vectors need to be defined by `range'");
        if ( opt.has_key("position") )
            throw InvalidParameter("cannot specify `position' if `range' is defined");
        Vector dAB = ( B - A ) / std::max(1UL, cnt-1);
        
        for ( size_t n = 0; n < cnt; ++n )
        {
            opt.define("position", A + n * dAB);
            res.append(new_object(set, prp, opt));
        }
    }
    // syntax sugar: positions specified for multiple objects
    else if ( opt.num_values("positions") > 0 )
    {
        for ( size_t n = 0; n < cnt; ++n )
        {
            opt.define("position", opt.least_used_value("positions"));
            res.append(new_object(set, prp, opt));
        }
    }
    else if ( set == &sim_->singles && opt.has_key("multi_base") )
    {
        // particular case: distribute singles onto many beads (31.01.2023):
        SingleProp const* sp = static_cast<SingleProp const*>(prp);
        std::string str;
        if ( opt.set(str, "multi_base") )
            sim_->singles.distributeWrists(res, sp, cnt, str);
    }
    else
    {
        // syntax sugar: specify the positions of the Fiber's ends
        if ( opt.has_key("position_ends") )
        {
            Vector A, B;
            if ( !opt.set(A, "position_ends") || !opt.set(B, "position_ends", 1) )
                throw InvalidParameter("two vectors need to be defined by `position_ends'");
            opt.define("length",    (A-B).norm());
            opt.define("position",  (A+B)*0.5);
            opt.define("direction", (B-A).normalized());
        }
        
        // normal pathway:
        for ( size_t n = 0; n < cnt; ++n )
            res.append(new_object(set, prp, opt));
    }
    //hold();
    
    /*
     Because the objects in ObjectList are not necessarily all of the same class,
     for example a Single can be created along with a Fiber in FiberSet::newObjects,
     we call here sim_->add() rather than directly set->add()
     */
    sim_->add(res);

    size_t required = 0;
    if ( opt.set(required, "required") )
    {
        size_t created = set->size() - amount;
        if ( created < required )
        {
            std::cerr << "created  = " << created << '\n';
            std::cerr << "required = " << required << '\n';
            throw InvalidParameter("could not create enough `"+name+"'");
        }
    }

    VLOG("+NEW " << cat << " `" << name << "' made " << set->size()-amount << " objects (total " << sim_->nbObjects() << ")");
    return res;
}


//------------------------------------------------------------------------------
/**
 Creates `cnt` objects of class `name`.
 The objects are distributed at the specified position in the given Space,
 with random orientations.
 
 This is meant to replace execute_new(cat, name, opt, cnt), when no fancy
 option were specified to the command.
 
 By default, this will only try one position for each object, and so it may
 create fewer objects than 'cnt'. 
 */
ObjectList Interface::execute_new(std::string const& name, size_t cnt, 
                                  Space const* spc, std::string const& position)
{
    Property const* prp = sim_->properties.find_or_die(name);
    ObjectSet * set = sim_->findSet(prp->category());
    if ( !set )
        throw InvalidSyntax("could not determine the class of `"+name+"'");

    Glossary opt;
    ObjectList res(cnt);
    set->reserve(cnt+set->inventory_.highest());
    for ( size_t n = 0; n < cnt; ++n )
    {
        ObjectList objs = set->newObjects(prp, opt);
        
        if ( objs.empty() )
            throw InvalidSyntax("could not create any `"+name+"'");
        
        Object * obj = nullptr;
        if ( objs.size() == 1 )
            obj = objs[0];

        if ( spc )
        {
            Vector pos;
            if ( position.empty() )
                pos = spc->place();
            else
                pos = Cytosim::readPosition(position, spc);
            
            if ( !pos.valid() )
            {
                objs.destroy();
                continue;
            }
            
            if ( obj )
            {
                // here the random rotation is only generated if needed:
                switch ( obj->mobile() )
                {
                    case 2: obj->rotate(Rotation::randomRotation()); break;
                    case 3: obj->rotate(Rotation::randomRotation());
                    case 1: obj->translate(pos);
                }
            }
            else
            {
                Isometry iso(Rotation::randomRotation(), pos);
                for ( Object * o : objs )
                    o->move(iso);
            }
        }
        
        /* Call sim_->add(), in case the list might contain heterogenous objects */
        if ( obj ) {
            set->add(obj);
            res.push_back(obj);
        } else {
            sim_->add(objs);
            res.append(objs);
        }
    }
    
    VLOG("-NEW " << cnt << " `" << name << "' at `" << position << "'");
    //hold();
    return res;
}

//------------------------------------------------------------------------------
#pragma mark -

/// holds a set of criteria used to select Objects
class Filter
{
public:

    Space const* ins;
    Space const* ous;
    Property const* prp;
    ObjectMark mrk;
    unsigned st1;
    unsigned st2;

    /// initialize
    void reset()
    {
        mrk = 0;
        st1 = ~0U;
        st2 = ~0U;
        prp = nullptr;
        ins = nullptr;
        ous = nullptr;
    }
    
    void set(Simul* sim, Property* pp, Glossary& opt)
    {
        prp = pp;
        std::string str;
        if ( opt.set(str, "position", 1) )
        {
            Space const* spc = sim->spaces.master();
            spc = sim->findSpace(str);
            if ( !spc )
                throw InvalidSyntax("unknown Space `"+str+"'");
            opt.set(str, "position");
            if ( str == "inside" )
                ins = spc;
            else if ( str == "outside" )
                ous = spc;
        }
        
        opt.set(mrk, "mark");
        opt.set(st1, "state1", "stateP") || opt.set(st1, "state");
        opt.set(st2, "state2", "stateM") || opt.set(st1, "state", 1);
    }

    /// initialize: will pass anything
    Filter()
    {
        reset();
    }
    
    /// initialize
    Filter(Simul* sim, Property* pp, Glossary& opt)
    {
        reset();
        set(sim, pp, opt);
    }

    /// return `true` if given object fulfills all the conditions specified
    bool pass(Object const* obj) const
    {
        if ( mrk > 0 && obj->mark() != mrk )
            return false;
        if ( ins && ins->outside(obj->position()) )
            return false;
        if ( ous && ous->inside(obj->position()) )
            return false;
        if ( prp && obj->property() != prp )
            return false;
        if ( st1 != ~0U )
        {
            if ( obj->tag()==Single::TAG && static_cast<Single const*>(obj)->attached() != st1 )
                return false;
            if ( obj->tag()==Couple::TAG && static_cast<Couple const*>(obj)->attached1() != st1 )
                return false;
            if ( obj->tag()==Fiber::TAG && static_cast<Fiber const*>(obj)->endStateP() != st1 )
                return false;
        }
        if ( st2 != ~0U )
        {
            if ( obj->tag()==Single::TAG )
                throw InvalidParameter("to select Single, `state2' is irrelevant");
            if ( obj->tag()==Couple::TAG && static_cast<Couple const*>(obj)->attached2() != st2 )
                return false;
            if ( obj->tag()==Fiber::TAG && static_cast<Fiber const*>(obj)->endStateM() != st2 )
                return false;
        }
        return true;
    }
};


bool pass_filter(Object const* obj, void const* val)
{
    return static_cast<Filter const*>(val)->pass(obj);
}


void Interface::execute_delete(std::string const& name, Glossary& opt, size_t cnt)
{
    if ( name == "objects" )
        return sim_->eraseObjects();
    Property * pp = nullptr;
    ObjectSet * set = findClass(name, pp);
    if ( set )
    {
        real rate = 0;
        if ( !pp && opt.empty() )
        {
            // all objects from a class are deleted:
            set->erase();
        }
        else if ( opt.set(rate, "rate_each") )
        {
            if ( cnt != ~0UL )
                throw InvalidParameter("an object count should not be specified with `rate_each'");
            ObjectList objs;
            if ( opt.num_keys() == 1 )
            {
                if ( pp )
                    objs = set->collect(match_property, pp);
                else
                    objs = set->collect();
            }
            else
            {
                Filter filter(sim_, pp, opt);
                objs = set->collect(pass_filter, &filter);
            }
            cnt = RNG.poisson(objs.size()*rate*sim_->time_step());
            if ( cnt > 0 )
            {
                //std::clog << "deleting " << cnt << " of " << objs.size() << " " << name << "\n";
                objs.shuffle_truncate(cnt);
                set->eraseObjects(objs);
            }
        }
        else
        {
            // multiple objects are specified, by matching the conditions
            Filter filter(sim_, pp, opt);
            ObjectList objs = set->collect(pass_filter, &filter, cnt);
            if ( objs.size() > 0 )
                set->eraseObjects(objs);
            else
                Cytosim::warn("found no `", name, "' to delete\n");
        }
    }
    else
    {
        // a single object is specified
        Object * obj = findObject(name, pp);
        if ( obj )
        {
            sim_->remove(obj);
            delete(obj);
        }
        else
            throw InvalidParameter("invalid object specified in command `delete'");
    }
}


/**
 This moves objects to a new position, or translates them by given vector
 */
size_t Interface::execute_move(std::string const& name, Glossary& opt, size_t cnt)
{
    Property * pp = nullptr;
    Space const* spc = sim_->spaces.master();
    ObjectList objs;
    bool detach = false;
    opt.set(detach, "detach");

    ObjectSet * set = findClass(name, pp);
    if ( set )
    {
        // multiple objects are specified, by matching the conditions
        Filter filter(sim_, pp, opt);
        objs = set->collect(pass_filter, &filter, cnt);
        if ( objs.empty() )
            throw InvalidParameter("no object found for command `move'");
    }
    else
    {
        // a single object is specified
        Object * obj = findObject(name, pp);
        if ( !obj )
            throw InvalidParameter("invalid object specified in command `move'");
        objs.push_back(obj);
    }
    
    VLOG("-MOVE " << objs.size() << " " << name << " "+opt.to_string());

    Vector vec;
    std::string str;
    for ( Object * obj : objs )
    {
        if ( detach )
            sim_->singles.detachWrists(obj);
        if ( opt.set_from_least_used_value(str, "position") )
        {
            vec = Cytosim::findPosition(str, spc);
            obj->setPosition(vec);
        }
        else if ( opt.set_from_least_used_value(str, "translation") )
        {
            vec = Cytosim::findPosition(str, spc);
            obj->translate(vec);
        }
#if NEW_SOLID_CLAMP
        else if ( opt.set(str, "clamp") )
        {
            Solid * sol = Solid::toSolid(obj);
            if ( sol )
            {
                vec = Cytosim::findPosition(str, spc);
                real val = sol->clampStiffness();
                opt.set(val, "clamp", 1);
                sol->setClamp(vec, val);
            }
            else
                throw InvalidParameter("invalid Solid for command `move'");
        }
#endif
        else
            throw InvalidParameter("unspecified position for command `move'");
    }
    return objs.size();
}


/**
 This can only mark one class of objects
 */
void Interface::execute_mark(std::string const& name, Glossary& opt, size_t cnt)
{
    Property * pp = nullptr;
    ObjectSet * set = findClass(name, pp);
    
    ObjectMark mk = 0;
    if ( ! opt.set(mk, "mark") )
        throw InvalidParameter("mark must be specified for command `mark'");
    opt.clear("mark");

    if ( set )
    {
        // multiple objects are specified, by matching the conditions
        Filter filter(sim_, pp, opt);
        ObjectList objs = set->collect(pass_filter, &filter, cnt);
        if ( objs.size() > 0 )
            ObjectSet::markObjects(objs, mk);
        else
            Cytosim::warn("found no `", name, "' to mark\n");
    }
    else
    {
        // a single object is specified
        Object * obj = findObject(name, pp);
        if ( !obj )
            throw InvalidParameter("invalid object specified in command `mark'");
        obj->mark(mk);
    }
}


void Interface::execute_cut(std::string const& name, Glossary& opt, size_t cnt)
{
    Vector n(1,0,0);
    real a = 0;
    real len = 0;
    
    if ( opt.value_is("plane", 0, "random") )
        n = Vector::randU();
    else
        opt.set(n, "plane");
    opt.set(a, "plane", 1);
    opt.set(len, "min_length");

    Glossary::dict_type<state_t> keys{{"white",     STATE_WHITE},
                                      {"green",     STATE_GREEN},
                                      {"yellow",    STATE_YELLOW},
                                      {"orange",    STATE_ORANGE},
                                      {"red",       STATE_RED},
                                      {"static",    STATE_WHITE},
                                      {"grow",      STATE_GREEN},
                                      {"growing",   STATE_GREEN},
                                      {"shrink",    STATE_RED},
                                      {"shrinking", STATE_RED}};

    state_t stateP = STATE_RED, stateM = STATE_WHITE;
    opt.set(stateP, "new_end_state", 0, keys);
    opt.set(stateM, "new_end_state", 1, keys);
    
    Property * pp = nullptr;
    ObjectSet * set = findClass(name, pp);
    if ( set != &sim_->fibers )
        throw InvalidSyntax("only `cut fiber' is supported");
    
    if ( pp )
    {
        ObjectList objs = set->collect(match_property, pp, cnt);
        VLOG("-CUT " << objs.size() << " " << name << " PLANE (" << n << ").x = " << -a);
        sim_->fibers.planarCut(objs, n, a, stateP, stateM, len);
    }
    else if ( cnt < sim_->fibers.size() )
    {
        ObjectList objs = set->collect(cnt);
        VLOG("-CUT " << objs.size() << " " << name << " PLANE (" << n << ").x = " << -a);
        sim_->fibers.planarCut(objs, n, a, stateP, stateM, len);
    }
    else
    {
        VLOG("-CUT all fibers PLANE (" << n << ").x = " << -a);
        sim_->fibers.planarCut(n, a, stateP, stateM, len);
    }
}


void Interface::execute_equilibrate(std::string const& name, Glossary& opt)
{
    int only_intersections = 0;
    opt.set(only_intersections, "only_intersections");
    
    if ( only_intersections )
    {
        if ( name == "couple" )
            sim_->couples.bindToIntersections();
        else
        {
            Property * pp = sim_->properties.find_or_die(name);
            if ( pp->category() == "couple" )
                sim_->couples.bindToIntersections(static_cast<CoupleProp*>(pp));
            else
                throw InvalidSyntax("can only equilibrate (only_intersections=1) `couple' or a class name");
        }
    }
    else
    {
        if ( name == "single" )
            sim_->singles.equilibrate();
        else if ( name == "couple" )
            sim_->couples.equilibrate();
        else
        {
            Property * pp = sim_->properties.find_or_die(name);
            if ( pp->category() == "single" )
                sim_->singles.equilibrate(sim_->fibers, static_cast<SingleProp const*>(pp));
            else if ( pp->category() == "couple" )
                sim_->couples.equilibrate(sim_->fibers, static_cast<CoupleProp const*>(pp));
            else
                throw InvalidSyntax("can only equilibrate `single', `couple' or a class name");
        }
    }
    
    VLOG("-CONNECT (" << name << ")");
}

//------------------------------------------------------------------------------
#pragma mark -

/** Using a static frame counter */
static void reportCPUtime(real t)
{
    static size_t frm = 1;
    static time_t nxt = 0;
    time_t now = TimeDate::seconds_since_1970();
    if ( now > nxt )
    {
        nxt = (nxt>0?nxt:now) + 3600;
        Cytosim::log("% ",TimeDate::date_string(),"\n");
    }
    static double clk = 0;
    double cpu = double(clock()) / CLOCKS_PER_SEC;
    Cytosim::log.print("F%-6lu  %7.2fs   CPU %10.3fs  %10.0fs\n", frm, t, cpu-clk, cpu);
    clk = cpu;
    ++frm;
}


template < Interface::SimulFuncPtr FUNC >
inline void Interface::step_simul()
{
    while ( sim_->incomplete() )
    {
        hold();
        //fprintf(stderr, "> step @%.12e\n", sim_->time());
        (sim_->*FUNC)();
        sim_->steps();
    }
}

/**
 Perform simulation steps. The accepted Syntax is:
 
     run POSITIVE_INTEGER SIMUL_NAME
     {
        duration   = POSITIVE_REAL
        solve      = SOLVE_MODE
        nb_frames  = INTEGER, ( CODE )
        prune      = BOOL
     }
 
 or
 
     run SIMUL_NAME
     {
        nb_steps   = POSITIVE_INTEGER
        ...
     }

 or, without specifying the Name of the Simul:
 
     run [POSITIVE_INTEGER] all simul
     {
        ...
     }

 
 The associated block can specify these parameters:
 
 Parameter    | Default | Description                                          |
 -------------|---------|-------------------------------------------------------
 `nb_steps`   |  1      | number of simulation steps
 `duration`   |  -      | when specified, `nb_steps` is set to `std::ceil(duration/time_step)`
 `solve`      |  `on`   | Define the type of method used for the mechanics
 `nb_frames`  |  0      | number of states written to trajectory file
 `prune`      |  `true` | Print only parameters that are different from default
 
 Set `nb_frames = 2` to save the initial and last time point of the run.
 Set `nb_frames = 1` to save only the last time point of the run.
 The parameter `solve` can be used to select alternative mechanical engines.
 The Monte-Carlo parts of the simulation is always done, which includes
 fiber assembly dynamics, binding/unbinding and diffusion of molecules.
 
 `solve`      | Result                                                         |
 -------------|-----------------------------------------------------------------
 `off`        | Objects are immobile and mechanical equations are not solved.
 `on`         | The mechanics is solved and applied to the objects (default).
 `auto`       | Same as 'on' but preconditionning method is set automatically.
 `force`      | Calculate forces given the current object's positions.
 `half`       | Solve mechanical system and calculate forces but do not apply movements.
 `horizontal` | The mechanics is solved only allowing motion in the X-direction.
 `flux`       | Fibers are translated at `flux_speed` according to their orientation.
 
 */
void Interface::execute_run(real sec, Glossary& opt, bool do_write)
{
    int solve = 1;
    int binary = 1;
    long frames = 0;
    bool prune = true;
    bool has_code = false;
    std::string code;
    
#if BACKWARD_COMPATIBILITY < 50
    // create an Event if 'event' is specified within the 'run' command:
    Event * evt = nullptr;
    if ( opt.has_key("event") )
    {
        evt = new Event();
        opt.set(evt->rate, "event");
        opt.set(evt->activity, "event", 1);
        evt->reload(sim_->time());
        sim_->events.add(evt);
    }
#endif
    opt.set(solve, "solve", {{"off",0}, {"on",1}, {"auto",2}, {"force", 3},
        {"uniaxial",4}, {"half",7} });
    opt.set(prune,  "prune");
    opt.set(binary, "binary");
    opt.set(frames, "nb_frames");
    has_code = opt.set(code, "nb_frames", 1);
    
    do_write &= ( frames > 0 );
    sim_->prepare();

    if ( do_write && frames > 1 )
    {
        // write initial state:
        sim_->writeProperties(prune);
        sim_->writeObjects(sim_->prop.system_file, true, binary);
    }
    
    VLOG("+RUN START +" << sec);
    double tau = sim_->time_step();
    // limit to one frame per time_step:
    long max = std::max(std::min(frames, std::lround(sec/tau)), 1L);
    // subtract half a time_step, to ensure we finish exactly on time!
    double start = sim_->time() - 0.5 * tau;
    double delta = sec / double(max);
    
    for ( int frm = 1; frm <= max; ++frm )
    {
        sim_->stop_at(start + delta * frm);
        switch ( solve )
        {
            case 0: step_simul<&Simul::solve_not>(); break;
            case 1: step_simul<&Simul::solve_meca>(); break;
            case 2: step_simul<&Simul::solve_auto>(); break;
            case 3: step_simul<&Simul::solve_force>(); break;
            case 4: step_simul<&Simul::solve_uniaxial>(); break;
            case 7: step_simul<&Simul::solve_half>(); break;
        }

        if ( do_write )
        {
            sim_->relax();
            if ( has_code )
                sim_->perform(code);
            sim_->writeObjects(sim_->prop.system_file, true, binary);
            reportCPUtime(sim_->time());
            sim_->sMeca.doNotify = 2;  // to print convergence parameters
            sim_->unrelax();
        }
        
        // check if ending the simulation was requested:
        if ( sim_->should_end() )
        {
            sim_->end_never(); // cleanup for next run
            break;
        }
    }

#if BACKWARD_COMPATIBILITY < 50
    if ( evt )
    {
        sim_->events.remove(evt);
        delete(evt);
    }
#endif
    
    sim_->relax();
    if ( frames )
    {
        sim_->writeProperties(prune);
        if ( frames < 0 )
            sim_->writeObjects(sim_->prop.system_file, true, binary);
    }
    
    VLOG("+RUN END t: " << sim_->time());
    hold();
}


/**
 Advance simulation, without any option, by alternating `step` and `solve`
*/
void Interface::execute_run(real sec)
{
    VLOG("-RUN START " << sec);
    sim_->prepare();
    // subtract half a time_step, to ensure we finish exactly on time!
    sim_->stop_at(sim_->time() + sec - 0.5 * sim_->time_step());

    while ( sim_->incomplete() )
    {
        hold();
        sim_->solve_meca();
        sim_->steps();
    }
    
    // reset termination time for next run:
    if ( sim_->should_end() )
        sim_->end_never();

    sim_->relax();
    VLOG("-RUN END");
    hold();
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 Import a simulation snapshot from a trajectory file
 
 The frame to be imported can be specified as an option: `frame=INTEGER`:
 
     import objects objects.cmi { frame = 10 }
 
 By default, this will replace the simulation state by the one loaded from file.
 To add the loaded objects to the simulation without deleting the current world,
 you should specify `append=1`:
 
     import objects objects.cmi { append = 1 }
 
 This will work however only if the ID of the objects are distinct, ie. are not
 already in use in the current world.
 In the examples, the `cmi` extension are like `cmo`. The extension is ignored.
 */
void Interface::execute_import(std::string const& file, std::string const& what, Glossary& opt)
{
    bool append = false;
    ObjectSet * subset = nullptr;
    
    if ( what != "all" && what != "objects" )
    {
        append = true;
        subset = sim_->findSet(what);
        if ( !subset )
            throw InvalidIO("expected a class to be specified (import fiber FILE)");
    }

    Inputter in(DIM, file.c_str(), true);

    if ( ! in.good() )
        throw InvalidIO("Could not open file `"+file+"'");
    
    size_t cnt = 0, frm = 0;

    opt.set(frm, "frame");
    opt.set(append, "append");

    VLOG("-IMPORT frame " << frm << " from " << file);

    while ( in.good() )
    {
        if ( append )
        {
            double t = sim_->time();
            sim_->reloadObjects(in, 0, subset);
            sim_->set_time(t);
        }
        else
            sim_->reloadObjects(in, 1, subset);
        if ( cnt >= frm )
            break;
        ++cnt;
    }
    
    if ( cnt < frm )
        throw InvalidIO("Could not import requested frame");
    
#if ( 0 )
    //unfinished code to mark imported objects
    int mrk;
    if ( opt.set(mrk, "mark") )
    {
        ObjectSet::markObjects(objs, mrk);
    }
#endif
    
    // set time
    double t;
    if ( opt.set(t, "time") )
        sim_->set_time(t);
}


/**
 see Parser::parse_export
 */
void Interface::execute_export(std::string const& name, std::string const& what, Glossary& opt)
{
    VLOG("-EXPORT " << what << " to " << name);

    // here '*' designates the standard output:
    if ( what == "all" || what == "objects" )
    {
        if ( name != "*" )
        {
            int binary = 1;
            bool append = true;
            opt.set(append, "append");
            opt.set(binary, "binary");
            sim_->writeObjects(name, append, binary);
        }
        else
        {
            Outputter out(stdout, false);
            sim_->writeObjects(out);
        }
    }
    else if ( what == "properties" )
    {
        bool prune = true;
        opt.set(prune, "prune");
        std::ofstream ofs;
        std::ostream out(std::cout.rdbuf());
        // a STAR designates the standard output:
        if ( name != "*" )
        {
            std::ios_base::openmode mode = std::ios_base::app;
            opt.set(mode, "append", {{"0", std::ios_base::out}, {"1", std::ios_base::app}});
            ofs.open(name.c_str(), mode);
            if ( ofs.is_open() )
                out.rdbuf(ofs.rdbuf());
            else
                throw InvalidIO("cannot open `"+name+"' for export");
        }
        sim_->writeProperties(out, prune);
    }
    else
        throw InvalidIO("only `objects' or `properties' can be exported");
}


/**
 see Parser::parse_report
 */
void Interface::execute_report(std::string const& name, std::string const& what, Glossary& opt)
{
    VLOG("-REPORT " << what << " to " << name);
    
    std::ofstream ofs;
    std::ostream out(std::cout.rdbuf());

    // a STAR designates the standard output:
    if ( name != "*" )
    {
        std::ios_base::openmode mode = std::ios_base::app;
        opt.set(mode, "append", {{"0", std::ios_base::out}, {"1", std::ios_base::app}});
        ofs.open(name.c_str(), mode);
        if ( ofs.is_open() )
            out.rdbuf(ofs.rdbuf());
    }
    int ver = 1;
    opt.set(ver, "verbose");

    sim_->mono_report(out, what, opt, ver);
}


void Interface::execute_call(std::string& str, Glossary& opt)
{
#if BACKWARD_COMPATIBILITY <= 50
    if ( str == "couple:equilibrate" )
        sim_->couples.equilibrate();
    else if ( str == "couple:connect" )
        sim_->couples.bindToIntersections();
    else if ( str == "single:equilibrate" )
        sim_->singles.equilibrate();
    else
#endif
    if ( str == "custom0" )
        sim_->custom0(opt);
    else if ( str == "custom1" )
        sim_->custom1(opt);
    else if ( str == "custom2" )
        sim_->custom2(opt);
    else if ( str == "custom3" )
        sim_->custom3(opt);
    else if ( str == "custom4" )
        sim_->custom4(opt);
    else if ( str == "custom5" )
        sim_->custom5(opt);
    else if ( str == "custom6" )
        sim_->custom6(opt);
    else if ( str == "custom7" )
        sim_->custom7(opt);
    else if ( str == "custom8" )
        sim_->custom8(opt);
    else if ( str == "custom9" )
        sim_->custom9(opt);
    else
        throw InvalidSyntax("unknown command `"+str+"' called");
}


void Interface::execute_dump(std::string const& path, unsigned mode)
{
    sim_->sMeca.doNotify = 1;
    sim_->solve_half();
    
    size_t dim = sim_->sMeca.dimension();
    Cytosim::log("Cytosim is dumping a system of size ", dim, " in `", path, "'...");
    int cwd = FilePath::change_dir(path, true);
    
    if ( mode & 1 ) sim_->sMeca.saveSystem();
    if ( mode & 2 ) sim_->sMeca.dumpSystem();
    if ( mode & 4 ) sim_->sMeca.exportSystem();
    if ( mode & 8 ) sim_->sMeca.saveMatrixBitmaps("");
    if ( mode & 16 ) sim_->sMeca.saveConnectivityBitmap();

    FilePath::change_dir(cwd);
    Cytosim::log("done\n");
}
