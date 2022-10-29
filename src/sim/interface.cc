// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include <fstream>
#include <cstdio>

#include "cymdef.h"
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

void Interface::erase_simul(bool arg) const
{
    return sim_->eraseObjects(arg);
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
    
    /* We do not allow for using the class name to name a property,
    as this should create confusion in the config file */
    
    Property* pp = sim_->newProperty(cat, name, def);
    
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


void Interface::change_simul_property(Glossary& opt)
{
    change_property(&sim_->prop, opt);
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

using StreamFunc::has_trail;

/// report a warning if some text was ignored
void warn_trail(std::istream& is)
{
    std::string str;
    std::streampos pos = is.tellg();
    std::getline(is, str);
    throw InvalidSyntax("unexpected `"+str+"' in `"+StreamFunc::extract_line(is, pos)+"'");
}

/**
 Define a placement = ( position, orientation ) from the parameters set in `opt'
 */
void Interface::read_placement(Isometry& iso, Glossary& opt)
{
    std::string str;
    
    Space const* spc = sim_->spaces.master();
    
    // Space specified as second argument to 'position'
    if ( opt.set(str, "position", 1) )
        spc = sim_->findSpace(str);
    
    // Position
    if ( opt.set(str, "position") )
        iso.mov = Movable::readPosition(str, spc);
    else if ( spc )
        iso.mov = spc->place();
    
    // Rotation applied before the translation
    if ( opt.set(str, "direction") )
    {
        std::istringstream iss(str);
        Vector vec = Movable::readDirection(iss, iso.mov, spc);
        if ( has_trail(iss) ) warn_trail(iss);
        iso.rot = Rotation::randomRotationToVector(vec);
    }
    else if ( opt.set(str, "rotation") )
    {
        std::istringstream iss(str);
        iso.rot = Movable::readRotation(iss);
        if ( has_trail(iss) ) warn_trail(iss);
    }
    else if ( opt.set(str, "orientation") )
    {
        std::istringstream iss(str);
        iso.rot = Movable::readOrientation(iss, iso.mov, spc);
        if ( has_trail(iss) ) warn_trail(iss);
    }
    else
        iso.rot = Rotation::randomRotation();
    
    // Second rotation applied after the translation
    if ( opt.set(str, "orientation", 1) )
    {
        std::istringstream iss(str);
        Rotation rot = Movable::readOrientation(iss, iso.mov, spc);
        if ( has_trail(iss) ) warn_trail(iss);
        iso.rotate(rot);
    }
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
bool Interface::find_placement(Isometry& iso, Glossary& opt, int placement, size_t& nb_trials)
{
    std::string str;

    Space const* spc = sim_->spaces.master();
    if ( opt.set(str, "placement", 1) )
        spc = sim_->findSpace(str);
    
    // define a condition:
    bool has_condition = opt.set(str, "placement", 2);

    while ( nb_trials > 0 )
    {
        --nb_trials;
        // generate a new position:
        read_placement(iso, opt);

        // check all conditions:
        bool condition = true;
        if ( has_condition )
        {
            Evaluator evaluator{{"X", iso.mov.x()}, {"Y", iso.mov.y()}, {"Z", iso.mov.z()},
                                {"R", iso.mov.norm()}, {"P", RNG.preal()}};
            try {
                condition = ( 0 != evaluator.eval(str) );
            }
            catch( Exception& e ) {
                e.message(e.message()+" in `"+str+"'");
                throw;
            }
        }
        
        if ( condition )
        {
            if ( !spc || placement == PLACE_ANYWHERE )
                return 0;
            
            if ( placement == PLACE_EDGE )
            {
                iso.mov = spc->project(iso.mov);
                return 0;
            }
            
            if ( spc->inside(iso.mov) )
            {
                if ( placement == PLACE_INSIDE || placement == PLACE_ALL_INSIDE )
                    return 0;
            }
            else
            {
                if ( placement == PLACE_OUTSIDE )
                    return 0;
            }
        }
    }
    
    //Cytosim::warn << "could not fulfill position=`" + opt.value("position") + "'\n";
    //throw InvalidParameter("could not fulfill position=`" + opt.value("position") + "'");
    return 1;
}


bool all_objects_inside(ObjectList const& objs, Space const* spc)
{
    for ( Object * i : objs )
    {
        Mecable * mec = Simul::toMecable(i);
        if ( mec && ! mec->allInside(spc) )
            return false;
    }
    return true;
}

/**
 This would usually create ONE object of type 'name', placed according to `opt`
 */
ObjectList Interface::new_object(ObjectSet* set, Property const* pp, Glossary& opt)
{
    size_t nb_trials = 1000;
    opt.set(nb_trials, "nb_trials");
    ObjectList objs;
    
    while ( nb_trials > 0 )
    {
        objs = set->newObjects(pp, opt);
        
#ifndef NDEBUG
        // check for `nullptr` in list, which should not happen:
        if ( objs.count(nullptr) )
        {
            std::clog << "Cytosim void slots in newObjects(" << pp->name() << ")\n";
            objs.remove_pack(nullptr);
        }
#endif
        
        // early bailout for immobile objects:
        if ( objs.size()==1 && !objs[0]->mobile() )
            break;
        
        PlacementType placement = PLACE_INSIDE;
        
        opt.set(placement, "placement",{{"off",       PLACE_NOT},
#if BACKWARD_COMPATIBILITY < 50
                                       {"none",       PLACE_NOT},
#endif
                                       {"anywhere",   PLACE_ANYWHERE},
                                       {"inside",     PLACE_INSIDE},
                                       {"all_inside", PLACE_ALL_INSIDE},
                                       {"outside",    PLACE_OUTSIDE},
                                       {"surface",    PLACE_EDGE}});
        
        if ( placement == PLACE_NOT )
            break;
        
        // find possible position & rotation:
        Isometry iso;
        bool err = find_placement(iso, opt, placement, nb_trials);
        if ( !err )
        {
            // place object at this position:
            ObjectSet::moveObjects(objs, iso);
            // special case for which we check all vertices:
            if ( placement == PLACE_ALL_INSIDE )
            {
                std::string str;
                Space const* spc = sim_->spaces.master();
                if ( opt.set(str, "placement", 1) )
                    spc = sim_->findSpace(str);
                if ( ! all_objects_inside(objs, spc) )
                {
                    objs.destroy();
                    continue;
                }
            }
            break;
        }
    }
    
    if ( objs.empty() )
    {
        std::string name = pp ? pp->name() : "object";
        Cytosim::log << "could not place `" << name << "' after " << nb_trials << " trials\n";
        return objs;
    }

    // optionally mark the objects:
    ObjectMark mk = 0;
    if ( opt.value("mark") == "random" )
        mk = RNG.pint32();
    if ( mk || opt.set(mk, "mark") )
    {
        for ( Object * i : objs )
            i->mark(mk);
    }
    
    // translation after placement
    Vector vec;
    if ( opt.set(vec, "translation") )
        ObjectSet::translateObjects(objs, vec);
    
    // set identity if specified
    ObjectID id = 0;
    if ( opt.set(id, "identity") )
    {
        if ( set->findID(id) )
            throw InvalidParameter("identity "+std::to_string(id)+" is already assigned");
        objs.front()->setIdentity(id);
    }
    
    /* 
     Because the objects in ObjectList are not necessarily all of the same class,
     we call sim_->add() rather than directly set->add()
     */
    sim_->add(objs);
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
    Property const* pp = sim_->properties.find(name);
    
    if ( cat.empty() && pp )
    {
        set = sim_->findSet(pp->category());
        if ( !set )
            throw InvalidSyntax("could not determine the class of `"+name+"'");
    }
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
            sim_->erase(objs);
            return res;
        }
        // create enough objects to reach target:
        cnt = target - amount;
    }

    // syntax sugar, to distribute objects regularly between two points:
    if ( opt.has_key("range") )
    {
        Vector A, B;
        if ( !opt.set(A, "range") || !opt.set(B, "range", 1) )
            throw InvalidParameter("two vectors need to be defined by `range'");
        if ( opt.has_key("position") )
            throw InvalidParameter("cannot specify `position' if `range' is defined");
        Vector dAB = ( B - A ) / std::max(1UL, cnt-1);
        
        for ( size_t n = 0; n < cnt; ++n )
        {
            opt.define("position", A + n * dAB);
            res.append(new_object(set, pp, opt));
        }
    }
    else
    {
        // syntax sugar, to specify the position of the Fiber ends
        if ( opt.has_key("position_ends") )
        {
            Vector A, B;
            if ( !opt.set(A, "position_ends") || !opt.set(B, "position_ends", 1) )
                throw InvalidParameter("two vectors need to be defined by `position_ends'");
            opt.define("length",    (A-B).norm());
            opt.define("position",  (A+B)*0.5);
            opt.define("direction", (B-A).normalized());
        }

        for ( size_t n = 0; n < cnt; ++n )
            res.append(new_object(set, pp, opt));
    }
    //hold();
    
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

    VLOG("+NEW `" << name << "' made " << set->size()-amount << " objects (total " << sim_->nbObjects() << ")");
    return res;
}


//------------------------------------------------------------------------------
/**
 Creates `cnt` objects of class `name`.
 The objects are distributed uniformly within the current Space, withrandom orientations.
 
 This is meant to replace execute_new(cat, name, opt, cnt), when no option were specified
 to the command.
 */
ObjectList Interface::execute_new(std::string const& name, size_t cnt)
{
    Property const* pp = sim_->properties.find_or_die(name);
    ObjectSet * set = sim_->findSet(pp->category());
    if ( !set )
        throw InvalidSyntax("could not determine the class of `"+name+"'");

    Space const* spc = sim_->spaces.master();

    Glossary opt;
    ObjectList res;
    set->reserve(cnt);
    for ( size_t n = 0; n < cnt; ++n )
    {
        ObjectList objs = set->newObjects(pp, opt);
        
        if ( objs.empty() )
            throw InvalidSyntax("could not create object class of `"+name+"'");

        if ( spc )
        {
            // This is a very common case, where we can skip Rotation::randomRotation():
            if ( objs.size() == 1 )
            {
                Object * obj = objs[0];
                switch ( obj->mobile() )
                {
                    case 2: obj->rotate(Rotation::randomRotation()); break;
                    case 3: obj->rotate(Rotation::randomRotation());
                    case 1: obj->translate(spc->place());
                }
            }
            else
            {
                Isometry iso(spc->place(), Rotation::randomRotation());
                ObjectSet::moveObjects(objs, iso);
            }
        }
        
        /* Call sim_->add() rather than directly set->add(), because the objects
         in ObjectList are not necessarily all of the same class */
        sim_->add(objs);
        res.append(objs);
        objs.clear();
    }
    
    VLOG("-NEW " << cnt << "`" << name << "' objects");
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
    Property * pp = sim_->properties.find(name);
    ObjectSet * set = nullptr;
    if ( pp )
        set = sim_->findSet(pp->category());
    else
        set = sim_->findSet(name);
    if ( !set )
    {
        if ( name == "objects" )
        {
            sim_->eraseObjects(0);
            return;
        }
        throw InvalidSyntax("could not determine the class of `"+name+"'");
    }
    
    Filter filter(sim_, pp, opt);
    ObjectList objs = set->collect(pass_filter, &filter, cnt);
    
    if ( objs.size() > 0 )
        sim_->erase(objs);
    else
        Cytosim::warn << "found no `" << name << "' to delete\n";
}


/**
 This moves objects to position `pos`
 */
void Interface::execute_move(std::string const& name, Glossary& opt, size_t cnt)
{
    Property * pp = sim_->properties.find(name);
    ObjectSet * set = nullptr;
    if ( pp )
        set = sim_->findSet(pp->category());
    else
        set = sim_->findSet(name);
    if ( !set )
        throw InvalidSyntax("could not determine the class of `"+name+"'");
    
    Filter filter(sim_, pp, opt);
    ObjectList objs = set->collect(pass_filter, &filter, cnt);
    
    Vector pos;
    if ( opt.set(pos, "position") )
    {
        for ( Object * obj : objs )
            obj->setPosition(pos);
    }
    else if ( opt.set(pos, "translation") )
    {
        for ( Object * obj : objs )
            obj->translate(pos);
    }
}


/**
 This can only mark one class of objects
 */
void Interface::execute_mark(std::string const& name, Glossary& opt, size_t cnt)
{
    Property * pp = sim_->properties.find(name);
    ObjectSet * set = nullptr;
    if ( pp )
        set = sim_->findSet(pp->category());
    else
        set = sim_->findSet(name);
    if ( !set )
        throw InvalidSyntax("could not determine the class of `"+name+"'");

    ObjectMark mk = 0;
    if ( ! opt.set(mk, "mark") )
        throw InvalidParameter("mark must be specified for command `mark'");
    opt.clear("mark");
    
    Filter filter(sim_, pp, opt);
    ObjectList objs = set->collect(pass_filter, &filter, cnt);
    
    sim_->mark(objs, mk);
}


void Interface::execute_cut(std::string const& name, Glossary& opt, size_t cnt)
{
    Vector n(1,0,0);
    real a = 0;
    
    opt.set(n, "plane");
    opt.set(a, "plane", 1);
    
    state_t stateP = STATE_RED, stateM = STATE_GREEN;
    opt.set(stateP, "new_end_state");
    opt.set(stateM, "new_end_state", 1);
    
    Property * pp = sim_->properties.find(name);
    ObjectSet * set = nullptr;
    if ( pp )
        set = sim_->findSet(pp->category());
    else
        set = sim_->findSet(name);
    if ( !set )
        throw InvalidSyntax("could not determine the class of `"+name+"'");

    if ( set != &sim_->fibers )
        throw InvalidSyntax("only `cut fiber' is supported");

    Filter filter(sim_, pp, opt);
    ObjectList objs = set->collect(pass_filter, &filter, cnt);
    
    VLOG("-CUT PLANE (" << n << ").x = " << -a);
    sim_->fibers.planarCut(objs, n, a, stateP, stateM);
}


void Interface::execute_connect(std::string const& name, Glossary& opt)
{
    ObjectList objs;

    if ( name == "couple" )
    {
        sim_->couples.bindToIntersections();
    }
    else
    {
        Property * pp = sim_->properties.find_or_die(name);
        if ( pp->category() != "couple" )
            throw InvalidSyntax("only `bind couple' or `bind couple_class' is supported");
        sim_->couples.bindToIntersections(static_cast<CoupleProp*>(pp));
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
        Cytosim::log << "% " << TimeDate::date_string() << '\n';
    }
    static double clk = 0;
    double cpu = double(clock()) / CLOCKS_PER_SEC;
    Cytosim::log("F%-6lu  %7.2fs   CPU %10.3fs  %10.0fs\n", frm, t, cpu-clk, cpu);
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
        sim_->step();
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
    int frames = 0;
    bool prune = true;
    bool binary = true;
    
#if BACKWARD_COMPATIBILITY < 50
    // check if 'event' is specified within the 'run' command,
    // and convert to a registered Event object
    Event * event = nullptr;
    if ( opt.has_key("event") )
    {
        event = new Event();
        opt.set(event->rate, "event");
        opt.set(event->activity, "event", 1);
        event->reload(sim_->time());
        sim_->events.add(event);
    }
#endif
    opt.set(solve, "solve", {{"off",0}, {"on",1}, {"auto",2}, {"force", 3},
        {"horizontal",4}, {"flux",5}, {"half",7}, {"separate", 8} });
    opt.set(prune,  "prune");
    opt.set(binary, "binary");
    opt.set(frames, "nb_frames");
    
    do_write &= ( frames > 0 );
    sim_->prepare();

    if ( do_write )
    {
        sim_->writeProperties(prune);
        if ( sim_->prop.clear_trajectory )
        {
            sim_->prop.clear_trajectory = false;
            // write initial state:
            if ( frames > 1 )
                sim_->writeObjects(sim_->prop.system_file, false, binary);
            else
                std::remove(sim_->prop.system_file.c_str());
        }
    }
    
    VLOG("+RUN START " << sec);
    int max = std::max(frames, 1);
    // subtract half a time_step, to ensure we finish exactly on time!
    double start = sim_->time() - 0.5 * sim_->time_step();
    double delta = sec / double(max);
    
    for ( int frm = 1; frm <= max; ++frm )
    {
        sim_->prop.end_time = start + delta * frm;
        switch ( solve )
        {
            case 0: step_simul<&Simul::solve_not>(); break;
            case 1: step_simul<&Simul::solve>(); break;
            case 2: step_simul<&Simul::solve_auto>(); break;
            case 3: step_simul<&Simul::solve_force>(); break;
            case 4: step_simul<&Simul::solve_onlyX>(); break;
            case 5: step_simul<&Simul::solve_flux>(); break;
            case 7: step_simul<&Simul::solve_half>(); break;
            case 8: step_simul<&Simul::solve_separate>(); break;
        }

        if ( do_write )
        {
            sim_->relax();
            sim_->writeObjects(sim_->prop.system_file, true, binary);
            reportCPUtime(sim_->time());
            sim_->sMeca.doNotify = 2;  // to print convergence parameters
            sim_->unrelax();
        }
    }
    
#if BACKWARD_COMPATIBILITY < 50
    if ( event )
        sim_->events.erase(event);
#endif
    
    sim_->relax();
    if ( frames )
    {
        sim_->writeProperties(prune);
        if ( frames < 0 )
            sim_->writeObjects(sim_->prop.system_file, true, binary);
    }
    
    VLOG("+RUN END");
    hold();
}


/**
 Advance simulation, without any option, by alternating step() and solve()
*/
void Interface::execute_run(real sec)
{
    VLOG("-RUN START " << sec);
    sim_->prepare();
    // subtract half a time_step, to ensure we finish exactly on time!
    sim_->prop.end_time = sim_->time() + sec - 0.5 * sim_->time_step();

    while ( sim_->incomplete() )
    {
        hold();
        sim_->solve();
        sim_->step();
    }
    
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
    // we could use the 'tag' to select a certain class of object
    ObjectSet * subset = nullptr;
    
    if ( what != "all" && what != "objects" )
    {
        subset = sim_->findSet(what);
        if ( !subset )
            throw InvalidIO("expected class specifier (eg. `import all FILE' or `import fiber FILE')");
    }

    Inputter in(DIM, file.c_str(), true);

    if ( ! in.good() )
        throw InvalidIO("Could not open file `"+file+"'");
    
    bool append = false;
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
            sim_->time(t);
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
        sim_->mark(objs, mrk);
    }
#endif
    
    // set time
    double t;
    if ( opt.set(t, "time") )
        sim_->time(t);
}


/**
 see Parser::parse_export
 */
void Interface::execute_export(std::string const& name, std::string const& what, Glossary& opt)
{
    bool append = true;
    bool binary = true;
    bool prune = true;

    VLOG("-EXPORT " << what << " to " << name);

    // here '*' designates the standard output:
    if ( what == "all" || what == "objects" )
    {
        if ( name != "*" )
        {
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
        opt.set(prune, "prune");
        std::ofstream ofs;
        std::ostream out(std::cout.rdbuf());
        // a STAR designates the standard output:
        if ( name != "*" )
        {
            opt.set(append, "append");
            ofs.open(name.c_str(), append ? std::ios_base::app : std::ios_base::out);
            out.rdbuf(ofs.rdbuf());
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
        bool append = true;
        opt.set(append, "append");
        ofs.open(name.c_str(), append ? std::ios_base::app : std::ios_base::out);
        out.rdbuf(ofs.rdbuf());
    }
    int ver = 1;
    opt.set(ver, "verbose");

    sim_->mono_report(out, what, opt, ver);
}


void Interface::execute_call(std::string& str, Glossary& opt)
{
    if ( str == "equilibrate" )
        sim_->couples.equilibrate();
    else if ( str == "connect" )
        sim_->couples.bindToIntersections();
    else if ( str == "custom0" )
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
        throw InvalidSyntax("called unknown command");
}


void Interface::execute_dump(std::string const& path, int mode)
{
    Cytosim::log("Cytosim is dumping a system of size %lu in `%s'...", sim_->sMeca.dimension(), path.c_str());
    sim_->sMeca.doNotify = 1;
    sim_->solve_half();
    
    int cwd = FilePath::change_dir(path, true);
    
    if ( mode & 1 ) sim_->sMeca.saveSystem();
    if ( mode & 2 ) sim_->sMeca.dumpSystem();
    if ( mode & 4 ) sim_->sMeca.exportSystem();
    if ( mode & 8 ) sim_->sMeca.saveMatrixBitmaps();
    if ( mode & 16 ) sim_->sMeca.saveConnectivityBitmap();
    
    FilePath::change_dir(cwd);
    Cytosim::log("done\n");
}
