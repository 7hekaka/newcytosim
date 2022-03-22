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
#include "simul.h"
#include "event.h"


// Use the second definition to get some verbose reports:
#define VLOG(ARG) ((void) 0)
//#define VLOG(ARG) std::clog << ARG << '\n';

//------------------------------------------------------------------------------

Interface::Interface(Simul& s)
: simul_(s)
{
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
    
    Property* pp = simul_.newProperty(cat, name, def);
    
    if ( !pp )
        throw InvalidSyntax("failed to create property of class `"+cat+"'");
    
    pp->read(def);
    pp->complete(simul_);
    
    return pp;
}


void Interface::change_property(Property * pp, Glossary& def)
{
    pp->read(def);
    pp->complete(simul_);
    
    /*
     Specific code to make 'change space:dimension' work.
     This is needed as dimensions are specified in Space hold, and not SpaceProp
     */
    if ( pp->category() == "space" )
    {
        // update any Space with this property:
        for ( Space * s = simul_.spaces.first(); s; s=s->next() )
        {
            if ( s->prop == pp )
            {
                s->resize(def);
                // allow Simul to update periodic:
                if ( s == simul_.spaces.master() )
                    simul_.spaces.setMaster(s);
            }
        }
    }
}


void Interface::change_simul_property(Glossary& opt)
{
    change_property(&simul_.prop, opt);
}


// in this form, 'name' designates the property name
Property * Interface::execute_change(std::string const& name, Glossary& def, bool strict)
{
    Property * pp = simul_.findProperty(name);
    
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
            e << simul_.properties.all_names(PREF);
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
    PropertyList plist = simul_.findAllProperties(cat);
    
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
    throw InvalidSyntax("unexpected `"+str+"' in `"+StreamFunc::get_line(is, pos)+"'");
}

/**
 Define a placement = ( position, orientation ) from the parameters set in `opt'
 */
Isometry Interface::read_placement(Glossary& opt)
{
    Isometry iso;
    std::string str;
    
    Space const* spc = simul_.spaces.master();
    
    // Space specified as second argument to 'position'
    if ( opt.set(str, "position", 1) )
        spc = simul_.findSpace(str);
    
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

    return iso;
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
 - if placement = `inside` (default), it tries to find a place inside the Space
 - if placement = `anywhere`, the position is returned
 - if placement = `outside`, the object is created only if it is outside the Space
 - if placement = `surface`, the position is projected on the edge of current Space
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
Isometry Interface::find_placement(Glossary& opt, int placement, size_t nb_trials)
{
    size_t ouf = 0;
    std::string str;

    Space const* spc = simul_.spaces.master();
    if ( opt.set(str, "placement", 1) )
        spc = simul_.findSpace(str);
    
    // define a condition:
    bool has_condition = opt.set(str, "placement", 2);
    
    Isometry iso;

    while ( ++ouf < nb_trials )
    {
        // generate a new position:
        iso = read_placement(opt);

        // check all conditions:
        bool condition = true;
        if ( has_condition )
        {
            Evaluator evaluator{{'X', iso.mov.x()}, {'Y', iso.mov.y()}, {'Z', iso.mov.z()},
                                {'R', iso.mov.norm()}, {'P', RNG.preal()}};
            try {
                condition = ( 0 != evaluator.eval(str.c_str()) );
            }
            catch( Exception& e ) {
                e.message(e.message()+" in `"+str+"'");
                throw;
            }
        }
        
        if ( condition )
        {
            if ( !spc || placement == PLACE_ANYWHERE )
                return iso;
            
            if ( placement == PLACE_EDGE )
            {
                iso.mov = spc->project(iso.mov);
                return iso;
            }
            
            if ( spc->inside(iso.mov) )
            {
                if ( placement == PLACE_INSIDE || placement == PLACE_ALL_INSIDE )
                    return iso;
            }
            else
            {
                if ( placement == PLACE_OUTSIDE )
                    return iso;
            }
        }
    }
    
    //Cytosim::warn << "could not fulfill position=`" + opt.value("position") + "'\n";
    throw InvalidParameter("could not fulfill position=`" + opt.value("position") + "'");
    return iso;
}


/**
 This would usually create ONE object of type 'name', placed according to `opt`
 */
void Interface::new_object(ObjectList& objs, std::string const& name, ObjectSet* set, Glossary& opt)
{
    size_t ouf = 0, nb_trials = 1<<14;
    opt.set(nb_trials, "nb_trials");

    objs.clear();
    do {
        
        // create the objects:
        set->newObjects(objs, name, opt);
        
#if ( 0 )
        // check for `nullptr` in list, which should not happen:
        if ( objs.count(nullptr) )
        {
            std::clog << "cytosim found empty slots in newObjects(" << name << ")\n";
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
        
        if ( placement != PLACE_NOT )
        {
            // find a position:
            Isometry iso = find_placement(opt, placement, nb_trials);
            // place object at this position:
            ObjectSet::moveObjects(objs, iso);
            // special case for which we check all vertices:
            if ( placement == PLACE_ALL_INSIDE )
            {
                std::string str;
                Space const* spc = simul_.spaces.master();
                if ( opt.set(str, "placement", 1) )
                    spc = simul_.findSpace(str);
                for ( Object * i : objs )
                {
                    Mecable * mec = Simul::toMecable(i);
                    if ( mec && ! mec->allInside(spc) )
                    {
                        objs.destroy();
                        break;
                    }
                }
            }
        }
        if ( ++ouf > nb_trials )
        {
            Cytosim::log << "could not place `" << name << "' after " << nb_trials << " trials\n";
            break;
        }
    } while ( objs.empty() );
    
    // optionally mark the objects:
    ObjectMark mk = 0;
    if ( opt.set(mk, "mark") )
    {
        for ( Object * i : objs )
            i->mark(mk);
    }
    
    // translation after placement
    Vector vec;
    if ( opt.set(vec, "translation") )
        ObjectSet::translateObjects(objs, vec);

    /* 
     Because the objects in ObjectList are not necessarily all of the same class,
     we call simul_.add() rather than directly set->add()
     */
    simul_.add(objs);
}


/**
 Create `cnt` objects of type 'name', according to specifications.
 */
void Interface::execute_new(std::string const& name, Glossary& opt, size_t cnt)
{
    ObjectSet * set = nullptr;
    Property * pp = simul_.properties.find(name);
    // Allows to make an object without an associated Property
    if ( pp )
        set = simul_.findSet(pp->category());
    else
        set = simul_.findSet(name);
    if ( !set )
    {
        if ( pp )
            throw InvalidSyntax("could not determine the class of `"+name+"'");
        throw InvalidSyntax("undefined class `"+name+"'");
    }
    size_t amount = set->size();
    
    /// allow to set a desired number of objects:
    size_t target = 0;
    if ( opt.set(target, "nb_objects") )
    {
        if ( target < amount )
        {
            ObjectList objs = set->collect();
            simul_.erase(objs, amount-target);
            return;
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
        
        ObjectList objs(4, 4);
        for ( size_t n = 0; n < cnt; ++n )
        {
            opt.define("position", A + n * dAB);
            new_object(objs, name, set, opt);
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

        ObjectList objs(4, 4);
        for ( size_t n = 0; n < cnt; ++n )
            new_object(objs, name, set, opt);
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

    VLOG("+NEW `" << name << "' made " << set->size()-amount << " objects (total " << simul_.nbObjects() << ")");
}


//------------------------------------------------------------------------------
/**
 Creates `cnt` objects of class `name`.
 The objects are distributed uniformly within the current Space, also with
 random orientations.
 
 This is meant to replace execute_new(name, opt, cnt), when no option was given
 to the command.
 */
void Interface::execute_new(std::string const& name, size_t cnt)
{
    Property * pp = simul_.properties.find_or_die(name);
    ObjectSet * set = simul_.findSet(pp->category());
    if ( !set )
        throw InvalidSyntax("could not determine the class of `"+name+"'");

    Space const* spc = simul_.spaces.master();

    Glossary opt;

    ObjectList objs(4, 4);
    for ( size_t n = 0; n < cnt; ++n )
    {
        set->newObjects(objs, name, opt);
        
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
        
        /* Call simul_.add() rather than directly set->add(), because the objects
         in ObjectList are not necessarily all of the same class */
        simul_.add(objs);
        objs.clear();
    }
    
    VLOG("-NEW " << cnt << "`" << name << "' objects");
    //hold();
}

//------------------------------------------------------------------------------
#pragma mark -

/// holds a set of criteria used to select Objects
class Filter
{
public:

    ObjectMark      mrk;
    unsigned        st1;
    unsigned        st2;
    Space    const* ins;
    Space    const* ous;
    Property const* prp;

    /// initialize
    Filter()
    {
        mrk = 0;
        st1 = ~0U;
        st2 = ~0U;
        prp = nullptr;
        ins = nullptr;
        ous = nullptr;
    }
    
    void set(Simul& sim, Property* pp, Glossary& opt)
    {
        prp = pp;
        
        std::string str;
        if ( opt.set(str, "position") )
        {
            Space const* spc = nullptr;
            std::string spn;
            if ( opt.set(spn, "position", 1) )
                spc = sim.findSpace(spn);
            else
                spc = sim.spaces.master();
            if ( !spc )
                throw InvalidSyntax("unknown Space `"+spn+"'");
            
            if ( str == "inside" )
                ins = spc;
            else if ( str == "outside" )
                ous = spc;
            else
                throw InvalidSyntax("unknown specification `"+str+"'");
        }
        
        opt.set(mrk, "mark");
        opt.set(st1, "state1", "stateP") || opt.set(st1, "state");
        opt.set(st2, "state2", "stateM") || opt.set(st1, "state", 1);
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
    Property * pp = simul_.properties.find(name);
    ObjectSet * set = nullptr;
    if ( pp )
        set = simul_.findSet(pp->category());
    else
        set = simul_.findSet(name);
    if ( !set )
    {
        if ( name == "objects" )
        {
            simul_.erase_all(0);
            return;
        }
        throw InvalidSyntax("could not determine the class of `"+name+"'");
    }
    
    Filter filter;
    filter.set(simul_, pp, opt);
    ObjectList objs = set->collect(pass_filter, &filter);
    
    if ( objs.size() > 0 )
        simul_.erase(objs, cnt);
    else
        Cytosim::warn << "found no `" << name << "' to delete\n";
}

/**
 This can only mark one class of objects
 */
void Interface::execute_mark(std::string const& name, Glossary& opt, size_t cnt)
{
    Property * pp = simul_.properties.find(name);
    ObjectSet * set = nullptr;
    if ( pp )
        set = simul_.findSet(pp->category());
    else
        set = simul_.findSet(name);
    if ( !set )
        throw InvalidSyntax("could not determine the class of `"+name+"'");

    ObjectMark mrk;
    if ( ! opt.set(mrk, "mark") )
        throw InvalidParameter("mark must be specified for command `mark'");
    opt.clear("mark");
    
    Filter filter;
    filter.set(simul_, pp, opt);
    ObjectList objs = set->collect(pass_filter, &filter);
    
    // optionally limit the list to a random subset
    if ( cnt < objs.size() )
    {
        objs.shuffle();
        objs.truncate(cnt);
    }
    
    simul_.mark(objs, mrk);
}


void Interface::execute_cut(std::string const& name, Glossary& opt)
{
    Vector n(1,0,0);
    real a = 0;
    
    opt.set(n, "plane");
    opt.set(a, "plane", 1);
    
    state_t stateP = STATE_RED, stateM = STATE_GREEN;
    opt.set(stateP, "new_end_state");
    opt.set(stateM, "new_end_state", 1);
    
    ObjectList objs;

    if ( name == "all" )
    {
        objs = simul_.fibers.collect();
    }
    else
    {
        Property * pp = simul_.properties.find_or_die(name);
        if ( pp->category() != "fiber" )
            throw InvalidSyntax("only `cut fiber' is supported");
        
        Filter filter;
        filter.set(simul_, pp, opt);
        objs = simul_.fibers.collect(pass_filter, &filter);
    }
    
    VLOG("-CUT PLANE (" << n << ").x = " << -a);
    simul_.fibers.planarCut(objs, n, a, stateP, stateM);
}


void Interface::execute_connect(std::string const& name, Glossary& opt)
{
    ObjectList objs;

    if ( name == "couple" )
    {
        simul_.couples.bindToIntersections();
    }
    else
    {
        Property * pp = simul_.properties.find_or_die(name);
        if ( pp->category() != "couple" )
            throw InvalidSyntax("only `bind couple' or `bind couple_class' is supported");
        simul_.couples.bindToIntersections(static_cast<CoupleProp*>(pp));
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
inline void Interface::do_steps(size_t& sss, size_t cnt)
{
    while ( sss < cnt )
    {
        hold();
        //fprintf(stderr, "> step %6zu\n", sss);
        (simul_.*FUNC)();
        simul_.step();
        ++sss;
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
void Interface::execute_run(size_t nb_steps, Glossary& opt, bool do_write)
{
    int  frames = 0;
    int  solve = 1;
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
        event->reload(simul_.time());
        simul_.events.add(event);
    }
#endif
    opt.set(solve, "solve", {{"off",0}, {"on",1}, {"auto",2}, {"force", 3},
        {"horizontal",4}, {"flux",5}, {"half",7}, {"separate", 8} });
    opt.set(prune,  "prune");
    opt.set(binary, "binary");
    opt.set(frames, "nb_frames");
    
    do_write &= ( frames > 0 );

    real   delta = (real)nb_steps;
    size_t check = nb_steps;
    
    simul_.prepare();

    if ( do_write )
    {
        simul_.writeProperties(prune);
        // write initial state:
        if ( simul_.prop.clear_trajectory )
        {
            simul_.prop.clear_trajectory = false;
            if ( frames > 1 )
                simul_.writeObjects(simul_.prop.system_file, false, binary);
            else
                std::remove(simul_.prop.system_file.c_str());
        }
        delta = real(nb_steps) / real(frames);
        check = size_t(delta);
    }
    
    VLOG("+RUN START " << nb_steps);
    
    size_t frm = 0;
    size_t sss = 0;
    do {
        switch ( solve )
        {
            case 0: do_steps<&Simul::solve_not>(sss, check); break;
            case 1: do_steps<&Simul::solve>(sss, check); break;
            case 2: do_steps<&Simul::solve_auto>(sss, check); break;
            case 3: do_steps<&Simul::solve_force>(sss, check); break;
            case 4: do_steps<&Simul::solve_onlyX>(sss, check); break;
            case 5: do_steps<&Simul::solve_flux>(sss, check); break;
            case 7: do_steps<&Simul::solve_half>(sss, check); break;
            case 8: do_steps<&Simul::solve_separate>(sss, check); break;
        }
        ++frm;
        check = size_t(delta*(frm+1));

        if ( do_write )
        {
            simul_.relax();
            simul_.writeObjects(simul_.prop.system_file, true, binary);
            reportCPUtime(simul_.time());
            simul_.sMeca.doNotify = 2;  // to print convergence parameters
            simul_.unrelax();
        }
    } while ( sss < nb_steps );
    
#if BACKWARD_COMPATIBILITY < 50
    if ( event )
        simul_.events.erase(event);
#endif
    
    simul_.relax();
    if ( frames )
    {
        simul_.writeProperties(prune);
        if ( frames < 0 )
            simul_.writeObjects(simul_.prop.system_file, true, binary);
    }
    
    VLOG("+RUN END");
}


/**
 Perform simulation steps, without any option, by alternating step() and solve()
*/
void Interface::execute_run(size_t nb_steps)
{
    VLOG("-RUN START " << nb_steps);
    simul_.prepare();
    
    for ( size_t sss = 0; sss < nb_steps; ++sss )
    {
        hold();
        //fprintf(stderr, "> step %6zu\n", sss);
        simul_.solve();
        simul_.step();
    }
    
    simul_.relax();
    VLOG("-RUN END");
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
        subset = simul_.findSet(what);
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
            real t = simul_.prop.time;
            simul_.reloadObjects(in, 0, subset);
            simul_.prop.time = t;
        }
        else
            simul_.reloadObjects(in, 1, subset);
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
         simul_.mark(objs, mrk);
    }
#endif
    
    // set time
    real t;
    if ( opt.set(t, "time") )
        simul_.prop.time = t;
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
            simul_.writeObjects(name, append, binary);
        }
        else
        {
            Outputter out(stdout, false);
            simul_.writeObjects(out);
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
        simul_.writeProperties(out, prune);
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

    simul_.mono_report(out, what, opt, ver);
}


void Interface::execute_call(std::string& str, Glossary& opt)
{
    if ( str == "equilibrate" )
        simul_.couples.equilibrate();
    else if ( str == "connect" )
        simul_.couples.bindToIntersections();
    else if ( str == "custom0" )
        simul_.custom0(opt);
    else if ( str == "custom1" )
        simul_.custom1(opt);
    else if ( str == "custom2" )
        simul_.custom2(opt);
    else if ( str == "custom3" )
        simul_.custom3(opt);
    else if ( str == "custom4" )
        simul_.custom4(opt);
    else if ( str == "custom5" )
        simul_.custom5(opt);
    else if ( str == "custom6" )
        simul_.custom6(opt);
    else if ( str == "custom7" )
        simul_.custom7(opt);
    else if ( str == "custom8" )
        simul_.custom8(opt);
    else if ( str == "custom9" )
        simul_.custom9(opt);
    else
        throw InvalidSyntax("called unknown command");
}

