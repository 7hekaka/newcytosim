// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.
#ifndef SIMUL_H
#define SIMUL_H

#include "assert_macro.h"
#include <iostream>
#include <string>
#include <stack>
#include <map>

#include "simul_part.h"
#include "simul_prop.h"
#include "single_set.h"
#include "couple_set.h"
#include "fiber_set.h"
#include "fiber_grid.h"
#include "bead_set.h"
#include "solid_set.h"
#include "sphere_set.h"
#include "organizer_set.h"
#include "field_set.h"
#include "space_set.h"
#include "event_set.h"
#include "property_list.h"
#include "meca.h"

class Meca1D;
class Parser;
class Aster;


/// Simulator class containing all Objects
class Simul
{
public:
    
    /// Meca used to set and integrate the equations of motion of Mecables
    mutable Meca sMeca;
    
    /// grid used for attachment of Hand to Fiber
    mutable FiberGrid fiberGrid;

    /// Meca used to solve the system with option 'solve=horizontal'
    Meca1D * pMeca1D;
    
    /// member used to indicate that the simulation state was renewed
    int fresh_;

    /// text for miscellaneous operations
    mutable std::string text_;
    
private:
    
    /// signals that Simul is ready to perform a Monte-Carlo step
    mutable int primed_;

    /// current parser object
    Parser * parser_;

    /// preconditionning method used/determined by `solve_auto()`
    unsigned autoPrecond;
    
    /// counter used by `solve_auto()`
    unsigned autoCounter;
    
    /// record of CPU time for `solve_auto()`
    float    autoCPU[8];
    
    /// record of iteration count for `solve_auto()`
    unsigned autoCNT[8];

#if POOL_UNATTACHED > 1
    /// counter for occasional Hand's attachment
    unsigned doAttachCounter;
#endif

    /// a copy of the properties as they were stored to file
    mutable std::string properties_saved;

public:

    /// Global cytosim parameters
    mutable SimulProp prop;
    
    /// list of all Object Property, except SimulProp
    PropertyList properties;
    
    /// list of Space in the Simulation
    SpaceSet spaces;
    
    /// list of Field in the Simulation
    FieldSet fields;
    
    /// list of Fiber in the Simulation
    FiberSet fibers;
    
    /// list of Sphere in the Simulation
    SphereSet spheres;
    
    /// list of Bead in the Simulation
    BeadSet beads;
    
    /// list of Solid in the Simulation
    SolidSet solids;
    
    /// list of Single in the Simulation
    SingleSet singles;
    
    /// list of Couple in the Simulation
    CoupleSet couples;
    
    /// list of Organizer in the Simulation
    OrganizerSet organizers;
    
    /// list of Events in the Simulation
    EventSet events;

    //--------------------------------------------------------------------------
    
    /// constructor
    Simul();
    
    /// destructor
    virtual ~Simul();
        
    //----------------------------- POPULATING ---------------------------------
    
    /// add Object to Simulation
    void add(Object *);

    /// add all Objects from given list
    void add(ObjectList const&);

    /// remove Object from Simulation
    void remove(Object *);

    /// remove all Objects from given list
    void remove(ObjectList const&);
    
    /// remove and delete object
    void eraseObject(Object *);
    
    /// remove and delete all objects from given list
    void eraseObjects(ObjectList const&);

    /// reset simulation world (clear all sub-lists and variables)
    void eraseObjects();

    /// reset simulation world (clear all sub-lists and variables)
    void eraseProperties();

    /// total number of objects in the Simulation
    size_t nbObjects() const;

    /// maximum Hand::binding_range
    real maxBindingRange() const;
    
    //----------------------------- SIMULATING ---------------------------------
   
    /// perform basic initialization; register callbacks
    void initCytosim();
    
    /// ready the engine for a subsequent call to `step` and `solve`
    void prepare();
    
    /// perform one Monte-Carlo step, corresponding to the time step
    void steps();

    /// time step (shortcut to `SimulProp::time_step`)
    double time_step() const { return prop.time_step; }

    /// time in the simulated world (shortcut to `SimulProp::time`)
    double time() const { return prop.time; }

    /// change time in the simulated world
    void set_time(double t) { prop.time = t; }

    /// true if `SimulProp::time < SimulProp::stop_time`
    bool incomplete() const { return prop.time < prop.stop_time; }

    /// ask to halt the current series of 'step' once time exceeds `t`
    void stop_at(double t) const { prop.stop_time = std::min(t, prop.end_time); }
    
    /// ask to stop the `run` at the next recorded frame, once time exceeds 't'
    void end_at(double t) const { prop.end_time = std::min(t, prop.end_time); }
    
    /// ask to stop the `run` after `t` seconds beyond current time
    void end_now(double t = 0.0) const { t += time(); end_at(t); if ( t < prop.stop_time ) stop_at(t); }

    /// true if `time >= end_time`
    bool should_end() const { return prop.time >= prop.end_time; }

    /// reset `end_time` such as to never terminate
    void end_never() const { prop.end_time = INFINITY; }

    /// this is called after a sequence of `step()` have been done
    void relax();
    
    /// this is called to signal that engine is ready to proceed
    void unrelax() { primed_ = 2; }

    /// true if engine is ready to go (between `prepare()` and `relax()`)
    int primed() const { return primed_; }

    
    /// call setInteractions(Meca) for all objects (this is called before `solve`
    void setAllInteractions(Meca&) const;

    /// display Meca's links
    void drawLinks() const;
    
    /// bring all objects to centered image using periodic boundary conditions
    void foldPositions() const;

    /// simulate the mechanics of the system and move Mecables accordingly
    void solve_meca();
    
    /// prepares system for solving
    void prepare_meca() { sMeca.getReady(*this); };

    /// same as `solve_meca` but skipping the `prepare_meca` part.
    void prepared_solve();

    /// calculate forces given the current positions
    void solve_force();

    /// solve mechanical system and calculate forces but do not apply movements
    void solve_half();

    /// replace coordinates of Mecables by the ones calculated in `solve`
    void apply() const { sMeca.apply(); }
    
    /// like 'solve' but automatically selects the fastest preconditionning method
    void solve_auto();

    /// do nothing
    void solve_not() {};

    /// calculate the motion of objects, but only in the X-direction
    void solve_uniaxial();
    
    /// calculate Forces and Lagrange multipliers on the Mecables
    void computeForces() const;
    
    /// calculate clusters based on Couples' connectivity and order Fibers accordingly
    ObjectFlag orderClustersCouple(Object**, ObjectFlag);
    
    /// this is used for development
    void addExperimentalInteractions(Meca&) const;
    
    /// set FiberGrid and StericGrid over the given space
    static void setFiberGrid(FiberGrid&, Space const*, real& grid_step);
    
    /// a Map to be displayed
    Map<DIM> const& visibleMap() const { return sMeca.locusGrid.map(); }
    
    /// return estimate for the cell size of the PointGrid used for steric interactions
    real minimumStericRange() const;

private:
    
    /// return estimate of the cell size of the Grid used for Hand attachment
    real estimateFiberGridStep() const;
    
    //----------------------------- PARSING ------------------------------------
    
    /// return the ObjectSet corresponding to this Tag in the simulation (used for IO)
    ObjectSet * findSetT(const ObjectTag);

public:
    
    /// return the ObjectSet corresponding to a class
    ObjectSet * findSet(const std::string& cat);
    
    /// convert Object to Mecable* if the conversion seems valid; returns 0 otherwise
    static Mecable * toMecable(Object *);

    /// find a Mecable from a string specifying name and inventory number (e.g. 'fiber1')
    Mecable * pickMecable(const std::string&) const;
    
    /// find a Movable from a string specifying name and inventory number (e.g. 'fiber1')
    Object * pickMovable(const std::string&) const;

    /// find a Solid Prop
    SolidProp* findSolidProp(std::string s) const { return static_cast<SolidProp*>(findProperty("solid", s)); }

    /// find a Solid by name
    Solid * pickSolid(std::string s) { return Solid::toSolid(solids.pickObject("solid", s)); }
    
    /// find a Fiber by name
    Fiber * pickFiber(std::string s) { return Fiber::toFiber(fibers.pickObject("fiber", s)); }
    
    /// find a Sphere by name
    Sphere * pickSphere(std::string s) { return Sphere::toSphere(spheres.pickObject("sphere", s)); }
    
    /// first Space with this Property
    Space * pickSpace(const Property * p) const { return Space::toSpace(spaces.pickObject(p)); }

    /// return first Space with given name, or return nullptr
    Space const* findSpace(std::string) const;
    
    /// first Field with this Property
    Field * pickField(const Property * p) const;
    
    //------------------------- CODE EVALUATOR ---------------------------------

    /// set Parser to be used in `perform()` below
    void parser(Parser * p) { parser_ = p; }
    
    /// return current Parser
    Parser * parser() const { return parser_; }

    /// Parse a text containing cytosim commands
    void perform(std::string const&);

    //---------------------------- PROPERTIES ----------------------------------

    /// change the name of the simulation
    void rename(std::string const&);
    
    /// read an Object reference and return the corresponding Object (`tag` is set)
    Object * readReference(Inputter&, ObjectTag& tag);

    /// read a Fiber reference and return the corresponding Object (`tag` is set)
    Fiber * readFiberReference(Inputter&, ObjectTag& tag, ObjectID&);
    
    
    /// check if `name` corresponds to a property class, but excluding 'simul'
    bool isCategory(const std::string& name) const;
    
    /// return existing property of given class and name, or return zero
    Property * findProperty(const std::string&, const std::string&) const;
    
    /// return existing property of given name, or return zero
    Property * findProperty(const std::string&) const;

    /// return all existing properties of requested class
    PropertyList findAllProperties(const std::string&) const;
    
    /// return Property in the requested type, or throw an exception
    template < typename T >
    T * findProperty(std::string const& cat, PropertyID id) const
    {
        Property * p = properties.find(cat, id);
        if ( !p )
        {
#if ( 0 )
            p = properties.find(cat, 1);
            if ( p )
                std::cerr << "Substituting " << cat << " ID 1 for ID " << id << "\n";
            else
#endif
                throw InvalidIO("undefined `"+cat+"' class with ID "+std::to_string(id));
        }
        return static_cast<T*>(p);
    }
    
    /// return Property in the requested type, or throw an exception
    template < typename T >
    T * findProperty(std::string const& cat, std::string const& nom) const
    {
        Property * p = properties.find(cat, nom);
        if ( !p )
            throw InvalidIO("could not find "+cat+" class `"+nom+"'");
        return static_cast<T*>(p);
    }

    /// create a new property
    Property* makeProperty(const std::string&, const std::string&, Glossary&);
    
    /// export all Properties to speficied file
    void writeProperties(std::ostream&, bool prune) const;
    
    /// export all Properties using default file name
    void writeProperties(bool prune) const;
    
    /// load the properties contained in given file
    void loadProperties(char const* filename, int verbose);

    /// load the properties contained in the standard output property file
    void loadProperties();
   
    /// read and set parameter for some object; syntax is `CLASS:PARAMETER=VALUE`
    bool readParameter(const char*) const;

    //---------------------------- LOAD OBJECTS --------------------------------
    
    /// class for reading trajectory file
    class ImportLock;

    /// default name for output trajectory file
    static const char TRAJECTORY[];

    /// current file format (check history in `simul_file.cc`)
    static constexpr unsigned currentFormatID = 60U;
    
    /// read objects from file, and add them to the simulation state
    int readMetaData(Inputter&, std::string& section, ObjectSet*& objset, ObjectSet* subset);

    /// read objects from file, and add them to the simulation state
    int readObjects(Inputter&, ObjectSet* subset);

    /// load sim-world from the named file
    int loadObjects(char const* filename);
    
    /// import objects from file, and delete objects that were not referenced in the file
    int reloadObjects(Inputter&, bool prune = 1, ObjectSet* subset = nullptr);

    /// write sim-world to specified file
    void writeObjects(Outputter&) const;
    
    /// write sim-world in binary or text mode, appending to existing file or creating new file
    void writeObjects(std::string const& filename, bool append, int binary) const;
    
    //----------------------------- REPORTING ----------------------------------
    
    /// calls report_one
    void poly_report(std::ostream&, std::string, Glossary&, int) const;
    
    /// calls report
    void mono_report(std::ostream&, std::string const&, Glossary&, int) const;

    /// calls report_one
    void report_one(std::ostream&, std::ostream&, std::string const&, Glossary&) const;
    
    /// call one of the report function
    void report_one(std::ostream&, std::ostream&, std::string const&, Property const*, std::string const&, Glossary&) const;

    
    /// print time
    void reportTime(std::ostream&) const;
    
    /// give a short inventory of the simulation state, obtained from ObjectSet::report()
    void reportInventory(std::ostream&) const;
    
    /// give a summary of the Simul
    void reportSimul(std::ostream&, std::ostream&) const;
    
    /// print the length and the points of each fiber
    void reportFiber(std::ostream&, Fiber const*) const;

    /// print the length and the points of each fiber
    void reportFibers(std::ostream&, std::ostream&, Property const*) const;
    
    /// print the length and the points of each fiber, sorted from longest to shortest
    void reportFibersSorted(std::ostream&, std::ostream&, Property const*);

    /// print the coordinates of the vertices of each fiber
    void reportFiberPoints(std::ostream&, std::ostream&, Property const*) const;
    
    /// print the coordinates of the vertices of each fiber
    void reportFiberDisplacement(std::ostream&, std::ostream&, Property const*) const;
    
    /// print the coordinates of the vertices of each fiber
    void reportFiberDirections(std::ostream&, std::ostream&, Property const*) const;

    /// print the positions and the states for one type of end of all fibers
    void reportFiberEnds(std::ostream&, std::ostream&, FiberEnd, Property const*) const;
    
    /// print number of fibers in each state of specified end
    void reportFiberEndState(std::ostream&, std::ostream&, FiberEnd, Property const*) const;

    /// print the mean and standard deviation of vertices for each class of fiber
    void reportFiberMoments(std::ostream&, std::ostream&) const;

    /// print average age and standard deviation for each class of fiber
    void reportFiberAge(std::ostream&, std::ostream&) const;
    
    /// print average length and standard deviation for each class of fiber
    void reportFiberMarkedLengths(std::ostream&, std::ostream&, Property const*) const;

    /// print average length and standard deviation for each class of fiber
    void reportFiberLengths(std::ostream&, std::ostream&, Property const*) const;
    
    /// print length distribution for each class of fiber
    void reportFiberLengthHistogram(std::ostream&, std::ostream&, Glossary&) const;
    
    /// print coordinates of speckles that follow a frozen random sampling
    void reportFiberSpeckles(std::ostream&, std::ostream&, Glossary&) const;
    
    /// print coordinates of points randomly and freshly distributed
    void reportFiberSamples(std::ostream&, std::ostream&, Glossary&) const;
    
    /// print the coordinates and forces on the vertices of each fiber
    void reportFiberForces(std::ostream&, std::ostream&) const;

    /// print Fiber tensions along certain planes defined in Glossary
    void reportFiberTension(std::ostream&, std::ostream&, Glossary&) const;
    
    
    /// print number of kinks in each class of Fiber
    void reportFiberSegments(std::ostream&, std::ostream&) const;

    /// document bending energy
    void reportFiberBendingEnergy(std::ostream&, std::ostream&) const;
    
    /// document end-to-end distance in each class of Fiber
    void reportFiberExtension(std::ostream&, std::ostream&) const;
    
    /// document nematic order (alignement) of Fiber
    void reportFiberNematic(std::ostream&, std::ostream&, Glossary& opt) const;

    /// document nematic order (alignement) of Fiber
    void reportFiberNematic(std::ostream&, std::ostream&, FiberProp const*, Space const*) const;
    
    /// print component of forces experienced by Fibers due to confinement
    void reportFiberConfineForce(std::ostream&, std::ostream&) const;

    /// print radial component of forces experienced by Fibers due to confinement
    real reportFiberConfinement(std::ostream&, std::ostream&) const;

    /// print summary of Fiber's lattice quantities
    void reportFiberLattice(std::ostream&, std::ostream&, Property const*) const;
    
    /// print values of Fiber's densities
    void reportFiberDensityTotal(std::ostream&, std::ostream&, bool density, Property const*) const;
    
    /// print summary of Fiber's densities
    void reportFiberDensity(std::ostream&, std::ostream&, bool density, Property const*) const;
    
    
    /// print position of hands bound to fibers
    void reportFiberHands(std::ostream&, std::ostream&) const;
    
    /// print position of bound hands that are associated with stiffness
    void reportFiberLinks(std::ostream&, std::ostream&) const;

    /// print interection abscissa between fibers
    void reportFiberConnectors(std::ostream&, std::ostream&, Glossary&) const;

    /// print interection abscissa between fibers
    void reportNetworkBridges(std::ostream&, std::ostream&, Glossary&) const;
 
    /// print network surface area
    void reportNetworkSize(std::ostream&, std::ostream&) const;

    /// print positions of interection between two fibers
    void reportFiberIntersections(std::ostream&, std::ostream&, Glossary&) const;


    /// print Organizer positions
    void reportOrganizer(std::ostream&, std::ostream&) const;

    /// print Aster positions
    void reportAster(std::ostream&, std::ostream&) const;
    
    /// print Bead positions 
    void reportBeadSingles(std::ostream&, std::ostream&) const;

    /// print Bead positions
    void reportBeadPosition(std::ostream&, std::ostream&, Property const*) const;

    /// print Solid positions 
    void reportSolidPosition(std::ostream&, std::ostream&, Property const*) const;
    
    /// print Solid positions
    void reportSolidOrientation(std::ostream&, std::ostream&, Property const*) const;

    /// print Solid's anchored Hands
    void reportSolidHands(std::ostream&, std::ostream&, Property const*) const;

    /// print state of Couples 
    void reportCouple(std::ostream&, std::ostream&, Property const*) const;
 
    /// print state of Couples
    void reportCoupleList(std::ostream&, std::ostream&, Property const*) const;

    /// print state of Couples
    void reportCoupleAnatomy(std::ostream&, std::ostream&, Property const*) const;
    
    /// print position of Couples of a certain kind
    void reportCoupleState(std::ostream&, std::ostream&, Property const*) const;
    
    /// print position of active Couples of a certain kind
    void reportCoupleActive(std::ostream&, std::ostream&, Property const*) const;
    
    /// print position and forces of doubly-attached Couples
    void reportCoupleLink(std::ostream&, std::ostream&, Property const*) const;
    
    /// print configurations of doubly-attached Couples
    void reportCoupleConfiguration(std::ostream&, std::ostream&, Property const*, Glossary&) const;

    /// print agregate properties of Couple force
    void reportCoupleForce(std::ostream&, std::ostream&, Property const*) const;
    
    /// print histogram of Couples force
    void reportCoupleForceHistogram(std::ostream&, std::ostream&, Glossary&) const;

    /// print state of Singles
    void reportSingle(std::ostream&, std::ostream&, Property const*) const;
    
    /// print position of Singles of a certain kind
    void reportSingleState(std::ostream&, std::ostream&, Property const*) const;

    /// print position of Singles
    void reportSinglePosition(std::ostream&, std::ostream&, Property const*) const;
   
    /// print force of attached Singles
    void reportSingleLink(std::ostream&, std::ostream&, Property const*) const;
    
    /// print agregate properties of Singles force
    void reportSingleForce(std::ostream&, std::ostream&, Property const*) const;

    /// print state of Couples 
    void reportSpherePosition(std::ostream&, std::ostream&, Property const*) const;

    /// print something about Spaces
    void reportSpace(std::ostream&, std::ostream&) const;
  
    /// print force on Spaces
    void reportSpaceForce(std::ostream&, std::ostream&) const;

    /// print something about Fields
    void reportField(std::ostream&, std::ostream&) const;
    
    //------------------------- CLUSTER ANALYSIS -------------------------------
    
    /// set Mecable's flag() to unique values, return highest value attributed + 1
    ObjectFlag setUniqueFlags() const;

    /// call `flag(f)` for all Mecable
    void setFlags(ObjectFlag f) const;
    
    /// replace all occurence of flag `f` by `g`
    void changeFlags(ObjectFlag f, ObjectFlag g) const;

    /// calculate clusters of Mecable derived from all interactions
    void flagClustersMeca() const;
    
    /// flag fibers according to connectivity defined by Couple
    void flagClustersCouples() const;
    
    /// flag fibers according to connectivity defined by Solids
    void flagClustersSingles() const;
    
    /// flag fibers according to connectivity defined by Couple of given type
    void flagClustersCouples(Property const*) const;

    /// analyse the network connectivity to identify isolated sub-networks
    void flagClusters(bool cop, bool sol, bool mec) const;
    
    /// change flag for fibers belonging to largest cluster
    void flagLargestCluster(ObjectFlag) const;
    
    /// print size of clusters defined by connections with Couples
    void reportClusters(std::ostream&, std::ostream&, Glossary&) const;

    //------------------------- PROJECT SPECIFIC -------------------------------

    /// flag the fibers that appear to constitute a ring
    size_t flagRing() const;
    
    /// extract radius and length of ring
    void analyzeRing(ObjectFlag, real& length, real& radius) const;
    
    /// estimates if Fibers form a connected ring around the Z-axis
    void reportRing(std::ostream&, std::ostream&) const;

    /// custom report for Platelet project
    void reportPlatelet(std::ostream&, std::ostream&) const;
    
    /// print Aster & Spindle indices
    void reportSpindleIndices(std::ostream&, std::ostream&) const;

    /// report position of 'condensate' beads on left and right sides (06.2023)
    void reportSpindleLength(std::ostream&, std::ostream&, Glossary&) const;

    /// print number of Fibers pointing left and right that intersect plane YZ at different X positions
    void reportSpindleProfile(std::ostream&, std::ostream&, Glossary&) const;

    /// print some coefficients calculated from the distribution of fibers
    void reportSpindleFitness(std::ostream&, std::ostream&, Glossary&) const;

    /// report position of minus ends for fiber per marks on left and right sides (08.2023)
    void reportMarkedFiberEnds(std::ostream&, std::ostream&, Glossary&) const;

    /// a special print for Romain Gibeaux
    void reportAshbya(std::ostream&, std::ostream&) const;
    
    /// analysis of MT collisions in the plant cortex
    void reportFiberCollision(std::ostream&, Fiber*, Fiber const*, int) const;

    /// analysis of MT collisions in the plant cortex
    void reportFiberCollision(std::ostream&, std::ostream&, Property const*, Glossary&) const;

    /// print something
    void reportSomething(std::ostream&, std::ostream&) const;

    //------------------------------ CUSTOM ------------------------------------

    /// custom function
    void custom0(Glossary&);
    /// custom function
    void custom1(Glossary&);
    /// custom function
    void custom2(Glossary&);
    /// custom function
    void custom3(Glossary&);
    /// custom function
    void custom4(Glossary&);
    /// custom function
    void custom5(Glossary&);
    /// custom function
    void custom6(Glossary&);
    /// custom function
    void custom7(Glossary&);
    /// custom function
    void custom8(Glossary&);
    /// custom function
    void custom9(Glossary&);
};

#endif

