// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SIMUL_H
#define SIMUL_H

#include "assert_macro.h"
#include <iostream>
#include <string>
#include <stack>
#include <map>

#include "simul_prop.h"
#include "single_set.h"
#include "couple_set.h"
#include "fiber_set.h"
#include "fiber_grid.h"
#include "bead_set.h"
#include "solid_set.h"
#include "sphere_set.h"
#include "organizer_set.h"
#include "tubule_set.h"
#include "field_set.h"
#include "space_set.h"
#include "event_set.h"
#include "property_list.h"
#include "meca.h"

class Meca1D;
class Parser;


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
    
#if POOL_UNATTACHED > 1
    /// counter for occasional Hand's attachment
    unsigned doAttachCounter;
#endif
    
    /// flag to end the 'run' as soon as possible
    unsigned abortRun;

private:
    
    /// current parser object
    Parser * parser_;

    /// signals that Simul is ready to perform a Monte-Carlo step
    int      primed_;

    /// preconditionning method used/determined by `solve_auto()`
    unsigned autoPrecond;
    
    /// counter used by `solve_auto()`
    size_t   autoCounter;
    
    /// record of CPU time for `solve_auto()`
    float    autoCPU[8];
    
    /// record of iteration count for `solve_auto()`
    size_t   autoCNT[8];
    
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

    /// list of Tubules
    TubuleSet tubules;
    
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
    void erase(Object *);
    
    /// remove and delete all objects from given list
    void erase(ObjectList const&);
    
    /// remove and delete 'cnt' objects from given list
    void erase(ObjectList&, size_t);

    /// mark objects from given list
    static void mark(ObjectList const&, ObjectMark);

    /// reset simulation world (clear all sub-lists and variables)
    void eraseObjects(bool erase_properties);

    /// total number of objects in the Simulation
    size_t nbObjects() const;

    //----------------------------- SIMULATING ---------------------------------
   
    /// perform basic initialization; register callbacks
    void initialize(Glossary&);
    
    /// ready the engine for a subsequent call to `step()` and `solve()`
    void prepare();
    
    /// perform one Monte-Carlo step, corresponding to the time step
    void step();

    /// time in the simulated world (shortcut to `SimulProp::time`)
    double time() const { return prop.time; }
    
    /// shortcut to prop->time_step;
    real time_step() const { return prop.time_step; }

    /// this is called after a sequence of `step()` have been done
    void relax();
    
    /// this is called to signal that engine is ready to proceed
    void unrelax() { primed_ = 2; }

    /// true if engine is ready to go (between `prepare()` and `relax()`)
    int primed() const { return primed_; }

    
    /// call setInteractions(Meca) for all objects (this is called before `solve()`
    void setAllInteractions(Meca&) const;

    /// display Meca's links
    void drawLinks() const;
    
    /// bring all objects to centered image using periodic boundary conditions
    void foldPositions() const;

    /// simulate the mechanics of the system and move Mecables accordingly
    void solve();
    
    /// calculate forces given the current positions
    void solve_force();

    /// solve mechanical system and calculate forces but do not apply movements
    void solve_half();

    /// replace coordinates of Mecables by the ones calculated in solve()
    void apply() const { sMeca.apply(); }
    
    /// like 'solve' but automatically selects the fastest preconditionning method
    void solve_auto();
    
    /// simulate the mechanics of the system and move Mecables accordingly
    void solve_separate();

    /// do nothing
    void solve_not() {};

    /// calculate the motion of objects, but only in the X-direction
    void solve_onlyX();
    
    /// move every Fiber backward by `shift` (this is an extremely crude model)
    void solve_flux();
    
    /// calculate Forces and Lagrange multipliers on the Mecables
    void computeForces() const;
    
    /// calculate clusters based on Couples' connectivity and order Fibers accordingly
    ObjectFlag orderClustersCouple(Object**, ObjectFlag);
    
    /// this is used for development
    void addExperimentalInteractions(Meca&) const;
    
    /// set FiberGrid and StericGrid over the given space
    void setFiberGrid(Space const*, real& grid_step) const;
    
    /// a Map to be displayed
    Map<DIM> const& visibleMap() const { return sMeca.locusGrid.map(); }
    
    /// return estimate for the cell size of the PointGrid used for steric interactions
    real estimateStericRange() const;

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
    
    /// find a Mecable from a string specifying name and inventory number (e.g. 'fiber1')
    Object * pickMovable(const std::string&) const;

    /// find a Solid by name
    Solid * pickSolid(std::string s) { return Solid::toSolid(solids.pickObject("solid", s)); }
    
    /// find a Fiber by name
    Fiber * pickFiber(std::string s) { return Fiber::toFiber(fibers.pickObject("fiber", s)); }
    
    /// find a Sphere by name
    Sphere * pickSphere(std::string s) { return Sphere::toSphere(spheres.pickObject("sphere", s)); }
    
    /// first Space with this Property
    Space * pickSpace(const Property * p) const { return static_cast<Space*>(spaces.pickObject(p)); }

    /// return first Space with given name, or return nullptr
    Space const* findSpace(std::string const&) const;
    
    /// first Field with this Property
    Field * pickField(const Property * p) const;
    
    //---------------------------- PROPERTIES ----------------------------------

    /// set Parser used in' evaluate'
    void parser(Parser* p) { parser_ = p; }
    
    /// Parse a text containing cytosim commands
    void evaluate(std::string const&);

    //---------------------------- PROPERTIES ----------------------------------

    /// change the name of the simulation
    void rename(std::string const&);
    
    /// read an Object reference and return the corresponding Object (`tag` is set)
    Object * readReference(Inputter&, ObjectTag& tag);

    /// read a Fiber reference and return the corresponding Object (`tag` is set)
    Fiber * readFiberReference(Inputter&, ObjectTag& tag);
    
    
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
    Property* newProperty(const std::string&, const std::string&, Glossary&);
    
    /// export all Properties to speficied file
    void writeProperties(std::ostream&, bool prune) const;

    /// export all Properties to a new file with specified name
    void writeProperties(char const* filename, bool prune) const;
    
    /// export all Properties using default file name
    void writeProperties(bool prune) const;
    
    /// load the properties contained in the standard output property file
    void loadProperties(char const* filename);

    /// load the properties contained in the standard output property file
    void loadProperties();

    //---------------------------- LOAD OBJECTS --------------------------------

    /// default name for output trajectory file
    static const char TRAJECTORY[];

    /// current file format (check history in `simul_file.cc`)
    static constexpr unsigned currentFormatID = 56;
    
    /// class for reading trajectory file
    class InputLock;
    
    /// read objects from file, and add them to the simulation state
    int readObjects(Inputter&, ObjectSet* subset);

    /// load sim-world from the named file
    int loadObjects(char const* filename);
    
    /// import objects from file, and delete objects that were not referenced in the file
    int reloadObjects(Inputter&, bool prune = 1, ObjectSet* subset = nullptr);

    /// write sim-world to specified file
    void writeObjects(Outputter&) const;
    
    /// write sim-world in binary or text mode, appending to existing file or creating new file
    void writeObjects(std::string const& filename, bool append, bool binary) const;
    
    //----------------------------- REPORTING ----------------------------------
    
    /// calls report_one
    void poly_report(std::ostream&, std::string, Glossary&, int) const;
    
    /// calls report
    void mono_report(std::ostream&, std::string const&, Glossary&, int) const;

    /// calls report_one
    void report_one(std::ostream&, std::string const&, Glossary&) const;
    
    /// call one of the report function
    void report_one(std::ostream&, std::string const&, Property const*, std::string const&, Glossary&) const;

    
    /// print time
    void reportTime(std::ostream&) const;
    
    /// give a short inventory of the simulation state, obtained from ObjectSet::report()
    void reportInventory(std::ostream&) const;
    
    /// give a summary of the Simul
    void reportSimul(std::ostream&) const;
    
    /// print the length and the points of each fiber
    void reportFiber(std::ostream&, Fiber const*) const;

    /// print the length and the points of each fiber
    void reportFibers(std::ostream&, Property const*) const;
    
    /// print the length and the points of each fiber, sorted from longest to shortest
    void reportFibersSorted(std::ostream&, Property const*) const;

    /// print the coordinates of the vertices of each fiber
    void reportFiberPoints(std::ostream&, Property const*) const;
    
    /// print the coordinates of the vertices of each fiber
    void reportFiberDisplacement(std::ostream&, Property const*) const;

    /// print the positions and the states of ends `end` of all fibers
    void reportFiberEnds(std::ostream&, FiberEnd end, Property const*) const;
    
    /// print number of fibers in each state of specified end
    void reportFiberEndState(std::ostream&, FiberEnd end, Property const*) const;

    /// print the mean and standard deviation of vertices for each class of fiber
    void reportFiberMoments(std::ostream&) const;

    /// print average age and standard deviation for each class of fiber
    void reportFiberAge(std::ostream&) const;
    
    /// print average length and standard deviation for each class of fiber
    void reportFiberLengths(std::ostream&, FiberProp const*) const;

    /// print average length and standard deviation for each class of fiber
    void reportFiberLengths(std::ostream&, Property const*) const;
    
    /// print length distribution for each class of fiber
    void reportFiberLengthHistogram(std::ostream&, Glossary&) const;
    
    /// print coordinates of speckles that follow a frozen random sampling
    void reportFiberSpeckles(std::ostream&, Glossary&) const;
    
    /// print coordinates of points randomly and freshly distributed
    void reportFiberSamples(std::ostream&, Glossary&) const;
    
    /// print the coordinates and forces on the vertices of each fiber
    void reportFiberForces(std::ostream&) const;

    /// print Fiber tensions along certain planes defined in Glossary
    void reportFiberTension(std::ostream&, Glossary&) const;
    
    
    /// print number of kinks in each class of Fiber
    void reportFiberSegments(std::ostream&) const;

    /// document bending energy
    void reportFiberBendingEnergy(std::ostream&) const;
    
    /// document end-to-end distance in each class of Fiber
    void reportFiberExtension(std::ostream&) const;
    
    /// document nematic order (alignement) of Fiber
    void reportFiberNematic(std::ostream&) const;

    
    /// print component of forces experienced by Fibers due to confinement
    void reportFiberConfineForce(std::ostream& out) const;

    /// print radial component of forces experienced by Fibers due to confinement
    real reportFiberConfinement(std::ostream& out) const;

    /// print summary of Fiber's lattice quantities
    void reportFiberLattice(std::ostream&, Property const*) const;
    
    /// print summary of Fiber's lattice quantities
    void reportFiberMeshAverage(std::ostream&, bool density, Property const*) const;
    
    /// print summary of Fiber's lattice quantities
    void reportFiberMesh(std::ostream&, bool density, Property const*) const;
    
    
    /// print position of hands bound to fibers
    void reportFiberHands(std::ostream&) const;
    
    /// print position of bound hands that are associated with stiffness
    void reportFiberLinks(std::ostream&) const;

    /// print interection abscissa between fibers
    void reportFiberConnectors(std::ostream&, Glossary&) const;

    /// print interection abscissa between fibers
    void reportNetworkBridges(std::ostream&, Glossary&) const;
 
    /// print network surface area
    void reportNetworkSize(std::ostream&) const;

    /// print positions of interection between two fibers
    void reportFiberIntersections(std::ostream&, Glossary&) const;


    /// print Organizer positions
    void reportOrganizer(std::ostream&) const;

    /// print Aster positions
    void reportAster(std::ostream&) const;
    
    /// print Bead positions 
    void reportBeadSingles(std::ostream&) const;

    /// print Bead positions
    void reportBeadPosition(std::ostream&, Property const*) const;

    /// print Solid positions 
    void reportSolidPosition(std::ostream&, Property const*) const;
    
    /// print Solid positions
    void reportSolidOrientation(std::ostream&, Property const*) const;

    /// print Solid's anchored Hands
    void reportSolidHands(std::ostream&, Property const*) const;

    /// print state of Couples 
    void reportCouple(std::ostream&, Property const*) const;
 
    /// print state of Couples
    void reportCoupleList(std::ostream&, Property const*) const;

    /// print state of Couples
    void reportCoupleAnatomy(std::ostream&, Property const*) const;
    
    /// print position of Couples of a certain kind
    void reportCoupleState(std::ostream&, Property const*) const;
    
    /// print position of active Couples of a certain kind
    void reportCoupleActive(std::ostream&, Property const*) const;
    
    /// print position and forces of doubly-attached Couples
    void reportCoupleLink(std::ostream&, Property const*) const;
    
    /// print configurations of doubly-attached Couples
    void reportCoupleConfiguration(std::ostream&, Property const*, Glossary&) const;

    /// print agregate properties of Couple force
    void reportCoupleForce(std::ostream&, Property const*) const;
    
    /// print histogram of Couples force
    void reportCoupleForceHistogram(std::ostream&, Glossary&) const;

    /// print state of Singles
    void reportSingle(std::ostream&, Property const*) const;
    
    /// print position of Singles of a certain kind
    void reportSingleState(std::ostream&, Property const*) const;

    /// print position of Singles
    void reportSinglePosition(std::ostream&, Property const*) const;
   
    /// print force of attached Singles
    void reportSingleLink(std::ostream&, Property const*) const;
    
    /// print agregate properties of Singles force
    void reportSingleForce(std::ostream&, Property const*) const;

    /// print state of Couples 
    void reportSpherePosition(std::ostream&, Property const*) const;

    /// print something about Spaces
    void reportSpace(std::ostream&) const;
  
    /// print force on Spaces
    void reportSpaceForce(std::ostream&) const;

    /// print something about Fields
    void reportField(std::ostream&) const;
    
    //------------------------- CLUSTER ANALYSIS -------------------------------
    
    /// set Mecable's flag() to unique values, return highest value attributed + 1
    ObjectFlag setUniqueFlags() const;

    /// call `flag(f)` for all Mecable
    void setFlags(ObjectFlag f) const;
    
    /// replace all occurence of flag `f` by `g`
    void changeFlag(ObjectFlag f, ObjectFlag g) const;

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
    void reportClusters(std::ostream&, Glossary&) const;

    //------------------------- PROJECT SPECIFIC -------------------------------

    /// flag the fibers that appear to constitute a ring
    size_t flagRing() const;
    
    /// extract radius and length of ring
    void analyzeRing(ObjectFlag, real& length, real& radius) const;
    
    /// estimates if Fibers form a connected ring around the Z-axis
    void reportRing(std::ostream&) const;

    /// custom report for Platelet project
    void reportPlatelet(std::ostream&) const;
    
    /// print Aster & Spindle indices
    void reportIndices(std::ostream&) const;

    /// print number of Fibers pointing left and right that intersect plane YZ at different X positions
    void reportProfile(std::ostream&) const;

    /// a special print for Romain Gibeaux
    void reportAshbya(std::ostream&) const;
    
    /// analysis of MT collisions in the plant cortex
    void reportFiberCollision(std::ostream&, Property const*, Glossary&) const;

    /// print something
    void reportSomething(std::ostream&) const;

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

