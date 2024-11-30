// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#ifndef FIBER_SET_H
#define FIBER_SET_H

#include "dim.h"
#include "object_set.h"
#include "fiber.h"
#include "fiber_site.h"

class FiberProp;
class CoupleProp;


/// a list of Fiber
/**
 The FiberSet stores Fiber, and derived classes.
 It constains many algorithms that specifically deal with Fibers,
 including many function to calculate average length, nematic organization, etc.
 */
class FiberSet : public ObjectSet
{
private:
    
    FiberSet();

public:
    
    /// creator
    FiberSet(Simul& s) : ObjectSet(s) {}
    
    //--------------------------
    
    /// identifies the class
    static std::string title() { return "fiber"; }
    
    /// create a new property of category `cat` for a class `name`
    Property * newProperty(std::string const& cat, std::string const& name, Glossary&) const;
    
    /// create objects specified by Property, given options provided in `opt`
    Fiber * newFiber(ObjectList&, FiberProp const*, Glossary& opt) const;
    
    /// create objects specified by Property, given the options provided in `spec`
    Fiber * newFiber(ObjectList&, FiberProp const*, std::string const& spec) const;

    /// create objects specified by Property, given options provided in `opt`
    ObjectList newObjects(Property const* p, Glossary& opt);

    /// create a new object (used for reading trajectory file)
    Object * newObject(ObjectTag, PropertyID);
    
    /// write all Objects to file
    void writeSet(Outputter&) const;
        
    /// print a summary of the content (nb of objects, class)
    void report(std::ostream&) const;

    //--------------------------

    /// first Fiber
    Fiber * first() const { return static_cast<Fiber*>(pool_.front()); }
    
    /// last Fiber
    Fiber * last() const { return static_cast<Fiber*>(pool_.back()); }
    
    /// first Fiber in inventory
    Fiber * firstID() const { return static_cast<Fiber*>(inventory_.first()); }

    /// returns Fiber immediately following 'obj' in inventory
    Fiber * nextID(Fiber const* obj) const { return static_cast<Fiber*>(inventory_.next(obj)); }

    /// return pointer to the Object of given ID, or zero if not found
    Fiber * identifyObject(ObjectID n) const { return static_cast<Fiber*>(inventory_.get(n)); }

    /// Cut all segments intersecting the plane defined by <em> n.pos + a = 0 </em>
    void planarCut(Vector const& n, real a, state_t stateP, state_t stateM, real min_len);

    /// Cut fibers in the list that are intersecting the plane defined by <em> n.pos + a = 0 </em>
    void planarCut(ObjectList&, Vector const& n, real a, state_t stateP, state_t stateM, real min_len);
    
    /// Special code to trim a spindle from both ends
    void shortenSpindle(real, real) const;
    
    /// get ready to do steps()
    void prepare();

    /// Monte-Carlo step for every Fiber
    void steps();
    
    /// bring all objects to centered image using periodic boundary conditions
    void foldPositions(Modulo const*) const;
    
    /// find intersections between fibers in entire network, within given threshold
    void allIntersections(FiberSiteList&, FiberSiteList&, real max_distance) const;
    
    /// find intersections between fibers in entire network, within given threshold
    void allIntersections0(FiberSiteList&, FiberSiteList&, real max_distance) const;

    /// set random sites along the fibers, separated on average by `gap`
    void uniFiberSites(FiberSiteList&, real gap) const;
    
    /// a random site on the fibers, uniformly distributed over all Fibers
    FiberSite randomSite() const;
    
    /// a random site on the fibers of specified class, uniformly distributed
    FiberSite randomSite(FiberProp const*) const;

    /// a site on a fiber, as specified by Glossary[key]
    FiberSite someSite(std::string const& key, Glossary&) const;

    /// set random sites on newly polymerized Fiber sites at the plus end
    void newFiberSitesP(FiberSiteList&, real gap) const;
    
    /// set random sites on newly polymerized Fiber sites at the minus end
    void newFiberSitesM(FiberSiteList&, real gap) const;
    
    /// reverse the polarity of all fibers
    void flipFiberPolarity(FiberProp *);
    
    /// update object after import
    void updateFibers() const;
    
    //--------------------------------------------------------------------------
    
    /// total length of Fiber 
    real totalLength() const;

    /// total length of Fiber for Fibers with given FiberProp
    real totalLength(FiberProp const *) const;
    
    /// calculate: number of fibers, mean, variance, min and max of fiber length
    static void infoLength(ObjectList const&, size_t& cnt, real& avg, real& var, real& mn, real& mx, real& off);
    
    /// calculate: number of fibers, mean, variance, min and max of fiber length
    static void infoBirthtime(ObjectList const&, size_t& cnt, real& avg, real& var, real& mn, real& mx);

    /// calculate: number of fibers, number of joints and number of kinks
    static void infoSegments(ObjectList const&, size_t& cnt, size_t& points, real&, real&, real&);
    
    /// calculate: number of fibers, number of joints and number of kinks
    static size_t nbKinks(ObjectList const&);

    /// calculate center of gravity G, average of minus end and plus end
    static real infoPosition(ObjectList const& objs, Vector& M, Vector& G, Vector& P);

    /// calculate the nematic directors, return the nematic scalar order parameter
    static real infoNematic(ObjectList const&, real vec[9]);

    /// calculate the nematic directors, return the nematic scalar order parameter
    static real infoOrthoNematic(ObjectList const&, real vec[9], Space const*);

    /// calculate center of gravity G, and principal components axes
    static int infoComponents(ObjectList const&, real& sum, real avg[3], real mom[9], real vec[9]);

    /// Count Fibers intersecting the plane defined by <em> n.pos + a = 0 </em>
    void infoPlane(int& np, int& na, Vector const& n, real a) const;
    
    /// Calculate characteristics of bendingEnergy()
    static void infoBendingEnergy(ObjectList const&, size_t& cnt, real& avg, real& var);
    
    /// sum Lagrange multipliers for segments that intersect the plane <em> n.pos + a = 0 </em>
    void infoTension(size_t& cnt, real& sum, real& inf, real& sup, Vector const& n, real a) const;
    
    /// sum Lagrange multipliers for all fibers
    void infoTension(size_t& cnt, real& sum, real& inf, real& sup) const;

    /// Calculate spindle indices
    void infoSpindle(real& ixa, real& ixs, Vector const& n, real a, real m, real da) const;

    /// Calculate averaged distance from origin - for all vertices
    void infoRadius(size_t&, real& rad) const;
    
    /// Calculate averaged distance from origin - for fiber ends
    void infoRadius(size_t&, real& rad, FiberEnd) const;

};


#endif

