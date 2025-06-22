// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University
#ifndef SINGLE_SET_H
#define SINGLE_SET_H

#include <vector>

#include "object_set.h"
#include "single.h"
class PropertyList;

/// a list of pointers to Single
typedef Array<Single *> SingleList;


/// Set holding objects of class Single
/**
 A Single is stored in one of 2 ObjectPool, depending on its state:
 - fList = free,
 - aList = attached.
 .
 Each list is accessible via its head firstF() and firstA(),
 and subsequent objects obtained with next() are in the same state.
 This way, the state of the Single are known when accessing them.
 
 A Single is automatically transferred to the appropriate list,
 if its Hand binds or unbinds. This is done by the HandMonitor.
 */
class SingleSet: public ObjectSet
{
private:
    
    /// List for non-attached Singles (f=free)
    ObjectPool fList;
    
    /// List for attached Singles (a=attached)
    ObjectPool aList;
    
    /// list of Properties for which `fast_diffusion > 0`
    std::vector<SingleProp const*> uniSingles;

    /// gather all Single with `fast_diffusion` in reserve lists
    void uniStepCollect(Single*);

    /// ensures that `can` holds `cnt` Singles, creating them of specified SingleProp
    void uniRefill(SingleProp const*, size_t cnt);

    /// attach Singles from `can` on locations specified in `loc`
    void uniAttach(FiberSiteList& can, SingleStock& loc);
    
    /// `fast_diffusion` attachment assuming that free Singles are uniformly distributed
    void uniAttach(FiberSet const&);

    /// release Single from reserve lists
    void uniRelax();

    /// save free Single for which `save_unbound > 0`
    void writeSomeObjects(Outputter&) const;

public:

    /// creator
    SingleSet(Simul& s) : ObjectSet(s) {}
    
    //--------------------------

    /// identifies the class
    static std::string title() { return "single"; }
    
    /// create a new property of category `cat` for a class `name`
    Property * newProperty(std::string const& cat, std::string const& name, Glossary&) const;
    
    /// create objects specified by Property, given options provided in `opt`
    ObjectList newObjects(Property const*, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object *   newObject(ObjectTag, PropertyID);
    
    /// print a summary of the content (nb of objects, class)
    void report(std::ostream&) const;

    /// write objects
    void writeSet(Outputter&) const;

    //--------------------------
    
    /// return one sublist where Couple should be linked
    ObjectPool& sublist(Single const* obj) { return obj->attached()?aList:fList; }

    /// add object
    void link(Object *);

    /// remove object
    void unlink(Object *);

    /// reassign Single to different sublist following attachement of Hand
    void relinkA(Single *);
    
    /// reassign Single to different sublist following detachment of Hand
    void relinkD(Single *);

    /// return all Wrists anchored on given object
    SingleList collectWrists(Object const*) const;
    
    /// detach all Wrists anchored on given object
    void detachWrists(Object const*);

    /// delete all Wrists anchored on given object
    void deleteWrists(Object const*);

    /// create Wrists anchored on given Mecable
    void makeWrists(ObjectList&, Mecable const*, index_t, index_t, std::string const&);

    /// create Single attached to the beads
    void distributeWrists(ObjectList&, SingleProp const*, size_t cnt, std::string const&) const;
    
    /// return the first free Single
    Single * firstF() const { return static_cast<Single*>(fList.front()); }
    
    /// return the first bound Single
    Single * firstA() const { return static_cast<Single*>(aList.front()); }
    
    /// return pointer to the Object of given ID, or zero if not found
    Single * identifyObject(ObjectID n) const { return static_cast<Single*>(inventory_.get(n)); }
    
    /// first Single in inventory
    Single * firstID() const { return static_cast<Single*>(inventory_.first()); }
    
    /// returns Single immediately following 'obj' in inventory
    Single * nextID(Single const* obj) const { return static_cast<Single*>(inventory_.next(obj)); }

    /// collect all objects
    ObjectList collect() const;
    
    /// collect objects for which func(obj, val) == true
    ObjectList collect(bool (*func)(Object const*, void const*), void const*) const;

    /// collect objects for which func(obj, val) == true
    size_t count(bool (*func)(Object const*, void const*), void const*) const;

    /// erase all objects
    void erase();
    
    /// detach all Hands
    void detachAll();

    /// number of unattached Singles
    size_t sizeF() const { return fList.size(); }
    
    /// number of attached Singles
    size_t sizeA() const { return aList.size(); }

    /// number of elements
    size_t size()  const { return fList.size() + aList.size(); }
    
    /// mix order of elements
    void shuffle();
    
    /// prepare for steps()
    void prepare();
    
    /// Monte-Carlo step
    void steps();
    
    /// Monte-Carlo step without Hand attachment
    void stepsSkippingUnattached();
    
    /// cleanup at end of simulation period
    void relax() { uniRelax(); }
    
    /// attach Single to nearby fibers to approximate an equilibrated state
    void equilibrate(FiberSet const&, SingleProp const*);

    /// attach Single to nearby fibers to approximate an equilibrated state
    void equilibrate(SingleProp const*);

    /// attach Single to nearby fibers to approximate an equilibrated state
    void equilibrate();

    /// bring all objects to centered image using periodic boundary conditions
    void foldPositions(Modulo const*) const;
    
    //--------------------------

    /// initialize `fast_diffusion` attachment algorithm
    void uniPrepare(PropertyList const& properties);

    /// total count in reserves
    size_t countStocks() const;
    
    /// print number of elements in each reserve bin
    void infoStocks(std::ostream& os) const;

    //--------------------------
    
    /// return a Single from the reserve, or made by newSingle()
    Single * makeSingle(SingleProp const*);

    /// create unattached Singles
    void makeSingles(SingleProp const*, size_t cnt);

    /// create unattached Singles
    void makeSingles(size_t cnt[], PropertyID n_cnt);
    
    /// link unattached Single
    void addFreeSingle(Single*);
    
    /// register given Single, which must be unattached
    Single * addFreeSingle(SingleProp const*, Vector const&);

    //--------------------------

    /// prepare all objects before reading
    void freeze();
    
    /// move Singles into reserve lists, instead of deleting them
    void defrostStore();

    /// detach objects that were not updated during import
    void reheat(size_t cnt[], PropertyID n_cnt);
    
    /// detach objects that were not updated during import
    void reheat();
    
    /// check internal consistency, returns 0 if everything is OK
    int bad() const;

};


#endif

