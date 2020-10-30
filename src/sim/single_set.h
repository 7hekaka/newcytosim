// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SINGLE_SET_H
#define SINGLE_SET_H

#include "object_set.h"
#include "single.h"
class PropertyList;

/// a list of pointers to Single
typedef Array<Single *> SingleList;


/// Set for Single
/**
 A Single is stored in one of 2 ObjectPool, depending on its state:
 - fList = free,
 - aList = attached.
 .
 Each list is accessible via its head firstF() and firstA(),
 and subsequent objects obtained with next() are in the same state.
 This way, the state of the Single are known when accessing them.
 
 A Single is automatically transfered to the appropriate list,
 if its Hand binds or unbinds. This is done by the HandMonitor.
 */
class SingleSet: public ObjectSet
{
public:
    
    /// flags affecting how free Couples are read/written from file
    static bool prune_free, skip_free;
    
private:
    
    /// List for non-attached Singles (f=free)
    ObjectPool     fList;
    
    /// List for attached Singles (a=attached)
    ObjectPool     aList;
    
    /// holds the property and the list of singles
    typedef std::pair<SingleProp const*, SingleList> SingleReserve;

    /// an array of SingleReserveList
    typedef std::vector<SingleReserve> SingleReserveList;
    
    /// uniReserves[p] contains the Single with ( property()->number() == p ) that are diffusing
    SingleReserveList uniReserves;
    
    /// flag to enable `fast_diffusion` attachment algorithm
    bool          uniEnabled;

    /// initialize `fast_diffusion` attachment algorithm
    bool          uniPrepare(PropertyList const& properties);
    
    /// gather all Single with `fast_diffusion` in reserve lists
    Single*       uniCollect(Single*);

    /// ensures that `can` holds `cnt` Singles, creating them of specified SingleProp
    void          uniRefill(SingleList& can, size_t cnt, SingleProp const*);

    /// attach Singles from `can` on locations specified in `loc`
    void          uniAttach(Array<FiberSite>& can, SingleList& loc);
    
    /// `fast_diffusion` attachment assuming that free Singles are uniformly distributed
    void          uniAttach(FiberSet const&);
    
    /// release Single from reserve lists
    void          uniRelax();

public:
        
    ///creator
    SingleSet(Simul& s) : ObjectSet(s), uniEnabled(false) {}
    
    //--------------------------

    /// identifies the class
    static std::string title() { return "single"; }
    
    /// create a new property of category `cat` for a class `name`
    Property *  newProperty(const std::string& cat, const std::string& name, Glossary&) const;
    
    /// create objects of class `name`, given the options provided in `opt`
    ObjectList  newObjects(const std::string& name, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object *    newObject(ObjectTag, size_t);
    
    /// print a summary of the content (nb of objects, class)
    void        report(std::ostream&) const;

    /// write
    void        write(Outputter&) const;

    //--------------------------

    /// add object
    void          link(Object *);

    /// remove object
    void          unlink(Object *);
    
    /// reassign Single to different sublist following attachement of Hand
    void          relinkA(Single *);
    
    /// reassign Single to different sublist following detachment of Hand
    void          relinkD(Single *);
    
    /// delete an attached Single
    void          deleteA(Single *);
    
    /// detach an attached Single
    void          detachA(Single *);

    /// create Wrists anchored on given Mecable
    ObjectList    makeWrists(Mecable const*, size_t, size_t, std::string&);

    /// return all Wrists anchored on `obj`
    SingleList    collectWrists(Object const*) const;
    
    /// remove all Wrists anchored on `obj`
    void          removeWrists(Object const*);
    
    
    ///returns the first free Single
    Single *      firstF()       const { return static_cast<Single*>(fList.front()); }
    
    ///returns the first bound Single
    Single *      firstA()       const { return static_cast<Single*>(aList.front()); }
    
    /// return pointer to the Object of given ID, or zero if not found
    Single *      findID(ObjectID n) const { return static_cast<Single*>(inventory.get(n)); }
    
    /// first Single in inventory
    Single *      firstID()           const { return static_cast<Single*>(inventory.first()); }
    
    /// next Single in inventory
    Single *      nextID(Single const* obj) const { return static_cast<Single*>(inventory.next(obj)); }

    /// collect all objects
    ObjectList    collect() const;
    
    /// collect objects for which func(obj, val) == true
    ObjectList    collect(bool (*func)(Object const*, void const*), void const*) const;

    /// collect objects for which func(obj, val) == true
    size_t        count(bool (*func)(Object const*, void const*), void const*) const;

    /// erase all Object and all Property
    void          erase();
    
    /// number of unattached Simgles
    size_t        sizeF() const { return fList.size(); }
    
    /// number of attached Singles
    size_t        sizeA() const { return aList.size(); }

    /// number of elements
    size_t        size()  const { return fList.size() + aList.size(); }
    
    /// mix order of elements
    void          shuffle();

    /// prepare for step()
    void          prepare(PropertyList const& properties);
    
    /// Monte-Carlo step
    void          step();
    
    /// cleanup at end of simulation period
    void          relax() { uniRelax(); }
    
    /// modulo the position (periodic boundary conditions)
    void          foldPosition(Modulo const*) const;

    //--------------------------
    
    /// mark object before import
    void          freeze(ObjectFlag f);
    
    /// delete marked object after import
    void          prune(ObjectFlag f);
    
    /// unmark objects after import
    void          thaw();
    
    /// check internal consistency, returns 0 if everything is OK
    int           bad() const;

};


#endif

