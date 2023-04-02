// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University
#ifndef SINGLE_SET_H
#define SINGLE_SET_H

#include <vector>

#include "object_set.h"
#include "single.h"
class PropertyList;

/// a list of pointers to Single
typedef Array<Single *> SingleList;

/// holds a list of Single of the same type
/**
 This list is build using 'Object::next()' and thus can only hold
 objects that are not linked in SingleSet already (see below).
 It is a single-linked list and addition/removal is made on the same end.
 SingleReserve is used for the 'fast_diffusion' algorithm.
 */
class SingleReserve
{
    /// All Single in the list should be of the same type
    SingleProp const* property_;
    /// Pointer to first member in list
    Single * head_;
    /// Number of Couples in list
    size_t count_;
    
public:
    
    /// constructor
    SingleReserve() { count_ = 0; head_ = nullptr; property_ = nullptr; }
    
    /// number of objects stored
    size_t size() const { return count_; }
    
    /// return property
    SingleProp const* property() { return property_; }
    
    /// set Property
    void property(SingleProp const* p) { property_ = p; }
    
    /// first object
    Single * head() const { return head_; }
    
    /// add object
    void push(Single* arg) { arg->Object::next(head_); head_ = arg; ++count_; }
    
    /// remove first object in list
    void pop() { head_ = head_->next(); --count_; }
    
    /// delete all objects
    void erase()
    {
        Single * obj = head_;
        while ( obj )
        {
            pop();
            obj->objset(nullptr);
            delete(obj);
            obj = head();
        }
    }
};


/// Set for Single
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

    /// an array of SingleReserveList
    typedef std::vector<SingleReserve> SingleReserveList;
    
    /// uniReserves[p] holds Singles with ( property()->number() == p )
    SingleReserveList uniReserves;
    
    /// flag to enable `fast_diffusion` attachment algorithm
    bool uniEnabled;

    /// initialize `fast_diffusion` attachment algorithm
    bool uniPrepare(PropertyList const& properties);
    
    /// gather all Single with `fast_diffusion` in reserve lists
    void uniStepCollect(Single*);

    /// ensures that `can` holds `cnt` Singles, creating them of specified SingleProp
    void uniRefill(SingleReserve& can, size_t cnt, SingleProp const*);

    /// attach Singles from `can` on locations specified in `loc`
    void uniAttach(Array<FiberSite>& can, SingleReserve& loc);
    
    /// `fast_diffusion` attachment assuming that free Singles are uniformly distributed
    void uniAttach(FiberSet const&);

    /// release Single from reserve lists
    void uniRelax();

    /// save free Single for which `fast_diffusion == 0`
    void writeF_skip(Outputter&) const;
    
public:
    
    /// total count in reserves
    size_t reserved() const;
    
    /// total count in reserves
    size_t reserved(size_t i) const { return uniReserves[i].size(); }

    /// flags to skip unattached Single in trajectory file
    mutable int skip_now;
    
    /// creator
    SingleSet(Simul& s) : ObjectSet(s), uniEnabled(0), skip_now(0) {}
    
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

    /// add object
    void link(Object *);

    /// remove object
    void unlink(Object *);
    
    /// link unattached Single
    void linkF(Single * s) { assert_true(!s->attached()); fList.push_back(s); }

    /// reassign Single to different sublist following attachement of Hand
    void relinkA(Single *);
    
    /// reassign Single to different sublist following detachment of Hand
    void relinkD(Single *);

    /// create Wrists anchored on given Mecable
    void makeWrists(ObjectList&, Mecable const*, size_t, size_t, std::string&);

    /// return all Wrists anchored on `obj`
    SingleList collectWrists(Object const*) const;
    
    /// detach all Wrists anchored on `obj`
    void detachWrists(Object const*);

    /// delete all Wrists anchored on `obj`
    void deleteWrists(Object const*);
    
    /// create Single attached to the beads
    ObjectList distributeWrists(SingleProp const*, size_t cnt, std::string const&) const;
    
    ///returns the first free Single
    Single * firstF() const { return static_cast<Single*>(fList.front()); }
    
    ///returns the first bound Single
    Single * firstA() const { return static_cast<Single*>(aList.front()); }
    
    /// return pointer to the Object of given ID, or zero if not found
    Single * findID(ObjectID n) const { return static_cast<Single*>(inventory_.get(n)); }
    
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

    /// erase all Object and all Property
    void erase();
    
    /// number of unattached Simgles
    size_t sizeF() const { return fList.size(); }
    
    /// number of attached Singles
    size_t sizeA() const { return aList.size(); }

    /// number of elements
    size_t size()  const { return fList.size() + aList.size(); }
    
    /// mix order of elements
    void shuffle();

    /// prepare for step()
    void prepare(PropertyList const& properties);
    
    /// Monte-Carlo step
    void step();
    
    /// Monte-Carlo step without Hand attachment
    void stepSkipUnattached();
    
    /// cleanup at end of simulation period
    void relax() { uniRelax(); }
    
    /// bring all objects to centered image using periodic boundary conditions
    void foldPositions(Modulo const*) const;

    //--------------------------
    
    /// unlink all objects before import
    void freeze();
    
    /// detach objects that were not updated during import
    void reheat(size_t cnt[], size_t n_cnt);
    
    /// detach objects that were not updated during import
    void reheat();
    
    /// check the anchoring Mecables of Wrists
    void deleteInvalidWrists();
    
    /// check internal consistency, returns 0 if everything is OK
    int bad() const;

};


#endif

