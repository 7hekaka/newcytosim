// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef COUPLE_SET_H
#define COUPLE_SET_H

#include <vector>

#include "object_set.h"
#include "couple.h"
#include "couple_prop.h"

class PropertyList;

/// a list of pointers to Single
typedef Array<Couple *> CoupleList;


/// Set for Couple
/**
 A Couple is stored in one of 4 ObjectPool, depending on its state:
 - ffList = hand1 and hand2 unattached,
 - afList = hand1 attached, hand2 unattached,
 - faList = hand1 unattached, hand2 attached,
 - aaList = hand1 and hand2 attached [the couple makes a `link`].
 .
 The head of each list is accessible via firstFF() and firstFA(), firstAF() and firstAA(),
 and subsequent objects obtained with next() are in the same state.
 This makes iterating through the list of Couple more efficient.
 
 A Couple is automatically transferred to the appropriate list,
 if one of its Hand binds or unbinds. This is done by the HandMonitor.
 */
class CoupleSet: public ObjectSet
{
private:
    
    /// list of Couple which are not attached (f=free)
    ObjectPool ffList;
    
    /// list of Couple with only one side attached (a=attached, f=free)
    ObjectPool afList, faList;
    
    /// list of Couple with both sides attached (a=attached)
    ObjectPool aaList;
    
    /// return one of ffList, afList, faList, aaList, corresponding to given states
    ObjectPool& sublist(bool attached1, bool attached2)
    {
        if ( attached1 )
        {
            if ( attached2 ) return aaList; else return afList;
        }
        else
        {
            if ( attached2 ) return faList; else return ffList;
        }
    }

    /// holds the property and the list of singles
    typedef std::pair<CoupleProp const*, CoupleList> CoupleReserve;

    /// an array of SingleReserveList
    typedef std::vector<CoupleReserve> CoupleReserveList;
    
    /// uniReserves[p] holds Couples with ( property()->number() == p+1 )
    CoupleReserveList uniReserves;
    
    /// flag to enable `fast_diffusion` attachment algorithm
    bool uniEnabled;
    
    /// initialize `fast_diffusion` attachment algorithm
    bool uniPrepare(PropertyList const& properties);
    
    /// gather all Couple with `fast_diffusion` in reserve lists
    template <void (Couple::*FUNC)()> void step_collect(Couple*);

    /// ensures that `can` holds `cnt` Couple, creating them of specified CoupleProp
    void uniRefill(CoupleList& can, size_t cnt, CoupleProp const*);

    /// attach Hand1 of Couple from `can` on locations specified in `loc`
    void uniAttach1(Array<FiberSite>& loc, CoupleList& can);
    
    /// attach Hand2 of Couple from `can` on locations specified in `loc`
    void uniAttach2(Array<FiberSite>& loc, CoupleList& can);
    
    /// attach both Hands of `nb` Couple at crossing points specified by first argument
    void uniAttach12(Array<FiberSite>&, Array<FiberSite>&, CoupleList&, size_t nb);

    /// `fast_diffusion` attachment assuming that free Singles are uniformly distributed
    void uniAttach(FiberSet const&);
    
    /// release Couples from reserve lists
    void uniRelax();
    
    /// detach objects that were not updated during import
    void pruneDelete();
    
    /// detach objects that were not updated during import
    void pruneDetach();

public:
        
    /// flags to skip unattached Couple in trajectory file
    mutable int skip_now, prune_mode;

    ///creator
    CoupleSet(Simul& s) : ObjectSet(s), uniEnabled(0), skip_now(0), prune_mode(0) {}
    
    //--------------------------
    
    /// identifies the set
    static std::string title() { return "couple"; }
    
    /// create a new property of category `cat` for a class `name`
    Property * newProperty(const std::string& cat, const std::string& name, Glossary&) const;
    
    /// create objects of class `name`, given the options provided in `opt`
    ObjectList newObjects(const Property*, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object * newObject(ObjectTag, PropertyID);
    
    /// save free Couples
    void writeFF(Outputter&) const;
    
    /// save attached Couples
    void writeAF(Outputter&) const;
    
    /// save attached Couples
    void writeFA(Outputter&) const;

    /// save bridging Couples
    void writeAA(Outputter&) const;

    /// save objects
    void writeSet(Outputter&) const;

    //--------------------------

    /// add object (should be a Couple)
    void link(Object *);
    
    /// remove object (should be a Couple)
    void unlink(Object *);

    /// link unattached Couple into FF list
    void linkFF(Couple * c) { assert_true(!c->attached1()&&!c->attached2()); ffList.push_back(c); }

    /// reassign Couple to sublist following attachement of Hand 1
    void relinkA1(Couple *);
    /// reassign Couple to sublist following detachment of Hand 1
    void relinkD1(Couple *);
    /// reassign Couple to sublist following attachement of Hand 2
    void relinkA2(Couple *);
    /// reassign Couple to sublist following detachment of Hand 2
    void relinkD2(Couple *);

    /// first unattached Couple
    Couple * firstFF() const{ return static_cast<Couple*>(ffList.front()); }
    /// first Couple attached only by cHand1
    Couple * firstAF() const{ return static_cast<Couple*>(afList.front()); }
    /// first Couple attached only by cHand2
    Couple * firstFA() const{ return static_cast<Couple*>(faList.front()); }
    /// first Couple attached by both hands
    Couple * firstAA() const{ return static_cast<Couple*>(aaList.front()); }

    /// last unattached Couple
    Couple * lastFF()  const{ return static_cast<Couple*>(ffList.back()); }
    /// last Couple attached by cHand1
    Couple * lastAF()  const{ return static_cast<Couple*>(afList.back()); }
    /// last Couple attached by cHand2
    Couple * lastFA()  const{ return static_cast<Couple*>(faList.back()); }
    /// last Couple attached by both hands
    Couple * lastAA()  const{ return static_cast<Couple*>(aaList.back()); }

    /// number of free Couples
    size_t sizeFF()  const{ return ffList.size(); }
    /// number of Couples attached by cHand1 only
    size_t sizeAF()  const{ return afList.size(); }
    /// number of Couples attached by cHand2 only
    size_t sizeFA()  const{ return faList.size(); }
    /// number of Couples attached by one hand
    size_t sizeA()   const{ return faList.size() + afList.size(); }
    /// number of Couples attached by both hands
    size_t sizeAA()  const{ return aaList.size(); }
    /// total number of elements
    size_t size()    const{ return ffList.size() + faList.size() + afList.size() + aaList.size(); }

    /// erase all Object and all Property
    void erase();
    
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
    
    
    /// return pointer to the Object of given ID, or zero if not found
    Couple * findID(ObjectID n)        const { return static_cast<Couple*>(inventory_.get(n)); }
    
    /// first Couple in inventory
    Couple * firstID()                 const { return static_cast<Couple*>(inventory_.first()); }
    
    /// returns Couple immediately following 'obj' in inventory
    Couple * nextID(Couple const* obj) const { return static_cast<Couple*>(inventory_.next(obj)); }
    
    /// collect all objects
    ObjectList collect() const;
    
    /// collect objects for which func(this, val) == true
    ObjectList collect(bool (*func)(Object const*, void const*), void const*) const;
    
    /// number of objects for which func(this, val) == true
    size_t count(bool (*func)(Object const*, void const*), void const*) const;
    
    /// print a summary of the content (nb of objects, class)
    void report(std::ostream&) const;
    
    /// sum tension accross plane defined by `n.pos + a = 0`
    void infoTension(size_t& cnt, real& sum, real& inf, real& sup, Vector const& n, real a) const;

    //--------------------------

    /// distribute the Couple on the fibers to approximate an equilibrated state
    void equilibrateSym(FiberSet const&, CoupleList&, CoupleProp const*);

    /// distribute Couples of given class on the fibers to approximate an equilibrated state
    void equilibrate(FiberSet const&, CoupleList&, CoupleProp const*);
    
    /// distribute all Couple on the fibers to approximate an equilibrated state
    void equilibrate(FiberSet const&, PropertyList const&);
    
    /// distribute all Couple on the fibers to approximate an equilibrated state
    void equilibrate();

    /// distribute Couples on filament intersections
    void bindToIntersections(FiberSet const&, CoupleList&, PropertyList const&);
    
    /// distribute Couples on filament intersections
    void bindToIntersections(CoupleProp const*);
    
    /// distribute all free Couples on filament intersections
    void bindToIntersections();

    /// bring all objects to centered image using periodic boundary conditions
    void foldPositions(Modulo const*) const;

    //--------------------------

    /// unlink all objects before import
    void freeze();

    /// delete objects that were not updated during import
    void prune();
    
    ///debug function
    int bad() const;
};


#endif

