// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef COUPLE_SET_H
#define COUPLE_SET_H

#include <vector>

#include "object_set.h"
#include "couple.h"
#include "couple_prop.h"

class PropertyList;


/// Set holding objects of class Couple
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
        switch (( attached2 << 1 ) | attached1 )
        {
            case 0: return ffList;
            case 1: return afList;
            case 2: return faList;
            default: return aaList;
        }
    }
    
    /// list of Properties for which `fast_diffusion == true`
    std::vector<CoupleProp const*> uniCouples;

    /// gather all Couple with `fast_diffusion` in reserve lists
    void uniStepCollect(Couple*);

    /// ensures that `can` holds `cnt` Couple, creating them of specified CoupleProp
    void uniRefill(CoupleProp const*, size_t cnt);

    /// attach Hand1 of Couple from `can` on locations specified in `loc`
    void uniAttach1(FiberSiteList& loc, CoupleStock& can);
    
    /// attach Hand2 of Couple from `can` on locations specified in `loc`
    void uniAttach2(FiberSiteList& loc, CoupleStock& can);
    
    /// attach both Hands of `nb` Couple at crossing points specified by arguments 1 & 2
    void uniAttach12(FiberSiteList&, FiberSiteList&, CoupleStock&, unsigned nb);

    /// `fast_diffusion` attachment assuming that free Couples are uniformly distributed
    void uniAttach(FiberSet const&);

    /// release Couples from reserve lists
    void uniRelax();
    
    /// save free Couple for which `save_unbound > 0`
    void writeSomeObjects(Outputter&) const;

public:

    /// constructor
    CoupleSet(Simul& s) : ObjectSet(s) {}
    
    //--------------------------
    
    /// identifies the set
    static std::string title() { return "couple"; }
    
    /// create a new property of category `cat` for a class `name`
    Property * newProperty(std::string const& cat, std::string const& name, Glossary&) const;
    
    /// create objects specified by Property, given options provided in `opt`
    ObjectList newObjects(Property const*, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object * newObject(ObjectTag, PropertyID);

    /// save objects
    void writeSet(Outputter&) const;

    //--------------------------

    /// add object (should be a Couple)
    void link(Object *);
    
    /// remove object (should be a Couple)
    void unlink(Object *);

    /// reassign Couple to sublist following attachement of Hand 1
    void relinkA1(Couple *);
    /// reassign Couple to sublist following detachment of Hand 1
    void relinkD1(Couple *);
    /// reassign Couple to sublist following attachement of Hand 2
    void relinkA2(Couple *);
    /// reassign Couple to sublist following detachment of Hand 2
    void relinkD2(Couple *);

    /// first unattached Couple
    Couple * firstFF() const { return static_cast<Couple*>(ffList.front()); }
    /// first Couple attached only by cHand1
    Couple * firstAF() const { return static_cast<Couple*>(afList.front()); }
    /// first Couple attached only by cHand2
    Couple * firstFA() const { return static_cast<Couple*>(faList.front()); }
    /// first Couple attached by both hands
    Couple * firstAA() const { return static_cast<Couple*>(aaList.front()); }

    /// last unattached Couple
    Couple * lastFF() const { return static_cast<Couple*>(ffList.back()); }
    /// last Couple attached by cHand1
    Couple * lastAF() const { return static_cast<Couple*>(afList.back()); }
    /// last Couple attached by cHand2
    Couple * lastFA() const { return static_cast<Couple*>(faList.back()); }
    /// last Couple attached by both hands
    Couple * lastAA() const { return static_cast<Couple*>(aaList.back()); }

    /// number of free Couples
    size_t sizeFF() const { return ffList.size(); }
    /// number of Couples attached by cHand1 only
    size_t sizeAF() const { return afList.size(); }
    /// number of Couples attached by cHand2 only
    size_t sizeFA() const { return faList.size(); }
    /// number of Couples attached by one hand
    size_t sizeA()  const { return faList.size() + afList.size(); }
    /// number of Couples attached by both hands
    size_t sizeAA() const { return aaList.size(); }
    /// total number of elements
    size_t size()   const { return ffList.size() + faList.size() + afList.size() + aaList.size(); }

    /// erase all objects
    void erase();
    
    /// detach all Hands
    void detachAll();

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
    
    
    /// return pointer to the Object of given ID, or zero if not found
    Couple * identifyObject(ObjectID n) const { return static_cast<Couple*>(inventory_.get(n)); }
    
    /// first Couple in inventory
    Couple * firstID() const { return static_cast<Couple*>(inventory_.first()); }
    
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

    /// initialize `fast_diffusion` attachment algorithm
    void uniPrepare(PropertyList const& properties);
    
    /// total count in reserves
    size_t countStocks() const;

    /// print number of elements in each reserve bin
    void infoStocks(std::ostream& os) const;

    //--------------------------

    /// distribute the Couple on the fibers to approximate an equilibrated state
    void equilibrateSym(FiberSet const&, CoupleProp const*, size_t);

    /// distribute Couples of given class on the fibers to approximate an equilibrated state
    void equilibrate(FiberSet const&, CoupleProp const*, CoupleStock&, size_t);
    
    /// distribute all Couple on the fibers to approximate an equilibrated state
    void equilibrate(FiberSet const&, CoupleProp const*);
    
    /// distribute all Couple on the fibers to approximate an equilibrated state
    void equilibrate();

    /// max `binding_range` of all known Couple
    real maxBindingRange() const;
    
    /// distribute given Couples on filament intersections
    void bindToIntersections(FiberSet const&, CoupleStock&, real range);
    
    /// distribute Couples of a certain type on filament intersections
    void bindToIntersections(CoupleProp const*);
    
    /// distribute all free Couples on filament intersections
    void bindToIntersections();

    /// bring all objects to centered image using periodic boundary conditions
    void foldPositions(Modulo const*) const;

    //--------------------------
    
    /// return a Couple from the reserve, or made by newCouple()
    Couple * makeCouple(CoupleProp const*);
    
    /// create a Couple at given position
    Couple * addFreeCouple(CoupleProp const*, Vector const&);

    /// create unattached Couples
    void makeCouples(CoupleProp const*, size_t cnt);

    /// create unattached Couples
    void makeCouples(size_t cnt[], PropertyID n_cnt);
    
    /// register given Couple, which must be unattached
    void addFreeCouple(Couple*);

    //--------------------------

    /// unlink all objects before import
    void freeze();
    
    /// move Couples into reserve lists, instead of deleting them
    void defrostStore();

    /// detach objects that were not updated during import
    void reheat(size_t cnt[], PropertyID n_cnt);
    
    /// detach objects that were not updated during import
    void reheat();

    ///debug function
    int bad() const;
};


#endif

