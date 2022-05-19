// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#ifndef OBJECT_SET_H
#define OBJECT_SET_H

#include "object.h"
#include "object_pool.h"
#include "inventory.h"
#include "property.h"

class Outputter;
class PropertyList;
class Glossary;
class Simul;

/// A set of Object
/**
 Encapsulates the different functions used to manage Objects.
 Pointers to the Objects are stored in two lists:
 - a doubly linked list: pool
 - an array: inventory
 .
 The ObjectPool pool is mixed at every time step,
 and thus it can be used to access the objects in a random order,
 as necessary for Monte-Carlo. 
 
 The Inventory can be used to access objects directly.
 
 Functions are used to manage:
 - object creation: newProperty(), newObjects().
 - object lists: size(), add(), remove(), erase().
 - object access: first(), find().
 - simulation: step(), shuffle().
 - I/O: loadObject(), read(), write(), freeze(), thaw().
 .
 */
class ObjectSet
{
private:
    
    ObjectSet();

public:

    /// holds pointers to the Objects organized by ObjectID
    Inventory inventory_;
    
    /// holds pointers to the Objects in a doubly linked list
    ObjectPool pool_;
    
    /// holds pointers to the Objects in a doubly linked list
    ObjectPool ice_;

    /// the Simul containing this set
    Simul& simul_;
    
protected:
    
    /// mark all objects from given list with value `f`
    static void flag(ObjectPool const&, ObjectFlag f);
    
    /// collect objects from ObjectPool for which func(obj, val) == true
    static size_t count(ObjectPool const&, bool (*func)(Object const*, void const*), void const*);

    /// collect all objects
    static ObjectList collect(ObjectPool const&);

    /// collect objects from ObjectPool for which func(obj, val) == true
    static ObjectList collect(ObjectPool const&, bool (*func)(Object const*, void const*), void const*);

    /// write Object in ObjectPool to file
    static void writeObjects(Outputter&, ObjectPool const&);
    
    /// print a list of the content (nb of objects, class)
    void writeReport(std::ostream&, const std::string& title) const;
    
    /// delete  Objects from sub list
    static void erasePool(ObjectPool&);

public:
    
    /// unlink all objects before import
    virtual void freeze();
    
    /// delete objects that were not updated during import
    virtual void prune();
    
    /// relink all objects after import
    void thaw();
    
    /// apply translation to all Objects in ObjectList
    static void translateObjects(ObjectList const&, Vector const&);
    
    /// apply rotation to all Objects in ObjectList
    static void rotateObjects(ObjectList const&, Rotation const&);
    
    /// apply Isometry to all Objects in ObjectList
    static void moveObjects(ObjectList const&, Isometry const&);

    /// set flag of all Objects to `f`
    static void flagObjects(ObjectList const&, ObjectFlag f);

    /// apply translation to unflagged Objects in list
    static void translateObjects(ObjectList const&, Vector const&, ObjectFlag f);
    
    /// apply rotation to unflagged Objects in list
    static void rotateObjects(ObjectList const&, Rotation const&, ObjectFlag f);

    /// apply Isometry to unflagged Objects in list
    static void moveObjects(ObjectList const&, Isometry const&, ObjectFlag f);

public:
    
    /// creator
    ObjectSet(Simul& s) : simul_(s) { }
    
    /// destructor
    virtual ~ObjectSet() { erase(); }    
    
    //--------------------------

    /// create a new property of category `cat` for a class `name`
    virtual Property * newProperty(const std::string& cat, const std::string& name, Glossary&) const = 0;
    
    /// create Objects of class `name`, given the options provided in `opt`
    virtual void newObjects(ObjectList&, const std::string& name, Glossary& opt) = 0;
    
    /// create Object with given Tag and PropertyID (used for reading trajectory file)
    virtual Object * newObject(ObjectTag, PropertyID) = 0;

    //--------------------------
    
    /// link the object last in the list
    virtual void link(Object *);
    
    /// link the object last in the list
    virtual void unlink(Object *);
    
    /// remove Object
    virtual void remove(Object *);

    /// register Object, adding it at the end of the list
    void add(Object *);
    
    /// remove Object, and delete it
    void erase(Object *);

    /// add multiple Objects
    void add(ObjectList const&);

    /// remove all Objects in list
    void remove(ObjectList const&);

    /// delete all Objects in list and forget all serial numbers
    virtual void erase();
    
    /// number of elements
    virtual size_t size()       const { return pool_.size(); }

    /// mix the order of elements in the doubly linked list pool
    virtual void shuffle()            { pool_.shuffle(); }
    
    /// first Object in the list
    Object * first()            const { return static_cast<Object*>(pool_.front()); }
    
    /// last Object
    Object * last()             const { return static_cast<Object*>(pool_.back()); }
    
    /// find Object of given serial-number (see Inventory)
    Object * findID(ObjectID n) const { return static_cast<Object*>(inventory_.get(n)); }
    
    /// return an Object which has this property
    Object * pickObject(Property const*) const;

    /// return Object corresponding to specifications
    Object * findObject(const std::string& cat, std::string spec, long identity) const;
    
    /// return Object corresponding to a certain criteria (eg. 'first' or 'last')
    Object * findObject(const std::string& cat, std::string spec) const;
    
    //--------------------------
    
    /// number of objects for which ( func(obj, val) == true )
    virtual size_t count(bool (*func)(Object const*, void const*), void const*) const;

    /// collect all objects
    virtual ObjectList collect() const;
 
    /// collect objects for which ( func(obj, val) == true )
    virtual ObjectList collect(bool (*func)(Object const*, void const*), void const*) const;

    /// collect objects that have given Property
    ObjectList collect(Property const*) const;
    
    /// load one Object from file, or skip it if `skip==true`
    void loadObject(Inputter&, ObjectTag tag, bool fat, bool update);
    
    /// write all Objects to file
    virtual void write(Outputter&) const = 0;
    
    /// print a summary of the content (nb of objects, class)
    virtual void report(std::ostream&) const = 0;

};


// This is declared here rather than in object.cc to permit inlining
inline Simul & Object::simul() const { return set_->simul_; }

#endif
