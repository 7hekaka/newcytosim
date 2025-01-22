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
 - simulation: steps(), shuffle().
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

    /// collect all objects
    static ObjectList collect(ObjectPool const&);

    /// collect objects from ObjectPool for which func(obj, val) == true
    static ObjectList collect(ObjectPool const&, bool (*func)(Object const*, void const*), void const*);
    
    /// print summary: nb of objects by class
    void writeReport(std::ostream&, const std::string& title) const;

    /// record number of objects, and topmost identity
    void writeRecords(Outputter&, size_t tot, size_t sup) const;
    
    /// write Object in ObjectPool to file
    void writePool(Outputter&, ObjectPool const&) const;

    /// delete  Objects from sub list
    static void erasePool(ObjectPool&);

public:
    
    /// prepare all objects for reading                                   <
    virtual void freeze();
    
    /// delete objects that are still on 'ice' because they were not imported
    void defrost();
    
    /// relink objects that are still on 'ice' because they were not imported
    void thaw();
    
    /// set flag of all Objects to `f`
    static void flagObjects(ObjectList const&, ObjectFlag f);

    /// set mark of all Objects to `k`
    static void markObjects(ObjectList const&, ObjectMark k);

    /// apply translation to all Objects in ObjectList
    static void translateObjects(ObjectList const&, Vector const&);
    
    /// apply rotation to all Objects in ObjectList
    static void rotateObjects(ObjectList const&, Rotation const&);
    
    /// apply Isometry to all Objects in ObjectList
    static void moveObject(Object*, Isometry const&);

    /// apply Isometry to all Objects in ObjectList
    static void moveObjects(ObjectList const&, Isometry const&);

    /// apply translation to Objects with `flag() != f`
    static void translateObjects(ObjectList const&, Vector const&, ObjectFlag f);
    
    /// apply rotation to Objects with `flag() != f`
    static void rotateObjects(ObjectList const&, Rotation const&, ObjectFlag f);

    /// apply Isometry to Objects with `flag() != f`
    static void moveObjects(ObjectList const&, Isometry const&, ObjectFlag f);

public:
    
    /// creator
    ObjectSet(Simul& s) : simul_(s) { }
    
    /// destructor
    virtual ~ObjectSet() { erase(); }    
    
    //--------------------------

    /// create a new property of category `cat` for a class `name`
    virtual Property * newProperty(std::string const& cat, std::string const& name, Glossary&) const = 0;
    
    /// create Objects of class `name`, given the options provided in `opt`
    virtual ObjectList newObjects(Property const*, Glossary& opt) = 0;
    
    /// create Object with given Tag and PropertyID (used for reading trajectory file)
    virtual Object * newObject(ObjectTag, PropertyID) = 0;

    //--------------------------
    
    /// link the object last in the list
    virtual void link(Object *);
    
    /// unlink object from the pool
    virtual void unlink(Object *);
    
    /// remove Object, from the pool and the inventory
    virtual void remove(Object *);
    
    /// remove Object, and delete it
    void eraseObject(Object *);

    /// allocate to be ready to handle indentities up to `sup`
    void reserve(size_t sup) { inventory_.reserve(sup); }
    
    /// register Object, adding it at the end of the list
    void add(Object *);

    /// add multiple Objects
    void add(ObjectList const&);

    /// remove Objects from list
    void remove(ObjectList const&);

    /// delete Objects from list
    void eraseObjects(ObjectList const&);

    /// delete all Objects in list and forget all serial numbers
    virtual void erase();
    
    /// number of registered elements
    virtual size_t size() const { return pool_.size(); }

    /// mix the order of elements in the doubly linked list pool
    virtual void shuffle() { pool_.shuffle(); }
    
    /// first Object in the list
    Object * first() const { return static_cast<Object*>(pool_.front()); }
    
    /// last Object
    Object * last() const { return static_cast<Object*>(pool_.back()); }
    
    /// find Object of given serial-number (see Inventory)
    Object * identifyObject(ObjectID n) const { return static_cast<Object*>(inventory_.get(n)); }
    
    /// check if object's identity match the inventory record (for debugging)
    bool badIdentity(Object const* obj) const { return obj != inventory_.get(obj->identity()); }
    
    /// name by which given Space can be recovered
    std::string nameObject(Object const*) const;

    /// return Object corresponding to specifications
    Object * findObject(Property const*, long identity) const;

    /// return Object corresponding to specifications
    Object * findObject(const std::string& cat, std::string spec, long identity) const;
    
    /// return an Object which has this property
    Object * pickObject(Property const*) const;

    /// return Object corresponding to a certain criteria (eg. 'first' or 'last')
    Object * pickObject(const std::string& cat, std::string spec) const;
    
    //--------------------------
    
    /// number of objects for which ( func(obj, val) == true )
    virtual size_t count(bool (*func)(Object const*, void const*), void const*) const;
    
    /// number of objects for which ( property() == p )
    size_t count(Property const* p) const { return count(match_property, p); }

    /// collect all objects
    virtual ObjectList collect() const;

    /// collect all objects for which ( func(obj, val) == true )
    virtual ObjectList collect(bool (*func)(Object const*, void const*), void const*) const;
    
    /// collect `cnt` objects selected randomly
    ObjectList collect(size_t cnt) const;

    /// collect `cnt` objects for which ( func(obj, val) == true )
    ObjectList collect(bool (*func)(Object const*, void const*), void const*, size_t cnt) const;

    /// collect objects that have given Property
    ObjectList collect(Property const* p) const { return collect(match_property, p); }

    /// load one Object from file
    void loadObject(Inputter&, ObjectTag tag, int bin);
    
    /// read meta data from file
    virtual void readObjectTypes(Inputter&) {}

    /// print a summary of the content (nb of objects, class)
    virtual void report(std::ostream&) const = 0;

};


// This is declared here rather than in object.cc to permit inlining
inline Simul & Object::simul() const { return set_->simul_; }

#endif
