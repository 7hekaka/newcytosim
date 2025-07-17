// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University.

#ifndef INVENTORY_H
#define INVENTORY_H

#include "inventoried.h"
#include "assert_macro.h"
#include <ostream>

/// Attributes serial-numbers to Inventoried, and remember them in an list
/**
 The Inventory assigns serial-numbers (of type ObjectID) to Inventoried,
 and records a pointer to these objects in a table indexed by ObjectID.
 
 This permits pointers to the objects to be recovered from their 'ObjectID' in constant time.
 
 Note that a nullptr is placed at the end of the array, at [1+alloca_] as sentinel.

\author FJN, August 2003--2023.
*/
class Inventory
{
private:
    
    /// array of objects, stored at the index corresponding to their ObjectID
    /**
     This stores pointers to the objects, at index 'i' such that
         record_[i]->identity() == i
     Valid ObjectID are > 0, starting at 1, up to alloca_
     */
    Inventoried ** record_;
    
    /// size of memory allocated minus one
    size_t alloca_;
    
    /// lowest `i > 0` for which `record_[i] != nullptr`
    ObjectID lowest_;
    
    /// highest given identity
    ObjectID highest_;
    
    /// allocate memory to hold identities within [1, sup]
    void allocate(size_t sup);
    
    /// release memory
    void release();

    /// Disabled copy constructor
    Inventory(Inventory const&);
    
    /// Disabled copy assignment
    Inventory& operator = (Inventory const&);
    
public:
        
    /// Constructor
    Inventory();
    
    /// Destructor
    ~Inventory() { release(); }
    
    /// the smallest assigned ID
    ObjectID lowest() const { return lowest_; }

    /// the largest assigned ID
    ObjectID highest() const { return highest_; }
    
    /// lowest assigned ID strictly greater than `n`
    ObjectID nextIdentity(ObjectID n) const;
    
    /// the smallest unassigned ID, or max if all are assigned
    ObjectID firstUnassigned() const;

    /// current size of array
    size_t capacity() const { return alloca_; }

    /// allocate to handle ObjectID in [1, sup]
    void reserve(const size_t sup) { if ( sup > alloca_ ) allocate(sup); }
    
    /// remember `obj`, assign a new ID if necessary
    void assign(Inventoried * obj);
    
    /// forget the object and release its ID
    void unassign(const Inventoried * obj);
    
    /// reattribute all IDs consecutively, packing the array to remove any gap
    void reassign();

    /// number of non-zero IDs in the registry
    size_t count() const;

    /// clear all entries
    void clear();

    /// return the object with given ID or 0 if not found
    Inventoried * get(ObjectID) const;
    
    /// return object with given serial number
    Inventoried * operator[](ObjectID n) const { assert_true(n<=alloca_); return record_[n]; }

    /// object with the smallest ID
    Inventoried * first() const;
    
    /// object with the largest ID
    Inventoried * last() const;
    
    /// return object just before given object
    Inventoried * previous(ObjectID) const;
    
    /// return object just before given object
    Inventoried * previous(Inventoried const* obj) const { return previous(obj->identity()); }

    /// return object found just after given object
    Inventoried * next(ObjectID) const;
    
    /// return object found just after given object
    Inventoried * next(Inventoried const* obj) const { return next(obj->identity()); }

    /// Human friendly ouput
    void print(std::ostream&) const;
    
    /// debug checking
    int bad() const;
};


/// output of all values
std::ostream& operator << (std::ostream&, Inventory const&);


#endif
