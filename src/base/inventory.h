// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University.

#ifndef INVENTORY_H
#define INVENTORY_H

#include "inventoried.h"
#include "assert_macro.h"
#include <ostream>

/// Attributes and remember serial-numbers to Inventoried
/**
A Inventory assigns serial-numbers (of type ObjectID) to Inventoried, 
and it records a pointer to these objects.
 
Pointers to the objects can be recovered from their 'number' in constant time.

\author FJ Nedelec, August 2003.
*/
class Inventory
{
private:
    
    /// array of objects, stored at the index corresponding to their ObjectID
    /**
     This stores pointers to the objects, at index 'i' such that
         record_[i]->identity() == i
     Valid ObjectID are > 0, starting at 1
     */
    Inventoried ** record_;
    
    /// size of memory allocated
    size_t alloca_;
    
    /// lowest i > 0 for which `record_[i] != 0`
    ObjectID lowest_;
    
    /// highest i > 0 for which `record_[i] != 0`
    ObjectID highest_;
    
    /// memory allocation function
    void allocate(size_t size);
    
    /// memory allocation function
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
    ObjectID next_identity(ObjectID n) const;
    
    /// the smallest unassigned ID
    ObjectID first_unassigned() const;

    /// current size of array
    size_t capacity() const { return alloca_; }

    /// allocate to be ready to handle `cnt` references
    void reserve(const size_t cnt) { if ( cnt > alloca_ ) allocate(cnt); }
    
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
    Inventoried * operator[](ObjectID n) const { assert_true(n<alloca_); return record_[n]; }

    /// object with the smallest ID
    Inventoried * first() const;
    
    /// object with the largest ID
    Inventoried * last() const;
    
    /// return object just before given object
    Inventoried * previous(Inventoried const*) const;
    
    /// return object found just after given object
    Inventoried * next(Inventoried const*) const;
    
    /// Human friendly ouput
    void print(std::ostream&) const;
};


/// output of all values
std::ostream& operator << (std::ostream&, Inventory const&);


#endif
