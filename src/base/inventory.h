// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

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
    
    /// array of objects, stored at the index corresponding to their ID
    /**
     This stores pointers to the objects, at index 'i' such that
         byNames[i]->identity() == i
     Valid ID are > 0.
     */
    Inventoried ** byNames;
    
    
    /// size of memory allocated
    size_t allocated_;
    
    /// lowest i > 0 for which `byNames[i] != 0`
    ObjectID lowest_;
    
    /// highest i > 0 for which `byNames[i] != 0`
    ObjectID highest_;
    
    /// memory allocation function
    void allocate(size_t size);
    
    /// memory allocation function
    void release();

    /// Disabled copy constructor
    Inventory(Inventory const&);
    
    /// Disabled copy assignment
    Inventory& operator=(Inventory const&);
    
public:
        
    /// Constructor
    Inventory();
    
    /// Destructor
    ~Inventory() { release(); }
    
    /// the smallest assigned identity
    ObjectID first_identity() const { return lowest_; }

    /// the largest assigned identity
    ObjectID last_identity() const { return highest_; }
    
    /// lowest assigned number strictly greater than `n`
    ObjectID next_identity(ObjectID n) const;
    
    /// the smallest unassigned number
    ObjectID first_unassigned() const;

    /// current size of array
    size_t capacity() const { return allocated_; }
    
    /// remember `obj`, assign a new ObjectID if necessary
    void assign(Inventoried * obj);
    
    /// forget the object and release its serial number
    void unassign(const Inventoried * obj);
    
    /// reattribute all serial numbers consecutively and pack the array
    void reassign();

    /// number of non-zero entries in the registry
    size_t count() const;

    /// clear all entries
    void clear();

    /// return the object with given serial number, or 0 if not found
    Inventoried * get(size_t number) const;
    
    /// return object with given serial number
    Inventoried * operator[](ObjectID n) const { assert_true(n<allocated_); return byNames[n]; }

    /// object with the smallest inventory number
    Inventoried * first() const;
    
    /// object with the largest inventory number
    Inventoried * last() const;
    
    /// return object just before in the inventory
    Inventoried * previous(Inventoried const*) const;
    
    /// return object just after in the inventory
    Inventoried * next(Inventoried const*) const;
    
    /// Human friendly ouput
    void print(std::ostream&) const;
};


/// output of all values
std::ostream& operator << (std::ostream&, Inventory const&);


#endif
