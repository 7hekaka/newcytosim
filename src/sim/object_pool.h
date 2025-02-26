// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#ifndef OBJECT_POOL_H
#define OBJECT_POOL_H

#include <stddef.h>

class Object;

/// Doubly linked list of Objects
/**
 The ObjectPool is a doubly-linked list of objects. It holds pointers to the first
 and last elements of the list, and it keeps track of the number of objects linked.
 Functions are given to link and unlink Objects in constant time.\n

 This class is similar to the standard template library <std::list>
 and the naming of the functions is consistent with STL whenever possible.
 
 In addition, a function shuffle() randomizes the order of the Objects in the list, as
 necessary in a simulation to avoid any bias which could derive from fixed ordering.
 
 Function blinksort() can be used to sort the objects using a compare function.

 The list is null-terminated on both sides, and it can be traversed up or down:
 
     for ( Object * n = front(); n ; n = n->next() );
     for ( Object * n = back() ; n ; n = n->prev() ); 
 */

class ObjectPool
{

private:
    
    /// First Object in the list
    Object * frontO;
    
    /// Last Object in the list
    Object * backO;
    
    /// Number of Objects in the list
    size_t nSize;
    
    /// Disabled copy constructor
    ObjectPool(ObjectPool const&);
    
    /// Disabled copy assignment
    ObjectPool& operator = (ObjectPool const&);
    
public:
    
    /// default constructor
    ObjectPool() : frontO(nullptr), backO(nullptr), nSize(0) { }
    
    /// Destructor
    virtual ~ObjectPool()  { clear(); }
    
    /// First Object in list
    Object * front() const { return frontO; }

    /// Last Object in list
    Object * back() const { return backO; }
    
    /// pointer to first element
    Object * begin() const { return frontO; }
    
    /// pointer to a position just past the last element
    Object * end() const { return nullptr; }

    /// set first Object in list
    void front(Object * o) { frontO = o; }
    
    /// set last Object in list
    void back(Object * o) { backO = o; }

    /// Number of objects in the list
    size_t size() const { return nSize; }
    
    /// change number of objects in the list
    void size(size_t s) { nSize = s; }
 
    /// true if list has zero elements
    bool empty() const { return frontO == nullptr; }
    
    /// put Object first in the list
    void push_front(Object *);
    
    /// put Object last in the list
    void push_back(Object *);
    
    /// add `n` before `p`
    void push_before(Object * p, Object * n);

    /// add `n` after `p`
    void push_after(Object * p, Object * n);

    /// transfer object and following ones
    void grab(ObjectPool& list, Object*);
    
    /// import all objects from given list, emptying it
    void grab(ObjectPool& list);

    /// Remove Object `n` from list
    void pop(Object * n);
    
    /// Remove top Object from list, returning it
    void pop_front();
    
    /// Remove last Object from list
    void pop_back();
    
    /// remove objects located after `n`, including `n` itself
    size_t truncate(Object * n);

    /// clear the list
    void clear();
    
    /// delete all pool, clearing the list on the way
    void erase();
    
    /// sort according to given function
    void mergesort(int (*comp)(const Object*, const Object*))
    {
        if ( frontO != backO )
        {
            MergeSortJob<Object> job;
            job.sort(frontO, backO, comp);
        }
    }
    
    /// sort according to given function
    void blinksort(int (*comp)(const Object*, const Object*))
    {
        if ( frontO != backO )
        {
            std::clog << "blinksort " << size() << "\n";
            BlinkSortJob<Object> job;
            job.sort(frontO, backO, comp);
        }
    }

    /// Rearrange the list by exchanging the portions before and after `p`
    void permute(Object *);
    
    /// Rearrange the list by moving a central portion to the top
    void shuffle_up(Object *, Object *);
    
    /// Rearrange the list by moving a central portion to the bottom
    void shuffle_down(Object *, Object *);
    
    /// Mix list using permute() and shuffle() functions
    void shuffle();
    
    /// Mix list using shuffle() functions
    void shuffle(Object *);

    /// count number of elements in the list
    size_t count() const;
    
    /// returns 1 if element appears in the list
    size_t count(Object const* n) const;
    
    /// count objects from ObjectPool for which `func(obj, val) == true`
    size_t count(bool (*func)(Object const*, void const*), void const* val) const;

    /// test coherence of data structure
    int bad() const;

};


#endif
