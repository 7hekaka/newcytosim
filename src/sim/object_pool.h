// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#ifndef OBJECT_POOL_H
#define OBJECT_POOL_H

#include <stddef.h>
class Object;

/// Doubly linked list of Objects
/**
 
 This class is similar to the standard template library <std::list>
 and the naming of the functions is consistent with STL whenever possible.
 
 The ObjectPool holds pointers to the first and last elements of the list,
 and it keeps track of the number of objects linked.
 Functions are given to link and unlink Nodes in constant time.\n
 
 A function mix() randomize the order of the Nodes in the list, which is
 necessary in a simulation to avoid any bias which could derive from fixed ordering.
 
 The list is zero-terminated on both sides, and it can be traversed in either ways:
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
    ObjectPool& operator=(ObjectPool const&);
    
public:
    
    /// default constructor
    ObjectPool() : frontO(nullptr), backO(nullptr), nSize(0) { }
    
    /// Destructor
    virtual ~ObjectPool()  { clear(); }
    
    /// First Object in list
    Object * front() const { return frontO; }
    
    /// Last Object in list
    Object * back() const { return backO; }
    
    /// Number of objects in the list
    size_t size() const { return nSize; }
    
    /// true if list has zero elements
    bool empty() const { return frontO == nullptr; }
    
    /// put Object first in the list
    void push_front(Object *);
    
    /// put Object last in the list
    void push_back(Object *);
    
    /// import all objects from given list, emptying it
    void append(ObjectPool& list);
    
    /// link `n` after already linked `p`
    void push_after(Object * p, Object * n);
    
    /// link `n` before already linked `p`
    void push_before(Object * p, Object * n);
    
    /// Remove Object `n` from list
    void pop(Object * n);
    
    /// Remove top Object from list, returning it
    Object* pop_front();
    
    /// Remove last Object from list
    Object* pop_back();
    
    /// clear the list
    void clear();
    
    /// delete all pool, clearing the list on the way
    void erase();
    
    /// slowly sort according to given function
    void bubblesort(int (*comp)(const Object*, const Object*));
    
    /// sort according to given function
    void mergesort(int (*comp)(const Object*, const Object*));
    
    /// quicksort according to given function using std::qsort()
    void quicksort(int (*comp)(const void*, const void*));
    
    /// sort according to given function
    void blinksort(int (*comp)(const Object*, const Object*));
    
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
    bool count(Object const* n) const;
    
    /// test coherence of list
    int bad() const;
    
private:
    
    /// move sublist
    void move_behind(Object*&, Object*&, Object*, Object*);
    
    /// sort sublist
    void blinksort(int (*comp)(const Object*, const Object*), Object*, Object*);
    
    /// split linked list roughly in two
    static Object * split(Object*);
    
    /// merge two linked lists
    static Object * merge(int (*comp)(const Object*, const Object*), Object *first, Object *second);
    
    /// sort according to given function
    static Object * mergesort(int (*comp)(const Object*, const Object*), Object*);

};


#endif
