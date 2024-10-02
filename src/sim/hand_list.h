// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#ifndef HAND_LIST_H
#define HAND_LIST_H

#include "real.h"


class Hand;
class Fiber;


/// a list of Hands, used to keep track of Hands attached to a Fiber
class HandList
{
    
    /// Pointer to first attached Hand
    Hand * haFront;
    
    /// Pointer to last attached Hand
    Hand * haBack;
    
public:

    /// constructor
    HandList() : haFront(nullptr), haBack(nullptr) {}
    
    /// first in list
    Hand * front() { return haFront; }
    
    /// last in list
    Hand * back() { return haBack; }
    
    /// register a new Hands that attached to this Fiber
    void add(Hand*);
    
    /// unregister bound Hands (which has detached)
    void remove(Hand*);
    
    /// update all Hands bound to this
    void updateAll() const;
    
    /// detach all Hands
    void detachAll();
    
    /// sort Hands by order of increasing abscissa
    void sort();
    
    /// number of attached Hands
    size_t count() const;
    
    /// number of times Hand is registered (should be zero or 1)
    size_t count(Hand const*) const;

    /// a function to count Hands using a custom criteria
    long count(int (*func)(Hand const*)) const;
    
    /// count Couple attached to another Fiber
    size_t countLinks(Fiber const*) const;

    /// number of Hands attached within a range of abscissa
    size_t countInRange(real abs_min, real abs_max) const;
    
    /// check data integrity
    int bad(Fiber const*) const;
};


#endif

