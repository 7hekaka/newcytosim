// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef BUDDY_H
#define BUDDY_H

#include <vector>
#include <algorithm>


/// Maintains a list of mutual relationship between objects.
/**
 Buddy implements mutual relationship between objects.
 
 The class keeps track of a list of `buddies`.
 Relationship is established with hello().
 When an object is destroyed, it calls goodbye() for all its buddies.
 
 This class can be used when an object needs to know if another object is destroyed,
 and vice-versa.
 
 F. Nedelec 11.08.2012 -- 16.04.2020
 */
class Buddy
{

private:
    
    /// type for a list of buddies
    typedef std::vector<Buddy *> BuddyList;
    
    /// list of buddies
    BuddyList buddies_;
    
private:
    
    /// replace the buddy that may have been at index `ix`
    void enlist(Buddy * b, size_t ix)
    {
        if ( ix < buddies_.size() )
        {
            if ( buddies_[ix]  &&  buddies_[ix] != b )
            {
                goodbye(buddies_[ix]);
                buddies_[ix]->goodbye(this);
                buddies_[ix]->unlist(this);
            }
        }
        else
            buddies_.resize(ix+1, nullptr);
        
        buddies_[ix] = b;
    }
    
    /// add `b` into the list of buddies, or complain if already present
    void enlist(Buddy * b)
    {
#if ( 1 )
        // complain if buddy is known already:
        BuddyList::iterator i = std::find(buddies_.begin(), buddies_.end(), b);
        if ( i != buddies_.end() )
        {
            std::clog << " Warning: duplicate Buddy::enlist()\n";
            return;
        }
#endif
        
        // find an empty spot:
        i = std::find(buddies_.begin(), buddies_.end(), nullptr);
        if ( i != buddies_.end() )
            *i = b;
        else
            buddies_.push_back(b);
        //std::clog << this << " has " << buddies_.size() << " buddies\n";
    }
    
    /// removes `b` from the list of known buddy, do not call goodbye()
    Buddy * unlist(Buddy * b)
    {
        BuddyList::iterator i = std::find(buddies_.begin(), buddies_.end(), b);
        if ( i != buddies_.end() )
        {
            *i = nullptr;
            return b;
        }
        else
            std::clog << " Warning: Buddy::unlist(unlisted)\n";
        return nullptr;
    }

public:
    
    /// constructor
    Buddy() {}
    
    /// upon destruction, invoke `goodbye` for all buddies
    virtual ~Buddy()
    {
        //std::clog << this << " ::~Buddy()\n";
        for ( Buddy * b : buddies_ )
        {
            if ( b )
            {
                b->goodbye(this);
                b->unlist(this);
            }
        }
    }
    
    /// this is called everytime a known buddy is destroyed
    virtual void goodbye(Buddy const* b)
    {
        //std::clog << "Buddy " << this << "::goodbye(" << b << ")\n";
    }
    
    /// used as a signal from a buddy
    virtual void salute(Buddy const*)
    {
    }

    /// invoke `salute(this)` for all buddies
    void salute()
    {
        for ( Buddy * b : buddies_ )
            b->salute(this);
    }
    
    /// will make `this` and `guy` mutual buddies
    void connect(Buddy * guy)
    {
        if ( guy )
        {
            enlist(guy);
            guy->enlist(this);
        }
    }
    
    /// remove `this` and `guy` from each other lists, without calling goodbye()
    void disconnect(Buddy * guy)
    {
        if ( guy )
        {
            unlist(guy);
            guy->unlist(this);
        }
    }

    /// returns the number of registered buddies
    size_t nbBuddies() const
    {
        return buddies_.size();
    }
    
    /// return buddy at index `ix`
    Buddy * buddy(const size_t ix) const
    {
        if ( ix < buddies_.size() )
            return buddies_[ix];
        return nullptr;
    }
    
    /// returns true if `guy` is a buddy
    bool check(Buddy const* guy) const
    {
        BuddyList::const_iterator i = std::find(buddies_.begin(), buddies_.end(), guy);
        
        return ( i != buddies_.end() );
    }
    
};

#endif

