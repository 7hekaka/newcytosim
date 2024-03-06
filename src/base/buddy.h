// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef BUDDY_H
#define BUDDY_H


/// Establishes `circles of friends' within objects.
/**
 Buddy implements mutual, symmetric, relationships between objects.
 
 The class is used to record a circle of `buddies` using a simple linked list.
 Relationship is established with connect().
 When an Object is destroyed, goodbye() if called for all its buddies.
 
 This class is used when an object needs to know if another object is destroyed,
 for example to inform Aster that one of its Fiber was destroyed.
 It also implements a virtual 'salute()' function that can be used to implement
 some basic form of communication between objects within the circle of buddies.
 
 F. Nedelec 11.08.2012 -- 16.04.2020.
 21.09.2022: implemented circular list using one pointer.
 28.02.2023: fixed unlist()
 */
class Buddy
{

private:
    
    /// next buddy, making a circular single-linked list
    Buddy * buddy_;
    
public:
    
    /// search for 'guy' in circle of buddies, returning Buddy anterior to `guy` if found
    Buddy const* findBuddy(Buddy const* guy) const
    {
        Buddy const* a = this;
        Buddy const* b = buddy_;
        while ( b != this )
        {
            if ( b == guy )
                return a;
            a = b;
            b = b->buddy_;
        }
        return nullptr;
    }
    
    /// true if 'guy' is within the circle of buddies
    bool isBuddy(Buddy const* guy) const
    {
        return this->findBuddy(guy);
    }
    
    /// return number of buddies in circle of friends
    size_t nbBuddies()
    {
        size_t cnt = 0;
        Buddy * b = buddy_;
        while ( b != this )
        {
            ++cnt;
            b = b->buddy_;
        }
        return cnt;
    }

    /// add `guy` into the list of buddies, joining circles if needed
    void enlist(Buddy * guy)
    {
        assert_true( guy != this );
        if ( !findBuddy(guy) )
        {
            // join the two circles:
            Buddy * b = buddy_;
            buddy_ = guy->buddy_;
            guy->buddy_ = b;
            //std::clog << this << " has " << nbBuddies() << " buddies\n";
        }
    }
    
    /// remove `guy` from the list of known buddy, do not call goodbye()
    void unlist(Buddy * guy)
    {
        Buddy * b = const_cast<Buddy*>(findBuddy(guy));
        if ( b )
        {
            assert_true( b->buddy_ == guy );
            b->buddy_ = guy->buddy_;
            guy->buddy_ = guy;
        }
        else
            std::clog << " Warning: Buddy::unlist(unlisted)\n";
    }
    
    /// removes `self` from the list of known buddy, do not call goodbye()
    void unlist()
    {
        Buddy * b = buddy_;
        while ( b->buddy_ != this )
            b = b->buddy_;
        b->buddy_ = buddy_;
        buddy_ = this;
    }

public:
    
    /// constructor
    Buddy() { buddy_ = this; }
    
    /// upon destruction, invoke `goodbye` for all buddies
    virtual ~Buddy()
    {
        // std::clog << this << " ::~Buddy()\n";
        Buddy * b = buddy_;
        while ( b != this )
        {
            b->goodbye(this);
            b = b->buddy_;
        }
        unlist();
    }
    
    /// this is called everytime a known buddy is destroyed
    virtual void goodbye(Buddy const*)
    {
        //std::clog << "Buddy " << this << "::goodbye(" << b << ")\n";
    }
    
    /// used as a signal between buddies
    virtual void salute(Buddy const*)
    {
    }

    /// invoke `salute(this)` for all buddies
    void salute()
    {
        Buddy * b = buddy_;
        while ( b != this )
        {
            b->salute(this);
            b = b->buddy_;
        }
    }
    
    /// will make `this` and `guy` mutual buddies
    void connect(Buddy * guy)
    {
        if ( guy )
            enlist(guy);
    }
    
    /// remove `this` and `guy` from each other lists, without calling goodbye()
    void disconnect(Buddy * guy)
    {
        if ( guy )
            unlist(guy);
    }
    
    /// print list of buddies
    void printBuddies(std::ostream& os) const
    {
        os << "Buddies of " << this << " : ";
        Buddy * b = buddy_;
        while ( b != this )
        {
            os << "   " << b;
            b = b->buddy_;
        }
        os << "\n";
    }
};

#endif

