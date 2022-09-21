// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef BUDDY_H
#define BUDDY_H


#define MAKE_NO_FRIENDS 0

/// Establishes `circles of friends' within objects.
/**
 Buddy implements mutual relationship between objects.
 
 The class is used to keep track of a circle of `buddies` using a single-linked list.
 Relationship is established with connect().
 When an object is destroyed, it calls goodbye() for all its buddies.
 
 This class is used when an object needs to know if another object is destroyed.
 It also implements a virtual 'salute()' function that can be used
 to implement some basic form of communication between objects of the same circle.
 
 F. Nedelec 11.08.2012 -- 16.04.2020.
 21.09.2022: implemented circular list using one pointer.
 */
class Buddy
{

private:
    
    /// next buddy, making a circular single-linked list
    Buddy * buddy_;
    
public:
    
    /// search for 'guy' in circle of buddies, returning buddy anterior to `guy` if found
    Buddy * find(Buddy * guy)
    {
        Buddy * a = this;
        Buddy * b = buddy_;
        while ( b != this )
        {
            if ( b == guy )
                return a;
            a = b;
            b = b->buddy_;
        }
        return nullptr;
    }
    
    /// return number of buddies, in circle of friends
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

    /// add `b` into the list of buddies, or complain if already present
    void enlist(Buddy * guy)
    {
        assert_true( guy != this );
#if ( 1 )
        // complain if buddy is known already:
        if ( find(guy) )
        {
            std::clog << " Warning: duplicate Buddy::enlist()\n";
            return;
        }
#endif
        // assume new buddy is not linked in any other circle:
        assert_true( guy.buddy_ == guy );
        Buddy * b = buddy_;
        buddy_ = guy;
        guy->buddy_ = b;
        //std::clog << this << " has " << nbBuddies() << " buddies\n";
    }
    
    /// remove `b` from the list of known buddy, do not call goodbye()
    void unlist(Buddy * guy)
    {
        Buddy * b = find(guy);
        if ( b )
        {
            assert_true( b->buddy_ == guy );
            b->buddy_ = guy->buddy_;
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
    virtual void goodbye(Buddy const* b)
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
    
    /// return first buddy or *this
    Buddy const* buddy() const
    {
        if ( buddy_ != this )
            return buddy_;
        return nullptr;
    }
    
    /// print list of buddies
    void print(std::ostream& os) const
    {
        os << "Object " << this << " buddies are: ";
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

