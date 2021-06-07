// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "inventory.h"
#include "assert_macro.h"
#include "exceptions.h"
#include <limits>


Inventory::Inventory()
{
    allocated_ = 31;
    byNames    = new Inventoried*[1+allocated_];
    lowest_    = std::numeric_limits<ObjectID>::max();
    highest_   = 0;

    for ( ObjectID n = 0; n <= allocated_; ++n )
        byNames[n] = nullptr;
}


void Inventory::allocate(size_t sz)
{
    constexpr size_t chunk = 32;
    sz = ( sz + chunk ) & ~( chunk -1 );
    
    Inventoried ** ptr = new Inventoried*[sz];
    
    ObjectID n = 0;
    for ( ; n <= allocated_; ++n )
        ptr[n] = byNames[n];
    for ( ; n < sz; ++n )
        ptr[n] = nullptr;
    
    delete[] byNames;
    byNames    = ptr;
    allocated_ = sz-1;
}


void Inventory::release()
{
    delete[] byNames;
    byNames = nullptr;
    allocated_ = 0;
    lowest_    = std::numeric_limits<ObjectID>::max();
    highest_   = 0;
}


void Inventory::clear()
{
    for ( ObjectID n = lowest_; n <= highest_; ++n )
        byNames[n] = nullptr;
    lowest_ = std::numeric_limits<ObjectID>::max();
    highest_ = 0;
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 This will assign a new serial-number for `obj`, if it does not have one.
 */
void Inventory::assign(Inventoried * obj)
{
    ObjectID& n = obj->identity();
    
    if ( n == 0 )
        n = ++highest_;
    else
        highest_ = std::max(highest_, n);
    
    if ( n >= allocated_ )
        allocate(n+1);
    
    assert_true(!byNames[n]);
    byNames[n] = obj;
    
    lowest_ = std::min(lowest_, n);
    //std::clog << "Inventory::store() assigned " << n << " to " << obj << "\n";
}


void Inventory::unassign(const Inventoried * obj)
{
    ObjectID n = obj->identity();
    assert_true( n <= highest_ );
    assert_true( byNames[n] == obj );
    byNames[n] = nullptr;
    
    if ( n == lowest_ )
    {
        if ( ++lowest_ <= highest_ )
        {
            assert_true(byNames[highest_]);
            while ( !byNames[lowest_] )
                ++lowest_;
        }
        else
        {
            assert_true(n == highest_);
            lowest_ = ~0U; // max unsigned value
            highest_ = 0;
        }
    }
    else if ( n == highest_ )
    {
        --highest_;
        assert_true(lowest_ <= highest_);
        assert_true(byNames[lowest_]);
        while ( !byNames[highest_] )
            --highest_;
    }
}


/**
 This will pack the array by filling up the empty spots,
 without changing the order of the objects.
 */
void Inventory::reassign()
{
    ObjectID nxt = first_unassigned();
    ObjectID inf = nxt;
    ObjectID sup = last_identity();
    assert_true(!byNames[inf]);
    
    while ( inf <= sup )
    {
        assert_true(byNames[sup]);
        while ( !byNames[inf] )
            ++inf;
        //swap:
        assert_true(!byNames[nxt]);
        byNames[nxt] = byNames[inf];
        byNames[nxt]->identity(nxt);
        byNames[inf] = nullptr;
        ++nxt;
        ++inf;
    }
    
    lowest_ = 1;
    highest_ = nxt-1;
}

//------------------------------------------------------------------------------
#pragma mark -


ObjectID Inventory::first_unassigned() const
{
    assert_true(!byNames[allocated_]);
    ObjectID n = 1;
    while ( byNames[n] )
        ++n;
    assert_true(n <= allocated_);
    return n;
}


ObjectID Inventory::next_identity(ObjectID n) const
{
    assert_true(n < allocated_);
    if ( byNames[++n] )
        return n;
    if ( ++n <= highest_ )
    {
        while ( !byNames[n] )
            ++n;
        return n;
    }
    return 0;
}


Inventoried * Inventory::get(const size_t n) const
{
    if ( n <= highest_ )
    {
        assert_true(n < allocated_);
        assert_true(!byNames[n] || byNames[n]->identity()==n );
        return byNames[n];
    }
    return nullptr;
}


Inventoried* Inventory::first() const
{
    if ( lowest_ < allocated_ )
        return byNames[lowest_];
    return nullptr;
}


Inventoried* Inventory::last() const
{
    assert_true(highest_ < allocated_);
    assert_true(0==highest_ || byNames[highest_]);
    return byNames[highest_];
}


Inventoried* Inventory::next(Inventoried const* i) const
{
    ObjectID n = i->identity() + 1;
    assert_true(n <= allocated_);
    if ( byNames[n] )
        return byNames[n];
    if ( ++n <= highest_ )
    {
        while ( !byNames[n] )
            ++n;
        assert_true(n < allocated_);
        return byNames[n];
    }
    return nullptr;
}


Inventoried* Inventory::previous(Inventoried const* i) const
{
    ObjectID n = i->identity() - 1;
    if ( lowest_ <= n )
    {
        assert_true(byNames[lowest_]);
        while ( !byNames[n] )
            --n;
        return byNames[n];
    }
    return nullptr;
}


//------------------------------------------------------------------------------
#pragma mark -


size_t Inventory::count() const
{
    size_t cnt = 0;
    for ( ObjectID n = 0; n <= highest_; ++n )
        if ( byNames[n] ) ++cnt;
    return cnt;
}


void Inventory::print(std::ostream& os) const
{
    os << "Inventory " << this << "\n";
    for ( ObjectID n = lowest_; n <= highest_; ++n )
        os << n << " -> " << byNames[n] << "\n";
}


std::ostream& operator << (std::ostream& os, Inventory const& arg)
{
    arg.print(os);
    return os;
}

