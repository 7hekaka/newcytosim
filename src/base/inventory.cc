// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University.

#include "inventory.h"
#include "assert_macro.h"
#include "exceptions.h"
#include <iostream>
#include <limits>


Inventory::Inventory()
{
    alloca_ = 31;
    record_ = new Inventoried*[1+alloca_];
    lowest_ = ~0;
    highest_ = 0;

    for ( ObjectID n = 0; n <= alloca_; ++n )
        record_[n] = nullptr;
}


void Inventory::allocate(size_t sz)
{
    constexpr size_t chunk = 1024;
    sz = ( sz + chunk ) & ~( chunk -1 );
    
    Inventoried ** ptr = new Inventoried*[1+sz];
    
    ObjectID n = 0;
    for ( ; n <= alloca_; ++n )
        ptr[n] = record_[n];
    for ( ; n <= sz; ++n )
        ptr[n] = nullptr;
    
    delete[] record_;
    record_ = ptr;
    alloca_ = sz;
    //std::clog << "Inventory::allocated(" << sz << ")\n";
}


void Inventory::release()
{
    delete[] record_;
    record_ = nullptr;
    alloca_ = 0;
    lowest_ = ~0;
    highest_ = 0;
}


void Inventory::clear()
{
    for ( ObjectID n = lowest_; n <= highest_; ++n )
        record_[n] = nullptr;
    lowest_ = ~0;
    highest_ = 0;
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 This will assign a new serial-number for `obj`, if it does not have one.
 */
void Inventory::assign(Inventoried * obj)
{
    ObjectID n = obj->identity();
    
    if ( n <= 0 )
    {
        n = ++highest_;
        obj->setIdentity(n);
        //std::clog << "identity(" << obj << ") <- " << n << "\n";
    }
    else
    {
        // already allocated and registered
        if ( n < alloca_ && record_[n] == obj )
            return;
        highest_ = std::max(highest_, n);
    }
    
    if ( n >= alloca_ )
        allocate(n);
    
    assert_true(!record_[n]);
    record_[n] = obj;
    //std::clog << "identity(" << obj << ") = " << n << "\n";

    lowest_ = std::min(lowest_, n);
    //assert_false(bad());
}


void Inventory::unassign(const Inventoried * obj)
{
    ObjectID n = obj->identity();
    assert_true( n <= highest_ );
    assert_true( n == 0 || record_[n] == obj );
    record_[n] = nullptr;
    
    if ( n == lowest_ )
    {
        if ( ++lowest_ <= highest_ )
        {
            assert_true(record_[highest_]);
            while ( !record_[lowest_] )
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
        assert_true(record_[lowest_]);
        while ( !record_[highest_] )
            --highest_;
    }
}


/**
 This will pack the array by filling up the empty spots,
 by moving objects from the end of the list to fill in available slots.
 */
void Inventory::reassign()
{
    ObjectID nxt = first_unassigned();
    ObjectID inf = nxt;
    ObjectID sup = highest_;
    assert_true(!record_[inf]);
    assert_true(record_[sup]);
    while ( inf < sup )
    {
        ++inf;
        if ( record_[inf] )
        {
            while ( record_[nxt] )
                ++nxt;
            assert_true( nxt < alloca_ );
            // move record from `inf` to `nxt`:
            record_[nxt] = record_[inf];
            record_[nxt]->setIdentity(nxt);
            record_[inf] = nullptr;
            ++nxt;
        }
    }
    lowest_ = 1;
    highest_ = nxt-1;
}

//------------------------------------------------------------------------------
#pragma mark -


ObjectID Inventory::first_unassigned() const
{
    assert_true(!record_[alloca_]);
    ObjectID n = 1;
    while ( record_[n] )
        ++n;
    assert_true(n <= alloca_);
    return n;
}


ObjectID Inventory::next_identity(ObjectID n) const
{
    ++n;
    assert_true(n <= alloca_);
    if ( record_[n] )
        return n;
    if ( ++n <= highest_ )
    {
        while ( !record_[n] )
            ++n;
        return n;
    }
    return 0;
}


Inventoried * Inventory::get(const ObjectID n) const
{
    if ( n <= highest_ )
    {
        assert_true(!record_[n] || record_[n]->identity()==n );
        return record_[n];
    }
    return nullptr;
}


Inventoried* Inventory::first() const
{
    if ( lowest_ < alloca_ )
        return record_[lowest_];
    return nullptr;
}


Inventoried* Inventory::last() const
{
    assert_true(highest_ < alloca_);
    assert_true(0==highest_ || record_[highest_]);
    return record_[highest_];
}


Inventoried* Inventory::next(Inventoried const* i) const
{
    ObjectID n = i->identity() + 1;
    assert_true(n <= alloca_);
    if ( record_[n] )
        return record_[n];
    if ( ++n <= highest_ )
    {
        while ( !record_[n] )
            ++n;
        assert_true(n <= alloca_);
        return record_[n];
    }
    return nullptr;
}


Inventoried* Inventory::previous(Inventoried const* i) const
{
    ObjectID n = i->identity() - 1;
    if ( lowest_ <= n )
    {
        assert_true(record_[lowest_]);
        while ( !record_[n] )
            --n;
        return record_[n];
    }
    return nullptr;
}


//------------------------------------------------------------------------------
#pragma mark - I/O


size_t Inventory::count() const
{
    size_t cnt = 0;
    for ( ObjectID n = 0; n <= highest_; ++n )
        if ( record_[n] ) ++cnt;
    return cnt;
}


void Inventory::print(std::ostream& os) const
{
    os << "Inventory " << this << "\n";
    for ( ObjectID n = lowest_; n <= highest_; ++n )
        os << n << " -> " << record_[n] << "\n";
}


std::ostream& operator << (std::ostream& os, Inventory const& arg)
{
    arg.print(os);
    return os;
}


int Inventory::bad() const
{
    if ( lowest_ < 1 ) return 1;
    for ( ObjectID n = 0; n < lowest_ && n < alloca_; ++n )
    {
        if ( record_[n] ) return 2;
    }
    for ( ObjectID n = lowest_; n <= highest_; ++n )
    {
        Inventoried * i = record_[n];
        if ( i && n != i->identity() ) return 4;
    }
    for ( ObjectID n = highest_+1; n <= alloca_; ++n )
    {
        if ( record_[n] ) return 8;
    }
    return 0;
}
