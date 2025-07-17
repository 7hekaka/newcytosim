// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University.

#include "inventory.h"
#include "assert_macro.h"
#include "exceptions.h"
#include <iostream>
#include <limits>


Inventory::Inventory()
{
    alloca_ = 30;
    record_ = new Inventoried*[2+alloca_];
    lowest_ = ~0U;
    highest_ = 0U;

    for ( ObjectID n = 0; n <= alloca_; ++n )
        record_[n] = nullptr;
    // place sentinel:
    record_[1+alloca_] = nullptr;
}


// a sentinel `nullptr` is placed at the end of the array
void Inventory::allocate(size_t sup)
{
    constexpr size_t chunk = 1024;
    // we want sup+1 to be a multiple of chunk:
    sup = ( sup + 1 ) | ( chunk - 1 );
    //for ( size_t s = 0; s < 32; ++s ) printf("%lu %lu\n", s, 1+((s+1)|(chunk-1)));
    
    assert_true(sup > alloca_);
    //std::clog << this << "::allocated(" << alloca_ << " --> " << 1+sup << ")\n";
    Inventoried ** ptr = new Inventoried*[1+sup];
    
    ObjectID n = 0;
    for ( ; n <= alloca_; ++n )
        ptr[n] = record_[n];
    for ( ; n <= sup; ++n )
        ptr[n] = nullptr;

    delete[] record_;
    record_ = ptr;
    alloca_ = sup-1; // max Identity that can be handled
}


void Inventory::release()
{
    delete[] record_;
    record_ = nullptr;
    alloca_ = 0;
    lowest_ = ~0U;
    highest_ = 0U;
}


void Inventory::clear()
{
    for ( ObjectID n = lowest_; n <= highest_; ++n )
        record_[n] = nullptr;
    for ( ObjectID n = highest_; n <= alloca_; ++n )
        record_[n] = nullptr;
    lowest_ = ~0U;
    highest_ = 0U;
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
        //std::clog << this << ":highest <- " << n << "\n";
        //std::clog << "identity(" << obj << ") <- " << n << "\n";
    }
    else
    {
        // already allocated and registered
        if ( n <= alloca_ && record_[n] == obj )
            return;
        highest_ = std::max(highest_, n);
    }
    assert_true(!record_[1+alloca_]);

    if ( n > alloca_ )
        allocate(n);
    
    assert_true(!record_[n]);

    record_[n] = obj;
    //std::clog << "  identity(" << obj << ") = " << n << " ";

    lowest_ = std::min(lowest_, n);
    //assert_false(bad());
}


//
void Inventory::unassign(const Inventoried * obj)
{
    ObjectID n = obj->identity();
    assert_true( n <= highest_ );
    assert_true( n == 0 || record_[n] == obj );
    record_[n] = nullptr;
    
    if ( n == lowest_ )
    {
        assert_true(highest_ <= alloca_);
        // find next valid identity:
        for ( ++n; n <= highest_; ++n )
        {
            if ( record_[n] )
            {
                lowest_ = n;
                return;
            }
        }
        lowest_ = ~0U; // max unsigned value
    }
}


/**
 This will pack the array by filling up the empty spots,
 by moving objects from the end of the list to fill in available slots.
 */
void Inventory::reassign()
{
    ObjectID nxt = firstUnassigned();
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
            assert_true( nxt <= alloca_ );
            // move record from `inf` to `nxt`:
            record_[nxt] = record_[inf];
            record_[nxt]->setIdentity(nxt);
            record_[inf] = nullptr;
            ++nxt;
        }
    }
    lowest_ = 1U;
    highest_ = nxt-1;
}

//------------------------------------------------------------------------------
#pragma mark -


/** This relies on the sentinel `nullptr` at [1+alloca_] */
ObjectID Inventory::firstUnassigned() const
{
    assert_true(!record_[1+alloca_]);
    ObjectID n = 1;
    while ( record_[n] )
        ++n;
    assert_true(n <= 1+alloca_);
    if ( n == 1+alloca_ )
        return ~0U;
    return n;
}


ObjectID Inventory::nextIdentity(ObjectID n) const
{
    assert_true(n <= alloca_);
    assert_true(highest_ <= alloca_);
    for ( ++n; n <= highest_; ++n )
        if ( record_[n] )
            return n;
    return 0U;
}


Inventoried * Inventory::get(const ObjectID n) const
{
    if ( n <= highest_ )
    {
        Inventoried * rec = record_[n];
        assert_true(!rec || rec->identity()==n );
        return rec;
    }
    return nullptr;
}


Inventoried* Inventory::first() const
{
    if ( lowest_ <= alloca_ )
        return record_[lowest_];
    return nullptr;
}


Inventoried* Inventory::last() const
{
    assert_true(highest_ <= alloca_);
    assert_true(0==highest_ || record_[highest_]);
    return record_[highest_];
}


Inventoried* Inventory::next(ObjectID n) const
{
    assert_true(n <= alloca_);
    assert_true(highest_ <= alloca_);
    for ( ++n; n <= highest_; ++n )
        if ( record_[n] )
            return record_[n];
    return nullptr;
}


Inventoried* Inventory::previous(ObjectID n) const
{
    if ( lowest_ < n )
    {
        assert_true(record_[lowest_]);
        do
            --n;
        while ( !record_[n] );
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
    for ( ObjectID n = 0; n < lowest_ && n <= alloca_; ++n )
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
    if ( record_[1+alloca_] ) return 16;
    return 0;
}
