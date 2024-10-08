// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University
// doubly linked list, STL style, with acces by iterators,
// some additions to manipulate the list: sorting, unsorting, etc.

#include "object.h"
#include "object_pool.h"
#include "assert_macro.h"
#include "random.h"
#include <stdlib.h>


void ObjectPool::push_front(Object * n)
{
    //std::clog << this << " push_front " << n->reference() << "\n";
    assert_true(n->set_);
    n->prev(nullptr);
    n->next(frontO);
    if ( frontO )
    {
        frontO->prev(n);
        assert_true(backO->next()==nullptr);
    }
    else
        backO = n;
    frontO = n;
    ++nSize;
}


void ObjectPool::push_back(Object * n)
{
    //std::clog << "ObjectPool: push_back " << n->reference() << "\n";
    n->prev(backO);
    n->next(nullptr);
    if ( backO )
    {
        assert_true(frontO->prev()==nullptr);
        backO->next(n);
    }
    else
        frontO = n;
    backO = n;
    ++nSize;
}


void ObjectPool::push_before(Object * p, Object * n)
{
    n->next(p);
    n->prev(p->prev());
    if ( p->prev() )
        p->prev()->next(n);
    else
        frontO = n;
    p->prev(n);
    ++nSize;
}


void ObjectPool::push_after(Object * p, Object * n)
{
    n->prev(p);
    n->next(p->next());
    if ( p->next() )
        p->next()->prev(n);
    else
        backO = n;
    p->next(n);
    ++nSize;
}


/**
 Transfer obj and any following objects the end of `this`
 */
void ObjectPool::grab(ObjectPool& list, Object * obj)
{
    if ( backO )
        backO->next(obj);
    else
        frontO = obj;
    Object * tmp = obj->prev();
    if ( tmp )
        tmp->next(nullptr);
    else
        list.frontO = nullptr;
    obj->prev(backO);
    backO = list.backO;
    list.backO = tmp;
    while ( obj )
    {
        ++nSize;
        --list.nSize;
        obj = obj->next();
    }
    //assert_false(bad());
    //assert_false(list.bad());
}


void ObjectPool::grab(ObjectPool& list)
{
    Object * n = list.frontO;
    
    if ( n )
    {
        if ( backO )
            backO->next(n);
        else
            frontO = n;
        
        n->prev(backO);
        backO = list.backO;
        nSize += list.nSize;
        
        list.nSize  = 0;
        list.frontO = nullptr;
        list.backO  = nullptr;
    }
    //assert_false(bad());
}


void ObjectPool::pop_front()
{
    --nSize;
    Object * n = frontO;
    frontO = frontO->next();
    n->next(nullptr);  // unnecessary?

    if ( frontO )
    {
        frontO->prev(nullptr);
        assert_true(backO->next()==nullptr);
    }
    else
        backO = nullptr;
}


void ObjectPool::pop_back()
{
    --nSize;
    Object * n = backO;
    backO = backO->prev();
    n->prev(nullptr);  // unnecessary?

    if ( backO )
    {
        backO->next(nullptr);
        assert_true(frontO->prev()==nullptr);
    }
    else
        frontO = nullptr;
}


void ObjectPool::pop(Object * n)
{
    assert_true( nSize > 0 );
    Object * x = n->next();

    if ( n->prev() )
        n->prev()->next(x);
    else {
        assert_true( frontO == n );
        frontO = x;
    }
    
    if ( x )
        x->prev(n->prev());
    else {
        assert_true( backO == n );
        backO = n->prev();
    }
    
    n->prev(nullptr); // unnecessary?
    n->next(nullptr); // unnecessary?
    --nSize;
}


size_t ObjectPool::truncate(Object * obj)
{
    size_t cnt = 0;
    Object * n = frontO;
    Object * p = nullptr;
    while ( n && n != obj )
    {
        ++cnt;
        p = n;
        n = n->next();
    }
    if ( n == obj )
    {
        nSize = cnt;
        backO = p;
        if ( p )
            p->next(nullptr);
        else
            frontO = nullptr;
    }
    return cnt;
}


void ObjectPool::clear()
{
#if ( 0 )
    // thorough unnecessary cleanup?
    Object * p, * n = frontO;
    while ( n )
    {
        n->prev(nullptr);
        p = n->next();
        n->next(nullptr);
        n = p;
    }
#endif
    frontO = nullptr;
    backO  = nullptr;
    nSize  = 0;
}


void ObjectPool::erase()
{
    Object * n = frontO;
    Object * p;
    while ( n )
    {
        p = n->next();
        delete(n);
        n = p;
    }
    frontO = nullptr;
    backO  = nullptr;
    nSize  = 0;
}


//------------------------------------------------------------------------------
#pragma mark - Shuffle


/**
 Rearrange [F--Q][P--B] into [P--B][F--Q]
 */
void ObjectPool::permute(Object * p)
{
    if ( p != frontO )
    {
        // close list into a loop
        backO->next(frontO);
        frontO->prev(backO);
        
        // open loop at 'p'
        frontO = p;
        backO  = p->prev();
        
        backO->next(nullptr);
        frontO->prev(nullptr);
    }
    //assert_false(bad());
}


/**
 Rearrange [F--P][X--Y][Q--B] into [X--Y][F--P][Q--B]
 
 Q must be after P: Q > P
 If Q is between Front and P, this will destroy the list,
 but it would be costly to check this condition here.
 */
void ObjectPool::shuffle_up(Object * p, Object * q)
{
    assert_true( p != q );
    assert_true( p  &&  p->next() );
    assert_true( q  &&  q->prev() );
    
    if ( q != p->next() )
    {
        frontO->prev(q->prev());
        q->prev()->next(frontO);
        frontO = p->next();
        frontO->prev(nullptr);
        p->next(q);
        q->prev(p);
    }
    //assert_false(bad());
}


/**
 Rearrange [F--P][X--Y][Q--B] into [F--P][Q--B][X--Y]
 
 Q must be after P
 If Q is between Front and P, this will destroy the list,
 but it would be costly to check this condition here.
 */
void ObjectPool::shuffle_down(Object * p, Object * q)
{
    assert_true( p  &&  p->next() );
    assert_true( q  &&  q->prev() );
    
    if ( q != p->next() )
    {
        backO->next(p->next());
        p->next()->prev(backO);
        p->next(q);
        backO = q->prev();
        backO->next(nullptr);
        q->prev(p);
    }
    //assert_false(bad());
}


/**
 This could be improved, as we traverse the list from the root.
 We could instead pick a random node, using Inventory, and move from there
 Upon hitting the list end, one would restart, knowning the order of the nodes
 */
void ObjectPool::shuffle()
{
    size_t pp, qq;
    if ( nSize > UINT32_MAX )
    {
        pp = RNG.pint64(nSize);
        qq = RNG.pint64(nSize);
    }
    else
    {
        pp = RNG.pint32((uint32_t)nSize);
        qq = RNG.pint32((uint32_t)nSize);
    }

    size_t n = 0;
    Object *p = frontO, *q;

    if ( pp+1 < qq )
    {
        for ( ; n < pp; ++n )
            p = p->next();
        for ( q = p; n < qq; ++n )
            q = q->next();
        
        shuffle_up(p, q);
    }
    else if ( qq+1 < pp )
    {
        for ( ; n < qq; ++n )
            p = p->next();
        for ( q = p; n < pp; ++n )
            q = q->next();
        
        shuffle_down(p, q);
    }
    else
    {
        for ( ; n < qq; ++n )
            p = p->next();

        permute(p);
    }
}

/**
 This traverses a quarter of the list on average, starting from a random
 node in the list that is provided externally.
 */
void ObjectPool::shuffle(Object * p)
{
    assert_true(p);
    size_t i = RNG.pint32((uint32_t)nSize>>2);
    
    Object * q = p;
    while ( q->next() && --i > 0 )
        q = q->next();
    
    if ( p != q )
        shuffle_up(p, q);
}

//------------------------------------------------------------------------------
#pragma mark - Count

size_t ObjectPool::count() const
{
    size_t cnt = 0;
    Object * p = front();
    while ( p )
    {
        ++cnt;
        p = p->next();
    }
    return cnt;
}


size_t ObjectPool::count(Object const* n) const
{
    size_t res = 0;
    Object * p = front();
    while ( p )
    {
        res += ( p == n );
        p = p->next();
    }
    return res;
}


size_t ObjectPool::count(bool (*func)(Object const*, void const*), void const* arg) const
{
    size_t res = 0;
    Object const* n = front();
    while ( n )
    {
        res += func(n, arg);
        n = n->next();
    }
    return res;
}

//------------------------------------------------------------------------------
#pragma mark - Check

/**
 This traverses the entire list and checks every link, which is costly
 */
int ObjectPool::bad() const
{
    size_t cnt = 0;
    
    if ( frontO && frontO->prev() != nullptr )
        return 1;
    
    Object * p = frontO, * q;
    while ( p )
    {
        q = p->next();
        if ( !q )
        {
            if ( p != backO )
                return 2;
        }
        else
        {
            if ( q->prev() != p )
                return 4;
        }
        p = q;
        ++cnt;
    }
    
    if ( cnt != nSize )
    {
        std::clog << "ObjectPool believed " << nSize << " true " << cnt << "\n";
        return 8;
    }
    return 0;
}


