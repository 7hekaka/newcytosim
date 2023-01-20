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
    //std::clog << "ObjectPool: push_front " << n->reference() << "\n";
    assert_true(n->set_);
    n->prevO = nullptr;
    n->nextO = frontO;
    if ( frontO )
        frontO->prevO = n;
    else
        backO = n;
    frontO = n;
    ++nSize;
}


void ObjectPool::push_back(Object * n)
{
    //std::clog << "ObjectPool: push_back " << n->reference() << "\n";
    assert_true(n->set_);
    n->prevO = backO;
    n->nextO = nullptr;
    if ( backO )
        backO->nextO = n;
    else
        frontO = n;
    backO = n;
    ++nSize;
}


/**
 Transfer objects in `list` to the end of `this`, until `list` is empty.
 */
void ObjectPool::append(ObjectPool& list)
{
    Object * n = list.frontO;
    
    if ( n )
    {
        if ( backO )
            backO->nextO = n;
        else
            frontO = n;
        
        n->prevO = backO;
        backO = list.backO;
        nSize += list.nSize;
        
        list.nSize  = 0;
        list.frontO = nullptr;
        list.backO  = nullptr;
    }
}


Object* ObjectPool::pop_front()
{
    Object * n = frontO;
    if ( n )
    {
        --nSize;
        frontO = n->nextO;
        
        if ( frontO )
            frontO->prevO = nullptr;
        else
            backO = nullptr;
        
        n->nextO = nullptr;  // unnecessary?
    }
    return n;
}


Object* ObjectPool::pop_back()
{
    Object * n = backO;
    if ( n )
    {
        --nSize;
        backO = n->prevO;
    
        if ( backO )
            backO->nextO = nullptr;
        else
            frontO = nullptr;
        
        n->prevO = nullptr;  // unnecessary?
    }
    return n;
}


void ObjectPool::pop(Object * n)
{
    assert_true( nSize > 0 );
    Object * x = n->nextO;

    if ( n->prevO )
        n->prevO->nextO = x;
    else {
        assert_true( frontO == n );
        frontO = x;
    }
    
    if ( x )
        x->prevO = n->prevO;
    else {
        assert_true( backO == n );
        backO = n->prevO;
    }
    
    n->prevO = nullptr; // unnecessary?
    n->nextO = nullptr; // unnecessary?
    --nSize;
}


void ObjectPool::clear()
{
#if ( 0 )
    // thorough unnecessary cleanup?
    Object * p, * n = frontO;
    while ( n )
    {
        n->prevO = nullptr;
        p = n->nextO;
        n->nextO = nullptr;
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
        p = n->nextO;
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
        backO->nextO  = frontO;
        frontO->prevO = backO;
        
        // open loop at 'p'
        frontO = p;
        backO  = p->prevO;
        
        backO->nextO  = nullptr;
        frontO->prevO = nullptr;
    }
    //assert_false( bad() );
}


/**
 Rearrange [F--P][X--Y][Q--B] into [X--Y][F--P][Q--B]
 
 Q must be after P
 If Q is between Front and P, this will destroy the list,
 but it would be costly to check this condition here.
 */
void ObjectPool::shuffle_up(Object * p, Object * q)
{
    assert_true( p  &&  p->nextO );
    assert_true( q  &&  q->prevO );
    
    if ( q != p->nextO )
    {
        frontO->prevO   = q->prevO;
        q->prevO->nextO = frontO;
        frontO          = p->nextO;
        frontO->prevO   = nullptr;
        p->nextO        = q;
        q->prevO        = p;
    }
    //assert_false( bad() );
}


/**
 Rearrange [F--P][X--Y][Q--B] into [F--P][Q--B][X--Y]
 
 Q must be after P
 If Q is between Front and P, this will destroy the list,
 but it would be costly to check this condition here.
 */
void ObjectPool::shuffle_down(Object * p, Object * q)
{
    assert_true( p  &&  p->nextO );
    assert_true( q  &&  q->prevO );
    
    if ( q != p->nextO )
    {
        backO->nextO    = p->nextO;
        p->nextO->prevO = backO;
        p->nextO        = q;
        backO           = q->prevO;
        backO->nextO    = nullptr;
        q->prevO        = p;
    }
    //assert_false( bad() );
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
            p = p->nextO;
        for ( q = p; n < qq; ++n )
            q = q->nextO;
        
        shuffle_up(p, q);
    }
    else if ( qq+1 < pp )
    {
        for ( ; n < qq; ++n )
            p = p->nextO;
        for ( q = p; n < pp; ++n )
            q = q->nextO;
        
        shuffle_down(p, q);
    }
    else
    {
        for ( ; n < qq; ++n )
            p = p->nextO;

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
    size_t i = RNG.pint32((uint32_t)nSize>>1);
    
    Object * q;
    for ( q = p ; q && i > 0; --i )
        q = q->nextO;
    
    if ( q )
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
    
    if ( frontO && frontO->prevO != nullptr )
        return 1;
    
    Object * p = frontO, * q;
    while ( p )
    {
        q = p->nextO;
        if ( !q )
        {
            if ( p != backO )
                return 2;
        }
        else
        {
            if ( q->prevO != p )
                return 3;
        }
        p = q;
        ++cnt;
    }
    
    if ( cnt != nSize )
        return 4;
    return 0;
}


