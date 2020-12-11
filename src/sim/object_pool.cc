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
    //Cytosim::log("ObjectPool: push_front   %p in   %p\n", n, this);

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
    //Cytosim::log("ObjectPool: push_back   %p in   %p\n", n, this);
    
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
void ObjectPool::merge(ObjectPool& list)
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


void ObjectPool::push_after(Object * p, Object * n)
{
    n->prevO = p;
    n->nextO = p->nextO;
    if ( p->nextO )
        p->nextO->prevO = n;
    else
        backO = n;
    p->nextO = n;
    ++nSize;
}


void ObjectPool::push_before(Object * p, Object * n)
{
    n->nextO = p;
    n->prevO = p->prevO;
    if ( p->prevO )
        p->prevO->nextO = n;
    else
        frontO = n;
    p->prevO = n;
    ++nSize;
}


void ObjectPool::pop_front()
{
    assert_true( frontO );
 
    Object * n = frontO->nextO;
    frontO = n;
    n->nextO = nullptr;  // unnecessary?

    if ( frontO )
        frontO->prevO = nullptr;
    else
        backO = nullptr;
    --nSize;
}


void ObjectPool::pop_back()
{
    assert_true( backO );
    
    Object * n = backO->prevO;
    backO = n;
    n->prevO = nullptr;  // unnecessary?
    
    if ( backO )
        backO->nextO = nullptr;
    else
        frontO = nullptr;
    --nSize;
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
    Object * p, * n = frontO;
    while ( n )
    {
        n->prevO = nullptr; // unnecessary?
        p = n->nextO;
        n->nextO = nullptr; // unnecessary?
        n = p;
    }
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



size_t ObjectPool::count() const
{
    size_t cnt = 0;
    Object * p = frontO;
    while ( p )
    {
        ++cnt;
        p = p->nextO;
    }
    return cnt;
}

//------------------------------------------------------------------------------
#pragma mark - Sort


Object* ObjectPool::split(Object* head)
{
    Object* slow = head;
    Object* fast = head;
    
    while ( fast->nextO && fast->nextO->nextO )
    {
        fast = fast->nextO->nextO;
        slow = slow->nextO;
    }
    
    Object *temp = slow->nextO;
    slow->nextO = nullptr;
    return temp;
}


Object* ObjectPool::merge(int (*comp)(const Object*, const Object*), Object *first, Object *second)
{
    // If first linked list is empty
    if ( !first )
        return second;
    
    // If second linked list is empty
    if ( !second )
        return first;
    
    // Pick the smaller value
    if ( comp(second, first) > 0 )
    {
        first->nextO = merge(comp, first->nextO, second);
        first->nextO->prevO = first;
        first->prevO = nullptr;
        return first;
    }
    else
    {
        second->nextO = merge(comp, first, second->nextO);
        second->nextO->prevO = second;
        second->prevO = nullptr;
        return second;
    }
}


Object* ObjectPool::mergesort(int (*comp)(const Object*, const Object*), Object* head)
{
    if ( !head || !head->nextO )
        return head;
    
    Object* tail = split(head);
    
    // Recur for left and right halves
    head = mergesort(comp, head);
    tail = mergesort(comp, tail);
    
    // Merge the two sorted halves
    return merge(comp, head, tail);
}


/**
The merge sort should have O(N*log(N))
comp(a,b) = -1 if (a<b) and 1 if (a>b) or 0
*/
void ObjectPool::mergesort(int (*comp)(const Object*, const Object*))
{
    Object * n = mergesort(comp, frontO);
    frontO = n;
    while ( n->nextO )
        n = n->nextO;
    backO = n;
}


/**
This is a bubble sort, which scales like O(N^2)
comp(a,b) = -1 if (a<b) and 1 if (a>b) or 0
*/
void ObjectPool::bubblesort(int (*comp)(const Object*, const Object*))
{
    Object * ii = front();
    
    if ( ii == nullptr )
        return;
    
    ii = ii->nextO;
    
    while ( ii )
    {
        Object * kk = ii->nextO;
        Object * jj = ii->prevO;
        
        if ( comp(ii, jj) > 0 )
        {
            jj = jj->prevO;
            
            while ( jj && comp(ii, jj) > 0 )
                jj = jj->prevO;
            
            pop(ii);
            
            if ( jj )
                push_after(jj, ii);
            else
                push_front(ii);
        }
        ii = kk;
    }
}


void ObjectPool::move_behind(Object *& subF, Object *& subL, Object * pvt, Object * down)
{
    if ( down != subL )
        down->nextO->prevO = down->prevO;
    else {
        if ( down == backO )
            backO = down->prevO;
        else
            down->nextO->prevO = down->prevO;
        subL = down->prevO;
    }
    down->prevO->nextO = down->nextO;
    down->nextO = pvt;
    down->prevO = pvt->prevO;
    pvt->prevO = down;
    if ( pvt != subF )
        down->prevO->nextO = down;
    else {
        if ( pvt == frontO )
            frontO = down;
        else
            down->prevO->nextO = down;
        subF = down;
    }
}

/*
 From `Partition Algorithms for the Doubly Linked List`
 John M. Boyer, University of Southern Mississippi; ACM 1990
 */
void ObjectPool::blinksort(int (*comp)(const Object*, const Object*), Object * subF, Object * subL)
{
    Object * pvt2 = subF->nextO;
    if ( comp(subF, pvt2) > 0 )
        move_behind(subF, subL, subF, pvt2);
    Object * pvt1 = subF;
    while ((pvt2 != subL) and comp(pvt2->nextO, pvt2) >= 0 )
        pvt2 = pvt2->nextO;
    if ( pvt2 != subL )
    {
        Object * down = subL;
        while ( down != pvt2 )
        {
            Object * temp = down->prevO;
            if ( comp(pvt1, down) > 0 )
                move_behind(subF, subL, pvt1, down);
            else if ( comp(pvt2, down) > 0 )
                move_behind(subF, subL, pvt2, down);
            down = temp;
        }
        if ((pvt1 != subF) and (pvt1->prevO != subF))
            blinksort(comp, subF, pvt1->prevO);
        if ((pvt1->nextO != pvt2) and (pvt1->nextO != pvt2->prevO))
            blinksort(comp, pvt1->nextO, pvt2->prevO);
        if ((pvt2 != subL) and (pvt2->nextO != subL))
            blinksort(comp, pvt2->nextO, subL);
    }
}


void ObjectPool::blinksort(int (*comp)(const Object*, const Object*))
{
    if ( frontO != backO )
        blinksort(comp, frontO, backO);
}


/**
This copies the data to a temporary space to use the standard library qsort()
comp(Object** a, Object** b) = -1 if (a<b) and 1 if (a>b) or 0
*/
void ObjectPool::quicksort(int (*comp)(const void*, const void*))
{
    const size_t cnt = nSize;
    Object ** tmp = new Object*[cnt];
    
    size_t i = 0;
    Object * n = frontO;
    
    while( n )
    {
        tmp[i++] = n;
        n = n->nextO;
    }
    
    qsort(tmp, cnt, sizeof(Object*), comp);
    
    n = tmp[0];
    frontO = n;
    n->prevO = nullptr;
    for( i = 1; i < nSize; ++i )
    {
        n->nextO = tmp[i];
        tmp[i]->prevO = n;
        n = tmp[i];
    }
    n->nextO = nullptr;
    backO = n;
    
    delete[] tmp;
}


//------------------------------------------------------------------------------
#pragma mark - Shuffle


/**
 Rearrange [F--P][Q--L] into [Q--L][F--P]
 */
void ObjectPool::permute(Object * p)
{
    if ( p  &&  p->nextO )
    {
        backO->nextO  = frontO;
        frontO->prevO = backO;
        frontO        = p->nextO;
        backO         = p;
        backO->nextO  = nullptr;
        frontO->prevO = nullptr;
    }
    //assert_false( bad() );
}


/**
 Rearrange [F--P][X--Y][Q--L] into [X--Y][F--P][Q--L]
 
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
 Rearrange [F--P][X--Y][Q--L] into [F--P][Q--L][X--Y]
 
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
    if ( nSize < 2 )
        return;
    
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


//------------------------------------------------------------------------------
#pragma mark - Check


bool ObjectPool::count(Object const* n) const
{
    Object * p = frontO;
    while ( p )
    {
        if ( p == n )
            return 1;
        p = p->nextO;
    }
    return 0;
}

/**
 This traverses the entire list and checks every link, which is costly
 */
int ObjectPool::bad() const
{
    size_t cnt = 0;
    Object * p = frontO, * q;
    
    if ( p  &&  p->prevO != nullptr )
        return 1;
    while ( p )
    {
        q = p->nextO;
        if ( q == nullptr )
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


