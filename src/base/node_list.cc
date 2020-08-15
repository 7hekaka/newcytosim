// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// doubly linked list, STL style, with acces by iterators,
// some additions to manipulate the list: sorting, unsorting, etc.

#include "node_list.h"
#include "assert_macro.h"
#include "random.h"
#include <stdlib.h>


void NodeList::push_front(Node * n)
{
    //Cytosim::log("NodeList: push_front   %p in   %p\n", n, this);

    n->nPrev = nullptr;
    n->nNext = nFront;
    if ( nFront )
        nFront->nPrev = n;
    else
        nBack = n;
    nFront = n;
    ++nSize;
}


void NodeList::push_back(Node * n)
{
    //Cytosim::log("NodeList: push_back   %p in   %p\n", n, this);
    
    n->nPrev = nBack;
    n->nNext = nullptr;
    if ( nBack )
        nBack->nNext = n;
    else
        nFront = n;
    nBack = n;
    ++nSize;
}


/**
 Transfer objects in `list` to the end of `this`, until `list` is empty.
 */
void NodeList::merge(NodeList& list)
{
    Node * n = list.nFront;
    
    if ( n )
    {
        if ( nBack )
            nBack->nNext = n;
        else
            nFront = n;
        
        n->nPrev = nBack;
        nBack = list.nBack;
        nSize += list.nSize;
        
        list.nSize  = 0;
        list.nFront = nullptr;
        list.nBack  = nullptr;
    }
}


void NodeList::push_after(Node * p, Node * n)
{
    n->nPrev = p;
    n->nNext = p->nNext;
    if ( p->nNext )
        p->nNext->nPrev = n;
    else
        nBack = n;
    p->nNext = n;
    ++nSize;
}


void NodeList::push_before(Node * p, Node * n)
{
    n->nNext = p;
    n->nPrev = p->nPrev;
    if ( p->nPrev )
        p->nPrev->nNext = n;
    else
        nFront = n;
    p->nPrev = n;
    ++nSize;
}


void NodeList::pop_front()
{
    assert_true( nFront );
 
    Node * n = nFront->nNext;
    nFront = n;
    n->nNext = nullptr;  // unnecessary?

    if ( nFront )
        nFront->nPrev = nullptr;
    else
        nBack = nullptr;
    --nSize;
}


void NodeList::pop_back()
{
    assert_true( nBack );
    
    Node * n = nBack->nPrev;
    nBack = n;
    n->nPrev = nullptr;  // unnecessary?
    
    if ( nBack )
        nBack->nNext = nullptr;
    else
        nFront = nullptr;
    --nSize;
}


void NodeList::pop(Node * n)
{
    assert_true( nSize > 0 );
    Node * x = n->nNext;

    if ( n->nPrev )
        n->nPrev->nNext = x;
    else {
        assert_true( nFront == n );
        nFront = x;
    }
    
    if ( x )
        x->nPrev = n->nPrev;
    else {
        assert_true( nBack == n );
        nBack = n->nPrev;
    }
    
    n->nPrev = nullptr; // unnecessary?
    n->nNext = nullptr; // unnecessary?
    --nSize;
}


void NodeList::clear()
{
    Node * p, * n = nFront;
    while ( n )
    {
        n->nPrev = nullptr; // unnecessary?
        p = n->nNext;
        n->nNext = nullptr; // unnecessary?
        n = p;
    }
    nFront = nullptr;
    nBack  = nullptr;
    nSize  = 0;
}


void NodeList::erase()
{
    Node * n = nFront;
    Node * p;
    while ( n )
    {
        p = n->nNext;
        delete(n);
        n = p;
    }
    nFront = nullptr;
    nBack  = nullptr;
    nSize  = 0;
}



size_t NodeList::count() const
{
    size_t cnt = 0;
    Node * p = nFront;
    while ( p )
    {
        ++cnt;
        p = p->nNext;
    }
    return cnt;
}

//------------------------------------------------------------------------------
#pragma mark - Sort

/**
This is a bubble sort, which scales like O(N^2)
comp(a,b) = -1 if (a<b) and 1 if (a>b) or 0
*/
void NodeList::bubblesort(int (*comp)(const Node*, const Node*))
{
    Node * ii = front();
    
    if ( ii == nullptr )
        return;
    
    ii = ii->nNext;
    
    while ( ii )
    {
        Node * kk = ii->nNext;
        Node * jj = ii->nPrev;
        
        if ( comp(ii, jj) > 0 )
        {
            jj = jj->nPrev;
            
            while ( jj && comp(ii, jj) > 0 )
                jj = jj->nPrev;
            
            pop(ii);
            
            if ( jj )
                push_after(jj, ii);
            else
                push_front(ii);
        }
        ii = kk;
    }
}


/**
The merge sort should have O(N*log(N))
comp(a,b) = -1 if (a<b) and 1 if (a>b) or 0
*/
void NodeList::mergesort(int (*comp)(const Node*, const Node*))
{
    Node * n = Node::mergesort(comp, nFront);
    nFront = n;
    while ( n->nNext )
    {
        n = n->nNext;
    }
    nBack = n;
}


void NodeList::move_behind(Node *& subF, Node *& subL, Node * pvt, Node * down)
{
    if ( down != subL )
        down->nNext->nPrev = down->nPrev;
    else {
        if ( down == nBack )
            nBack = down->nPrev;
        else
            down->nNext->nPrev = down->nPrev;
        subL = down->nPrev;
    }
    down->nPrev->nNext = down->nNext;
    down->nNext = pvt;
    down->nPrev = pvt->nPrev;
    pvt->nPrev = down;
    if ( pvt != subF )
        down->nPrev->nNext = down;
    else {
        if ( pvt == nFront )
            nFront = down;
        else
            down->nPrev->nNext = down;
        subF = down;
    }
}

/*
 From `Partition Algorithms for the Doubly Linked List`
 John M. Boyer, University of Southern Mississippi; ACM 1990
 */
void NodeList::blinksort(int (*comp)(const Node*, const Node*), Node * subF, Node * subL)
{
    Node * pvt2 = subF->nNext;
    if ( comp(subF, pvt2) > 0 )
        move_behind(subF, subL, subF, pvt2);
    Node * pvt1 = subF;
    while ((pvt2 != subL) and comp(pvt2->nNext, pvt2) >= 0 )
        pvt2 = pvt2->nNext;
    if ( pvt2 != subL )
    {
        Node * down = subL;
        while ( down != pvt2 )
        {
            Node * temp = down->nPrev;
            if ( comp(pvt1, down) > 0 )
                move_behind(subF, subL, pvt1, down);
            else if ( comp(pvt2, down) > 0 )
                move_behind(subF, subL, pvt2, down);
            down = temp;
        }
        if ((pvt1 != subF) and (pvt1->nPrev != subF))
            blinksort(comp, subF, pvt1->nPrev);
        if ((pvt1->nNext != pvt2) and (pvt1->nNext != pvt2->nPrev))
            blinksort(comp, pvt1->nNext, pvt2->nPrev);
        if ((pvt2 != subL) and (pvt2->nNext != subL))
            blinksort(comp, pvt2->nNext, subL);
    }
}


void NodeList::blinksort(int (*comp)(const Node*, const Node*))
{
    if ( nFront != nBack )
        blinksort(comp, nFront, nBack);
}


/**
This copies the data to a temporary space to use the standard library qsort()
comp(Node** a, Node** b) = -1 if (a<b) and 1 if (a>b) or 0
*/
void NodeList::quicksort(int (*comp)(const void*, const void*))
{
    const size_t cnt = nSize;
    Node ** tmp = new Node*[cnt];
    
    size_t i = 0;
    Node * n = nFront;
    
    while( n )
    {
        tmp[i++] = n;
        n = n->nNext;
    }
    
    qsort(tmp, cnt, sizeof(Node*), comp);
    
    n = tmp[0];
    nFront = n;
    n->nPrev = nullptr;
    for( i = 1; i < nSize; ++i )
    {
        n->nNext = tmp[i];
        tmp[i]->nPrev = n;
        n = tmp[i];
    }
    n->nNext = nullptr;
    nBack = n;
    
    delete[] tmp;
}


//------------------------------------------------------------------------------
#pragma mark - Shuffle


/**
 Rearrange [F--P][Q--L] into [Q--L][F--P]
 */
void NodeList::permute(Node * p)
{
    if ( p  &&  p->nNext )
    {
        nBack->nNext   = nFront;
        nFront->nPrev  = nBack;
        nFront         = p->nNext;
        nBack          = p;
        nBack->nNext   = nullptr;
        nFront->nPrev  = nullptr;
    }
    assert_false( bad() );
}


/**
 Rearrange [F--P][X--Y][Q--L] into [X--Y][F--P][Q--L]
 
 If Q is between nFront and P, this will destroy the list,
 but there is no way to check such condition here.
 */
void NodeList::shuffle_up(Node * p, Node * q)
{
    assert_true( p  &&  p->nNext );
    assert_true( q  &&  q->nPrev );
    
    if ( q != p->nNext )
    {
        nFront->nPrev   = q->nPrev;
        q->nPrev->nNext = nFront;
        nFront          = p->nNext;
        nFront->nPrev   = nullptr;
        p->nNext        = q;
        q->nPrev        = p;
    }
    assert_false( bad() );
}


/**
 Rearrange [F--P][X--Y][Q--L] into [F--P][Q--L][X--Y]
 
 If Q is between nFront and P, this will destroy the list,
 but there is no way to check such condition here.
 */
void NodeList::shuffle_down(Node * p, Node * q)
{
    assert_true( p  &&  p->nNext );
    assert_true( q  &&  q->nPrev );
    
    if ( q != p->nNext )
    {
        nBack->nNext    = p->nNext;
        p->nNext->nPrev = nBack;
        p->nNext        = q;
        nBack           = q->nPrev;
        nBack->nNext    = nullptr;
        q->nPrev        = p;
    }
    assert_false( bad() );
}


void NodeList::shuffle()
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
    Node *p = nFront, *q;

    if ( pp+1 < qq )
    {
        for ( ; n < pp; ++n )
            p = p->nNext;
        for ( q = p; n < qq; ++n )
            q = q->nNext;
        
        shuffle_up(p, q);
    }
    else if ( qq+1 < pp )
    {
        for ( ; n < qq; ++n )
            p = p->nNext;
        for ( q = p; n < pp; ++n )
            q = q->nNext;
        
        shuffle_down(p, q);
    }
    else
    {
        for ( ; n < qq; ++n )
            p = p->nNext;

        permute(p);
    }
}


void NodeList::shuffle3()
{
    shuffle();
    shuffle();
    shuffle();
}


//------------------------------------------------------------------------------
#pragma mark - Check


bool NodeList::check(Node const* n) const
{
    Node * p = nFront;
    while ( p )
    {
        if ( p == n )
            return true;
        p = p->nNext;
    }
    return false;
}


int NodeList::bad() const
{
    size_t cnt = 0;
    Node * p = nFront, * q;
    
    if ( p  &&  p->nPrev != nullptr )
        return 1;
    while ( p )
    {
        q = p->nNext;
        if ( q == nullptr ) {
            if ( p != nBack )
                return 2;
        }
        else {
            if ( q->nPrev != p )
                return 3;
        }
        p = q;
        ++cnt;
    }
    
    if ( cnt != nSize )
        return 4;
    return 0;
}


