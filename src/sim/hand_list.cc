// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University


#include "hand_list.h"
#include "hand_monitor.h"
#include "hand.h"

/**
Link the hand at the front of the list.

Keeping track of the bound Hands is needed to run Cytosim if the filaments are
dynamic, so that a Fiber that has changed can update all the Hands bound to it.
The list is also used in reports, or to quickly count Hands bound to the fiber.
*/
void HandList::add(Hand * n)
{
    //assert_true(count(n) == 0);
    //std::clog << this << " add " << n->prop->name() << " " << n << '\n';
    n->prev(nullptr);
    n->next(haFront);
    if ( haFront )
        haFront->prev(n);
    else
        haBack = n;
    haFront = n;
    //assert_false(bad());
}


void HandList::remove(Hand * n)
{
    //assert_true(count(n) == 1);
    //std::clog << this << " rem " << n->prop->name() << " " << n << '\n';
    Hand * x = n->next();
    if ( n->prev() )
    {
        n->prev()->next(x);
    } else {
        assert_true( haFront == n );
        haFront = x;
    }
    
    if ( x ) {
        x->prev(n->prev());
    } else {
        assert_true( haBack == n );
        haBack = n->prev();
    }
    //assert_false(bad());
    n->prev(nullptr); // is this necessary?
    n->next(nullptr); // is this necessary?
}


void HandList::updateAll() const
{
    for ( Hand * h = haFront; h; h = h->next() )
        h->reinterpolate();
}


void HandList::detachAll()
{
    // we must iterate one step ahead, because detach() will unlink
    Hand * h = haFront;
    while ( h )
    {
        Hand * n = h->next();
        // update since upon detachment Hands may need their position?
        //h->reinterpolate();
        // no need to update Lattice here:
        h->Hand::detach();
        h = n;
    }
    assert_true(haFront==nullptr);
    assert_true(haBack==nullptr);
}


/// qsort function comparing Hand::abscissa()
static int compareHandAbscissa(const Hand* A, const Hand* B)
{
    real a = A->abscissa();
    real b = B->abscissa();
    return ( a > b ) - ( a < b );
}

/**
 This sorts the Hands in order of increasing abscissa
 Sorting is done by copying to temporary array space, using std::qsort
 */
void HandList::sort()
{
    if ( haFront != haBack )
    {
        BlinkSortJob<Hand> job;
        job.sort(haFront, haBack, compareHandAbscissa);
    }
}


size_t HandList::count() const
{
    size_t res = 0;
    
    for ( Hand const* h = haFront; h; h = h->next() )
        ++res;
    
    return res;
}


size_t HandList::count(Hand const* arg) const
{
    size_t res = 0;
    
    for ( Hand const* h = haFront; h; h = h->next() )
        res += ( h == arg );
    
    return res;
}


long HandList::count(int (*func)(Hand const*)) const
{
    long res = 0;
    
    for ( Hand const* h = haFront; h; h = h->next() )
        res += func(h);
    
    //printf("count(%p) = %u\n", this, res);
    return res;
}


size_t HandList::countLinks(Fiber const* fib) const
{
    size_t res = 0;
    
    for ( Hand const* h = haFront; h; h = h->next() )
    {
        Hand const* g = h->otherHand();
        if ( g )
            res += ( g->fiber() == fib );
    }
    //printf("countLinks(%p) = %u\n", fib, res);
    return res;
}


size_t HandList::countInRange(real i, real s) const
{
    size_t res = 0;
        
    for ( Hand const* h = haFront; h; h = h->next() )
        res += (( i <= h->abscissa()) && ( h->abscissa() <= s ));
    
    //printf("countInRange(%8.2f, %8.2f) = %u\n", i, s, res);
    return res;
}


/**
 This traverses the entire list, checking every link
 */
int HandList::bad(Fiber const* fib) const
{
    int res = 0;
    Hand * h = haFront;
    Hand * p = nullptr;
    while ( h )
    {
        if ( h->fiber() != fib )
            return 7;
        res |= ( h->prev() != p );
        p = h;
        h = h->next();
    }
    res |= 2 * ( p != haBack );
    return res;
}


