// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University


#include "hand_list.h"
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
    //std::clog << this << " add " << n->prop->name() << '\n';
    n->prev(nullptr);
    n->next(haFront);
    if ( haFront )
        haFront->prev(n);
    else
        haBack = n;
    haFront = n;
}


void HandList::remove(Hand * n)
{
    //assert_true(count(n) == 1);
    //std::clog << this << " rem " << n->prop->name() << '\n';
    Hand * x = n->next();
    if ( n->prev() )
        n->prev()->next(x);
    else {
        assert_true( haFront == n );
        haFront = x;
    }
    
    if ( x )
        x->prev(n->prev());
    else {
        assert_true( haBack == n );
        haBack = n->prev();
    }
}


void HandList::update() const
{
    for ( Hand * h = haFront; h; h = h->next() )
        h->update();
}


void HandList::detachAll() const
{
    // we must iterate one step ahead, because detach() will unlink
    Hand * h = haFront;
    while ( h )
    {
        Hand * n = h->next();
        h->detach();
        h = n;
    }
}

/**
Sort in ascending order
*/
static int compareAbscissa(const void* A, const void* B)
{
    real a = (*static_cast<Hand *const*>(A))->abscissa();
    real b = (*static_cast<Hand *const*>(B))->abscissa();
    return ( a > b ) - ( b > a );
}

/**
 This sorts the Hands in order of increasing abscissa
 Sorting is done by copying to temporary array space, using std::qsort
 */
void HandList::sort()
{
    if ( haFront != haBack )
    {
        size_t cnt = count();
        Hand ** tmp = new Hand*[cnt];
        
        size_t i = 0;
        Hand * n = haFront;
        
        while ( n )
        {
            tmp[i++] = n;
            n = n->next();
        }
        
        qsort(tmp, cnt, sizeof(Hand*), compareAbscissa);
        
        n = tmp[0];
        haFront = n;
        n->prev(nullptr);
        for ( i = 1; i < cnt; ++i )
        {
            n->next(tmp[i]);
            tmp[i]->prev(n);
            n = tmp[i];
        }
        n->next(nullptr);
        haBack = n;
        
        delete[] tmp;
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


size_t HandList::count(int (*func)(Hand const*)) const
{
    size_t res = 0;
    
    for ( Hand const* h = haFront; h; h = h->next() )
        res += func(h);
    
    //printf("nbHands(%p) = %u\n", count, res);
    return res;
}


size_t HandList::countInRange(real i, real s) const
{
    size_t res = 0;
        
    for ( Hand const* h = haFront; h; h = h->next() )
        res += (( i <= h->abscissa()) & ( h->abscissa() <= s ));
    
    //printf("nbHandsInRange(%8.2f, %8.2f) = %u\n", a, b, res);
    return res;
}
