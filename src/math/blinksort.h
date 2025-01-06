// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University.

#include "assert_macro.h"

/**
 Templated sorting methods adapted for doubly-linked lists, 19.01.2023, 5.01.2024
 Three methods are implemented:
 - bubblesort (slow)
 - mergesort
 - blinksort (best)
 */

//------------------------------------------------------------------------------
#pragma mark - Bubble sort

/**
This is a bubble sort, which scales like O(N^2)
To sort in increasing order, COMP(OBJECT * a, OBJECT * b) should return:
  -1 if (a<b)
   0 if (a=b)
  +1 if (a>b)
 .
*/
template<typename OBJECT>
class BubbleSortJob
{
private:

    void push_after(OBJECT *& back, OBJECT * p, OBJECT * n)
    {
        n->prev(p);
        n->next(p->next());
        if ( p->next() )
            p->next()->prev(n);
        else
            back = n;
        p->next(n);
    }
    
    void push_before(OBJECT *& front, OBJECT * p, OBJECT * n)
    {
        n->next(p);
        n->prev(p->prev());
        if ( p->prev() )
            p->prev()->next(n);
        else
            front = n;
        p->prev(n);
    }
    
public:
    
    void sort(OBJECT *& front, OBJECT *& back, int (*COMP)(const OBJECT *, const OBJECT *))
    {
        OBJECT * ii = front;
        
        if ( !ii )
            return;
        
        ii = ii->next();
        
        while ( ii )
        {
            OBJECT * kk = ii->next();
            OBJECT * jj = ii->prev();
            
            if ( COMP(ii, jj) > 0 )
            {
                jj = jj->prev();
                
                while ( jj && COMP(ii, jj) > 0 )
                    jj = jj->prev();
                
                //pop(ii):
                ii->prev_->next(ii->next());
                if ( ii->next() )
                    ii->next()->prev(ii->prev());
                else
                    back = ii->prev();
                
                if ( jj )
                    push_after(back, jj, ii);
                else {
                    ii->next(front);
                    ii->prev(nullptr);
                    front->prev(ii);
                    front = ii;
                }
            }
            ii = kk;
        }
    }
    
};

//------------------------------------------------------------------------------
#pragma mark - Merge sort

template<typename OBJECT>
class MergeSortJob
{
private:
    int (*COMP)(const OBJECT *, const OBJECT *);

    void splitlist_(OBJECT *& head, OBJECT *& tail)
    {
        assert_true( head && head != tail );
        while ( head != tail && head->next() != tail )
        {
            head = head->next();
            tail = tail->prev();
        }
        
        if ( head == tail )
            tail = head->next();
        
        head->next(nullptr);
        tail->prev(nullptr);
    }
    
    
    void mergelists_(OBJECT *& head, OBJECT * node, OBJECT * from, OBJECT *& tail)
    {
        assert_true( head && from );
        if ( COMP(head, from) < 0 )
        {
            OBJECT * temp = head->next();
            if ( temp )
            {
                // head remains on top of the list
                mergelists_(temp, node, from, tail);
                head->next(temp);
                temp->prev(head);
            }
            else
                head->next(from);
        }
        else
        {
            if ( from->next() )
                mergelists_(head, node, from->next(), tail);
            // add `from` on top of the list
            head->prev(from);
            from->next(head);
            from->prev(nullptr);
            head = from;
        }
    }
    
    void mergesort_(OBJECT *& head, OBJECT *& tail)
    {
        if ( head && head != tail )
        {
            OBJECT * node = head;
            OBJECT * next = tail;
            splitlist_(node, next);
            // Recur for left and right halves
            mergesort_(head, node);
            mergesort_(next, tail);
            // Merge the two sorted halves
            mergelists_(head, node, next, tail);
        }
    }

public:
    
    MergeSortJob() : COMP(nullptr) {}

    void sort(OBJECT *& head, OBJECT *& tail, int (*comp)(const OBJECT *, const OBJECT *))
    {
        COMP = comp;
        mergesort_(head, tail);
    }
    
};

//------------------------------------------------------------------------------
#pragma mark - Blink sort


/*
 From `Partition Algorithms for the Doubly Linked List`
 John M. Boyer, University of Southern Mississippi; ACM 1990
 Blink Sort is O(NlogN) for random data and 0(N) for reversed order and sorted lists.
 
 The comparison function must return an integer less than, equal to, or greater than zero
 if the first argument is considered to be respectively less than, equal to, or greater than the second.
 COMP(A, B) should return:
     - a negative integer if 'A < B'
     - a positive integer if 'B < A'
 for the list to be sorted in order of increasing values
 */
template<typename OBJECT>
class BlinkSortJob
{
private:
    
    OBJECT * FIRST;
    OBJECT * LAST;
    
    /// comparison function
    int (*COMP)(const OBJECT *, const OBJECT *);
    
    void movebehind_(OBJECT *& subF, OBJECT *& subL, OBJECT* pvt, OBJECT* down)
    {
        if ( down != subL )
            down->next_->prev_ = down->prev_;
        else {
            if ( down == LAST )
                LAST = down->prev_;
            else
                down->next_->prev_ = down->prev_;
            subL = down->prev_;
        }
        down->prev_->next_ = down->next_;
        down->next_ = pvt;
        down->prev_ = pvt->prev_;
        pvt->prev_ = down;
        if ( pvt != subF )
            down->prev_->next_ = down;
        else {
            if ( pvt == FIRST )
                FIRST = down;
            else
                down->prev_->next_ = down;
            subF = down;
        }
    }
    
    void blinksort_(OBJECT *& subF, OBJECT *& subL)
    {
        OBJECT * pvt2 = subF->next_;
        if ( COMP(subF, pvt2) > 0 )
            movebehind_(subF, subL, subF, pvt2);
        
        OBJECT * pvt1 = subF;
        while ((pvt2 != subL) and COMP(pvt2->next_, pvt2) >= 0 )
            pvt2 = pvt2->next_;
        
        if ( pvt2 != subL )
        {
            OBJECT * down = subL;
            while ( down != pvt2 )
            {
                OBJECT * temp = down->prev_;
                if ( COMP(down, pvt1) < 0 )
                    movebehind_(subF, subL, pvt1, down);
                else if ( COMP(down, pvt2) < 0 )
                    movebehind_(subF, subL, pvt2, down);
                down = temp;
            }
            if ((pvt1 != subF) and (pvt1->prev_ != subF))
                blinksort_(subF, pvt1->prev_);
            
            if ((pvt1->next_ != pvt2) and (pvt1->next_ != pvt2->prev_))
                blinksort_(pvt1->next_, pvt2->prev_);
            
            if ((pvt2 != subL) and (pvt2->next_ != subL))
                blinksort_(pvt2->next_, subL);
        }
    }
    
public:
    
    BlinkSortJob() : FIRST(nullptr), LAST(nullptr), COMP(nullptr) {}

    void sort(OBJECT *& front, OBJECT *& back, int (*comp)(const OBJECT *, const OBJECT *))
    {
        if ( front != back )
        {
            COMP = comp;
            FIRST = front;
            LAST = back;
            blinksort_(FIRST, LAST);
            front = FIRST;
            back = LAST;
        }
    }
};
