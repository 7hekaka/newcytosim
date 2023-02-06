// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University.
/**
 Templated sorting methods adapted for doubly-linked lists, 19.01.2023
 Three methods are implemented:
 - bubblesort (slow)
 - mergesort
 - blinksort (best)
 */

//------------------------------------------------------------------------------
#pragma mark - Bubble sort

template<typename OBJECT>
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


template<typename OBJECT>
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

/**
This is a bubble sort, which scales like O(N^2)
To sort in increasing order, comp(OBJECT* a, OBJECT* b) should return:
  -1 if (a<b)
   0 if (a=b)
  +1 if (a>b)
 .
*/
template<typename OBJECT>
void bubblesort(OBJECT* front, OBJECT* back, int (*comp)(const OBJECT*, const OBJECT*))
{
    OBJECT * ii = front;
    
    if ( !ii )
        return;
    
    ii = ii->next();
    
    while ( ii )
    {
        OBJECT * kk = ii->next();
        OBJECT * jj = ii->prev();
        
        if ( comp(ii, jj) > 0 )
        {
            jj = jj->prev();
            
            while ( jj && comp(ii, jj) > 0 )
                jj = jj->prev();
            
            pop(ii);
            
            if ( jj )
                push_after(back, jj, ii);
            else
                push_front(ii);
        }
        ii = kk;
    }
}

//------------------------------------------------------------------------------
#pragma mark - Merge sort

template<typename OBJECT>
void split_list(OBJECT*& head, OBJECT*& tail)
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


template<typename OBJECT>
void merge_lists(OBJECT*& head, OBJECT* node, OBJECT* from, OBJECT*& tail,
                    int (*comp)(const OBJECT*, const OBJECT*))
{
    assert_true( head && from );
    if ( comp(head, from) < 0 )
    {
        OBJECT * temp = head->next();
        if ( temp )
        {
            // head remains on top of the list
            merge_lists(temp, node, from, tail, comp);
            head->next(temp);
            temp->prev(head);
        }
        else
            head->next(from);
    }
    else
    {
        if ( from->next() )
            merge_lists(head, node, from->next(), tail, comp);
        // add `from` on top of the list
        head->prev(from);
        from->next(head);
        from->prev(nullptr);
        head = from;
    }
}


template<typename OBJECT>
void mergesort(OBJECT *& head, OBJECT *& tail,
                  int (*comp)(const OBJECT*, const OBJECT*))
{
    if ( head && head != tail )
    {
        OBJECT * node = head;
        OBJECT * next = tail;
        split_list(node, next);
        // Recur for left and right halves
        mergesort(head, node, comp);
        mergesort(next, tail, comp);
        // Merge the two sorted halves
        merge_lists(head, node, next, tail, comp);
    }
}

//------------------------------------------------------------------------------
#pragma mark - Blink sort

template<typename OBJECT>
void move_behind(OBJECT*& front, OBJECT*& back,
                 OBJECT*& subF, OBJECT*& subL, OBJECT* pvt, OBJECT* down)
{
    if ( down != subL )
        down->next()->prev(down->prev());
    else {
        if ( down == back )
            back = down->prev();
        else
            down->next()->prev(down->prev());
        subL = down->prev();
    }
    down->prev()->next(down->next());
    down->next(pvt);
    down->prev(pvt->prev());
    pvt->prev(down);
    if ( pvt != subF )
        down->prev()->next(down);
    else {
        if ( pvt == front )
            front = down;
        else
            down->prev()->next(down);
        subF = down;
    }
}

/*
 From `Partition Algorithms for the Doubly Linked List`
 John M. Boyer, University of Southern Mississippi; ACM 1990
 Blink Sort is O(NlogN) for random data and 0(N) for reversed order and sorted lists.
 */
template<typename OBJECT>
void blinksort(OBJECT*& front, OBJECT*& back,
               int (*comp)(const OBJECT*, const OBJECT*), OBJECT* subF, OBJECT* subL)
{
    OBJECT * pvt2 = subF->next();
    if ( comp(subF, pvt2) > 0 )
        move_behind(front, back, subF, subL, subF, pvt2);
    OBJECT * pvt1 = subF;
    while ((pvt2 != subL) and comp(pvt2->next(), pvt2) >= 0 )
        pvt2 = pvt2->next();
    if ( pvt2 != subL )
    {
        OBJECT * down = subL;
        while ( down != pvt2 )
        {
            OBJECT * temp = down->prev();
            if ( comp(pvt1, down) > 0 )
                move_behind(front, back, subF, subL, pvt1, down);
            else if ( comp(pvt2, down) > 0 )
                move_behind(front, back, subF, subL, pvt2, down);
            down = temp;
        }
        if ((pvt1 != subF) and (pvt1->prev() != subF))
            blinksort(front, back, comp, subF, pvt1->prev());
        if ((pvt1->next() != pvt2) and (pvt1->next() != pvt2->prev()))
            blinksort(front, back, comp, pvt1->next(), pvt2->prev());
        if ((pvt2 != subL) and (pvt2->next() != subL))
            blinksort(front, back, comp, pvt2->next(), subL);
    }
}

