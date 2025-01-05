// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University.
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
*/template<typename OBJECT>
class BubbleSortJob
{
private:

    void push_after(OBJECT *& back, OBJECT * p, OBJECT * n)
    {
        n->prev_ = p;
        n->next_ = p->next_;
        if ( p->next_ )
            p->next_->prev_ = n;
        else
            back = n;
        p->next_ = n;
    }
    
    void push_before(OBJECT *& front, OBJECT * p, OBJECT * n)
    {
        n->next_ = p;
        n->prev_ = p->prev_;
        if ( p->prev_ )
            p->prev_->next_ = n;
        else
            front = n;
        p->prev_ = n;
    }
    
public:
    
    void sort(OBJECT *& front, OBJECT *& back, int (*COMP)(const OBJECT *, const OBJECT *))
    {
        OBJECT * ii = front;
        
        if ( !ii )
            return;
        
        ii = ii->next_;
        
        while ( ii )
        {
            OBJECT * kk = ii->next_;
            OBJECT * jj = ii->prev_;
            
            if ( COMP(ii, jj) > 0 )
            {
                jj = jj->prev_;
                
                while ( jj && COMP(ii, jj) > 0 )
                    jj = jj->prev_;
                
                //pop(ii):
                ii->prev_->next_ = ii->next_;
                if ( ii->next_ )
                    ii->next_->prev_ = ii->prev_;
                else
                    back = ii->prev_;
                
                if ( jj )
                    push_after(back, jj, ii);
                else {
                    ii->next_ = front;
                    ii->prev_ = nullptr;
                    front->prev_ = ii;
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
        while ( head != tail && head->next_ != tail )
        {
            head = head->next_;
            tail = tail->prev_;
        }
        
        if ( head == tail )
            tail = head->next_;
        
        head->next_ = nullptr;
        tail->prev_ = nullptr;
    }
    
    
    void mergelists_(OBJECT *& head, OBJECT * node, OBJECT * from, OBJECT *& tail)
    {
        assert_true( head && from );
        if ( COMP(head, from) < 0 )
        {
            OBJECT * temp = head->next_;
            if ( temp )
            {
                // head remains on top of the list
                mergelists_(temp, node, from, tail);
                head->next_ = temp;
                temp->prev_ = head;
            }
            else
                head->next_ = from;
        }
        else
        {
            if ( from->next_ )
                mergelists_(head, node, from->next_, tail);
            // add `from` on top of the list
            head->prev_ = from;
            from->next_ = head;
            from->prev_ = nullptr;
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
 
 This uses the COMP() function provided,
 COMP(A, B) should return:
     - a positive integer if 'A > B'
     - a negative integer if 'B < A'
 for the list to be sorted in order of increasing values
 */
template<typename OBJECT>
class BlinkSortJob
{
private:
    
    OBJECT * front_;
    OBJECT * back_;
    
    int (*COMP)(const OBJECT *, const OBJECT *);
    
    void movebehind_(OBJECT *& subF, OBJECT *& subL, OBJECT * pvt, OBJECT * down)
    {
        if ( down != subL )
            down->next_->prev_ = down->prev_;
        else {
            if ( down == back_ )
                back_ = down->prev_;
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
            if ( pvt == front_ )
                front_ = down;
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
                if ( COMP(pvt1, down) > 0 )
                    movebehind_(subF, subL, pvt1, down);
                else if ( COMP(pvt2, down) > 0 )
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
    
    BlinkSortJob() : front_(nullptr), back_(nullptr), COMP(nullptr) {}

    void sort(OBJECT *& front, OBJECT *& back, int (*comp)(const OBJECT *, const OBJECT *))
    {
        if ( front != back )
        {
            COMP = comp;
            front_ = front;
            back_ = back;
            blinksort_(front_, back_);
            front = front_;
            back = back_;
        }
    }
};
