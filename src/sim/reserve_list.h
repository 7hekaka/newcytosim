// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University.
#ifndef RESERVE_LIST_H
#define RESERVE_LIST_H


/// holds a list of Single/Couple with identical Property.
/**
 This list is build using 'Object::next()' and thus can only hold
 objects that are not linked in SingleSet/CoupleSet already (see below).
 It is a single-linked list and addition/removal is made at the head only.
 The Reserve list are used to buffer creation/deletion of Single/Couple.
 */
template < typename OBJECT >
class ReserveList
{
    /// Pointer to first member in list
    OBJECT * head_;
    /// Number of elements in list
    size_t count_;

public:
    
    /// constructor
    ReserveList() : head_(nullptr), count_(0) { }
    
    /// number of objects stored
    size_t size() const { return count_; }
    
    /// first object
    OBJECT * head() const { return head_; }
    
    /// add object
    void push(OBJECT * obj)
    {
        obj->objset(nullptr);
        obj->next(head_);
        head_ = obj;
        ++count_;
    }
    
    /// remove first object in list
    void pop()
    {
        --count_;
        head_ = head_->next();
        assert_true( count_ > 0 || head_ == nullptr );
    }
    
    /// delete all objects
    void erase()
    {
        OBJECT * obj = head_;
        while ( obj )
        {
            pop();
            delete(obj);
            obj = head();
        }
    }
    
    size_t recount() const
    {
        size_t res = 0;
        OBJECT * obj = head_;
        while ( obj )
        {
            ++res;
            obj = obj->next();
        }
        return res;
    }
};

#endif
