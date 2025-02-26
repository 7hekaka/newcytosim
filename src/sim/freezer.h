// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University.
#ifndef FREEZER_H
#define FREEZER_H


/// holds a list of Single/Couple with identical Property.
/**
 The `freezer` list is built using `Object::next()` and thus can only hold
 objects that are not linked in SingleSet/CoupleSet already (see below).
 The single-linked list permits addition/removal at one end only: the head.
 These lists are used to buffer creation/deletion of Single/Couple.
 */
template < typename OBJECT >
class Freezer
{
    /// Pointer to first member in freezer
    OBJECT * head_;
    
    /// Number of elements in freezer
    size_t count_;

public:
    
    /// constructor
    Freezer() : head_(nullptr), count_(0) { }
    
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
    
    /// reset count and head
    void clear()
    {
        head_ = nullptr;
        count_ = 0;
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
    
    /// count number of objects in list
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
