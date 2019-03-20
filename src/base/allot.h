// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef ALLOT_H
#define ALLOT_H

#include <cstdio>
#include <cstdlib>


/// Allot manages a piece of memory (this class is deprecated)
/**
 Allot holds an array of type T and remembers the size of the array.
 The destructor will release this memory by calling delete[]
 */

template <typename VAL>
class Allot
{
public:
    
    /// copy of the first template argument
    typedef VAL value_type;

private:
    
    /// size currently allocated
    size_t  alc_;
    
    /// array allocated
    VAL *   val_;
    
    /// copy old memory to new one
    bool    cpy_;

    /// chunk size (a power of 2)
    size_t  chk_;
    
private:
    
    /// the integer above s that is a multiple of mChunk
    size_t chunked(size_t s)
    {
        return ( s + chk_ - 1 ) & ~( chk_ - 1 );
    }
    
    /// return smallest power of 2 that is greater or equal to `s`
    size_t next_power(size_t s)
    {
        if ( s & (s-1) )
        {
            do
                s &= s-1;
            while ( s & (s-1) );
            
            s <<= 1;
        }
        return s;
    }
    
    /// copy data
    inline void copy(VAL* to, VAL const* from, const size_t cnt)
    {
        for ( size_t ii=0; ii<cnt; ++ii )
            to[ii] = from[ii];
    }

public:
    
    /// constructor
    Allot() : alc_(0), val_(0), cpy_(0), chk_(8) { }
    
    
    /// Allocate size `s`, and set copy flag to `cop`
    Allot(bool c) : alc_(0), val_(0), cpy_(c), chk_(8)
    {
    }
    
    /// Allocate size `s`, set copy flag to `cop`, and chunk size to `chk`
    Allot(bool c, size_t k) : alc_(0), val_(0), cpy_(c)
    {
        if ( !k )
        {
            fprintf(stderr, "Allot::chunk must not be null");
            exit(1);
        }
        
        chk_ = next_power(k);
    }
    
    /// copy constructor
    Allot(Allot<VAL> const& o) : alc_(0), val_(0), cpy_(o.cpy_), chk_(o.chk_)
    {
        if ( o.alc_ )
        {
            allocate(o.alc_);
            copy(val_, o.val_, o.alc_);
        }
    }
    
    /// release memory
    void deallocate()
    {
        if ( val_ )
        {
            delete[] val_;
            val_ = 0;
        }
        alc_ = 0;
    }
    
    /// destructor
    ~Allot()
    {
        deallocate();
    }
    
    /// copy assignment operator
    Allot& operator = (Allot<VAL> const& o)
    {
        deallocate();
        if ( o.alc_ )
        {
            allocate(o.alc_);
            copy(val_, o.val_, o.alc_);
        }
        return *this;
    }
    
    /// change size of allocated memory
    void reallocate(size_t alc_new)
    {
        //std::clog << "Allot " << this << ":reallocate(" << alc_new << ")"<<std::endl;
        VAL * val_new = new VAL[alc_new];
        if ( val_ )
        {
            if ( cpy_ )
            {
                if ( alc_ < alc_new )
                    copy(val_new, val_, alc_);
                else
                    copy(val_new, val_, alc_new);
            }
            delete[] val_;
        }
        alc_ = alc_new;
        val_ = val_new;
    }
    
    /// allocate, but only if size increases
    size_t allocate(size_t s)
    {
        if ( s > alc_ )
        {
            reallocate(chunked(s));
            return s;
        }
        return 0;
    }
    
    /// forget current allocation
    VAL * release()
    {
        VAL * res = val_;
        val_ = 0;
        alc_ = 0;
        return res;
    }
    
    /// exchange the data between `this` and `o`
    void swap(Allot<VAL> o)
    {
        VAL * v = val_;
        size_t a = alc_;
        val_ = o.val_;
        alc_ = o.alc_;
        o.val_ = v;
        o.alc_ = a;
    }
    
    /// allocated size
    size_t capacity() const
    {
        return alc_;
    }

    /// pointer to data array
    VAL const* data() const
    {
        return val_;
    }

    /// pointer to data array
    VAL* data()
    {
        return val_;
    }
    
    /// conversion to array
    operator VAL * () { return val_; }
    
    /// conversion to const array
    operator VAL const* () const { return val_; }

    /// access to array
    VAL& operator [] (size_t i) { return val_[i]; }

};

#endif

