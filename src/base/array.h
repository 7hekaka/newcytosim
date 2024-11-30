// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef ARRAY_H
#define ARRAY_H

#include "assert_macro.h"
#include "random.h"
#include <iostream>
#include <algorithm>


/** 
 Array<typename VAL> stores objects of class VAL.
 
 This class is similar to std::array<VAL> and resembles std::vector<VAL>,
 but many functions of std::vector are missing, and some functions were
 added as needed:
 - remove_pack() will pack the array removing 'zero' values,
 - sort() will sort the array given a order function,
 - allocate() ensures that array[s-1] can then be accessed, as with C-arrays.
 - operator[](int) returns the object stored at a particular index.
 - shuffle() permutes the values to produce a random ordering.
 - data() returns a pointer to the underlying C-array.
 .
 
 VAL needs to have a public constructor without argument.
 New memory is allocated if necessary by allocate(), and the values
 from the old array are copied to the new memory space.
 
 Allocation when it is done exceeds what is requested, to ensure that allocation
 only occurs from time-to-time, even if objects are added one by one to the array.
 */

/// Dynamically allocated array of VAL
template <typename VAL, unsigned CHUNK = 4>
class Array final
{
public:

    /// typedef for the template argument
    typedef VAL value_type;

    /// iterator class type
    typedef value_type * iterator;
    
    /// const iterator class type
    typedef value_type const* const_iterator;

private:
    
    /// C-array holding the values
    VAL * val_;
    
    /// size of memory that was allocated for val_[]
    size_t alc_;
    
    /// number of objects currently present in the array
    size_t nbo_;
    
#pragma mark -
private:
    
    /// additional size allocated at the end of the array
    static constexpr size_t EXTRA = 0;

    /// the integer above `s` that is a multiple of CHUNK
    size_t chunked(size_t s)
    {
        if ( CHUNK > 0 )
            return ( s + CHUNK - 1 ) & ~( CHUNK - 1 );
        else if ( val_ && alc_ > 0 )
            return 2 * alc_;
        else
            return 4;
    }
    
    /// copy data
    inline void copy(VAL * dst, const VAL* src, const size_t cnt)
    {
        //fprintf(stderr, "Array::copy(%lu) ", cnt);
#if ( 1 )
        // to bypass the constructor/destructor
        memcpy(dst, src, cnt*sizeof(VAL));
#else
        for ( size_t n = 0; n < cnt; ++n )
            dst[n] = src[n];
#endif
    }
    
#pragma mark -
public:
        
    /// Default creator does not allocate
    Array() : val_(nullptr), alc_(0), nbo_(0)
    {
        assert_true(__builtin_popcount(CHUNK) <= 1);
    }

    /// allocate to hold `a` objects and set chunk size to `k`
    Array(size_t a) : alc_(0), nbo_(0)
    {
        assert_true(__builtin_popcount(CHUNK) <= 1);
        //printf("Array %p constructor size %i\n", this, alc_);
        alc_ = chunked(a);
        if ( alc_ > 0 )
            val_ = new VAL[alc_+EXTRA];
        else
            val_ = nullptr;
    }
    
    /// create list with only one entry
    Array(VAL arg) : alc_(1), nbo_(1)
    {
        //printf("Array %p constructor(obj)\n", this);
        val_ = new VAL[1+EXTRA];
        val_[0] = arg;
    }

    /// Copy constructor
    Array(const Array<VAL, CHUNK>& o) : val_(nullptr), alc_(0), nbo_(o.nbo_)
    {
        if ( o.alc_ )
        {
            //printf("Array %p copy constructor size %i\n", this, nbo_);
            allocate(o.alc_);
            copy(val_, o.val_, o.alc_);
        }
    }

    /// Destructor
    ~Array()
    {
        deallocate();
    }
    
    /// Copy assignment operator
    Array& operator = (Array<VAL, CHUNK> const& o)
    {
        if ( o.nbo_ > alc_ )
        {
            //printf("Array %p allocated %i = from size %i\n", this, alc_, o.nbo_);
            deallocate();
            allocate(o.nbo_);
        }
        nbo_ = o.nbo_;
        copy(val_, o.val_, nbo_);
        return *this;
    }
    
    /// Move assignment operator
    Array& operator = (Array<VAL, CHUNK> && o)
    {
        //printf("Array %p <---move--- %p\n", this, o);
        delete[] val_;
        nbo_ = o.nbo_;
        alc_ = o.alc_;
        val_ = o.val_;
        o.val_ = nullptr;
        o.alc_ = 0;
        o.nbo_ = 0;
        return *this;
    }

#pragma mark -
    
    /// Number of objects
    size_t size() const
    {
        return nbo_;
    }

    /// true if this Array holds no value
    bool empty() const
    {
        return ( nbo_ == 0 );
    }
    
    /// Currently allocated size
    size_t capacity() const
    {
        return alc_;
    }
    
    /// Address of the underlying C-array
    VAL * data()
    {
        return val_;
    }
    
    /// Address of the underlying C-array
    VAL const * data() const
    {
        return val_;
    }
    
    /// Address of the underlying C-array
    VAL const * addr(const size_t i) const
    {
        return val_+i;
    }

    /// pointer to first element
    iterator begin() const
    {
        return val_;
    }
    
    /// pointer to a position just past the last element
    iterator end() const
    {
        return val_+nbo_;
    }
    
    /// reference to Object at index i (val_[i])
    VAL& at(const size_t i) const
    {
        assert_true( i < alc_ );
        return val_[i];
    }
    
    /// reference to Object that is beyond last one
    VAL& past_back() const
    {
        return val_[nbo_];
    }
    
    /// increase size by `arg`
    void increase_size(const size_t arg)
    {
        assert_true( nbo_+arg < alc_ );
        nbo_ += arg;
    }
 
    /// reference to Object at index i (val_[i])
    VAL& operator[](const size_t i) const
    {
        assert_true( i < nbo_ );
        return val_[i];
    }
    
    /// return element at index 0
    VAL const& front() const
    {
        assert_true( 0 < nbo_ );
        return val_[0];
    }
    
    /// return element at index 0
    VAL& front()
    {
        assert_true( 0 < nbo_ );
        return val_[0];
    }

    /// return last element
    VAL const& back() const
    {
        assert_true( 0 < nbo_ );
        return val_[nbo_-1];
    }
    
    /// return last element
    VAL& back()
    {
        assert_true( 0 < nbo_ );
        return val_[nbo_-1];
    }

#pragma mark -
    
    /// Allocate to hold `s` objects: valid indices are 0 <= indx < max
    void reallocate(const size_t alc_new)
    {
        VAL * val_new = new VAL[alc_new+EXTRA];
        if ( val_ )
        {
            // copy over valid data
            copy(val_new, val_, std::min(nbo_, alc_new));
            delete[] val_;
        }
        //fprintf(stderr, "  reallocate(%lu -> %lu)\n", alc_, alc_new);
        alc_ = alc_new;
        val_ = val_new;
    }
    
    /// Allocate to hold at least `s` objects: valid indices are 0 <= indx < max
    size_t allocate(const size_t s)
    {
        if ( s > alc_ )
        {
            reallocate(chunked(s));
            assert_true( alc_ >= s );
            return s;
        }
        return 0;
    }
    
    /// Allocate to hold at least `s` objects: valid indices are 0 <= indx < max
    void reserve(const size_t s)
    {
        allocate(s);
    }
    
    /// Allocate and set new values to `zero`
    size_t allocate_zero(const size_t arg, VAL const& zero)
    {
        size_t res = allocate(arg);
        if ( res )
        {
            //set the newly allocated memory to zero
            for ( size_t i = nbo_; i < alc_; ++i )
                val_[i] = zero;
        }
        return res;
    }
    
    /// Reduce size to min(arg, actual_size), keeping elements starting from index 0
    void truncate(const size_t arg)
    {
        nbo_ = std::min(arg, nbo_);
    }
    
    /// Set size to `arg`, allocating if necessary
    void resize(const size_t arg)
    {
        if ( arg < nbo_ )
            nbo_ = arg;
        else if ( arg > nbo_ )
        {
            allocate(arg);
            nbo_ = arg;
        }
    }
    
    /// Release allocated memory
    void deallocate()
    {
        //printf("Array %p deallocate %i\n", this, allocated);
        if ( alc_ ) delete[] val_;
        val_ = nullptr;
        alc_ = 0;
        nbo_ = 0;
    }
    
    /// Set the number of objects to zero
    inline void clear()
    {
        nbo_ = 0;
    }
    
    /// Set the number of objects to zero
    inline void clear(VAL const& zero)
    {
        for ( size_t i=0; i < nbo_; ++i )
            val_[i] = zero;
        nbo_ = 0;
    }
    
    /// Set all values to `zero`
    void zero(VAL const& zero)
    {
        assert_true( val_ || alc_==0 );
        for ( size_t i=0; i < alc_; ++i )
            val_[i] = zero;
    }

    /// Delete all values as if they were pointers to Object
    void destroy()
    {
        assert_true( val_ || nbo_==0 );
        for ( size_t i=0; i < nbo_; ++i )
        {
            //std::clog << " delete " << val_[i] << '\n';
            delete(val_[i]);
            val_[i] = nullptr;
        }
        nbo_ = 0;
    }
    
#pragma mark -
    
    /// Increment the size of the array, and return new value at end of it
    VAL& new_val()
    {
        if ( nbo_ >= alc_ )
            reallocate(chunked(nbo_+1));
        VAL& res = val_[nbo_++];
        return res;
    }
    
    /// Add `np` at the end of this Array
    void push_back(VAL const& v)
    {
        if ( nbo_ >= alc_ )
            reallocate(chunked(nbo_+1));
        val_[nbo_++] = v;
    }
    
    template < typename... Args >
    void emplace(Args&&... args)
    {
        if ( nbo_ >= alc_ )
            reallocate(chunked(nbo_+1));
        val_[nbo_++] = VAL(args...);
    }

    /// remove last element
    void pop_back()
    {
        --nbo_;
    }

    /// Add the elements of `array` at the end of this Array
    void append(const Array<VAL, CHUNK> array)
    {
        allocate(nbo_+array.nbo_);
        for ( size_t i = 0; i < array.nbo_; ++i )
            val_[i+nbo_] = array.val_[i];
        nbo_ += array.nbo_;
    }
    
    /// Add the elements of `array` at the end of this Array
    void append_except(const Array<VAL, CHUNK> array, VAL const& v)
    {
        allocate(nbo_+array.nbo_);
        for ( size_t i = 0; i < array.nbo_; ++i )
            if ( array.val_[i] != v )
                val_[nbo_++] = array.val_[i];
    }

    /// Return index of `obj`, or ~0 if not found in the list (linear search)
    size_t index(const VAL obj) const
    {
        assert_true( val_ || nbo_==0 );
        for ( size_t i = 0; i < nbo_; ++i )
            if ( val_[i] == obj )
                return i;
        return ~0UL;
    }
    
    /// Replace `old_value` by `new_value`, or return false if `old_value` is not found
    bool replace(VAL const& old_value, VAL const& new_value)
    {
        assert_true( val_ || nbo_==0 );
        for ( size_t i=0; i < nbo_; ++i )
        {
            if ( val_[i] == old_value )
            {
                val_[i] = new_value;
                return true;
            }
        }
        return false;
    }
    
#pragma mark -
    
    /// set `np` at any position equal to `zero`, or at the end of the array
    void push_pack(const VAL np, VAL const& zero)
    {
        assert_true( val_ || nbo_==0 );
        for ( size_t i = 0; i < nbo_; ++i )
        {
            if ( val_[i] == zero )
            {
                val_[i] = np;
                return;
            }
        }
        push_back(np);
    }

    
    /// Returns the number of occurence of 'val' in the array
    size_t count(VAL const& v) const
    {
        if ( !val_ || nbo_==0 )
            return 0;
        size_t res = 0;
        for ( size_t i = 0; i < nbo_; ++i )
            if ( val_[i] == v ) ++res;
        return res;
    }

    /// Number of values which are different from `val` in the array
    size_t count_except(VAL const& v) const
    {
        if ( val_ == 0 || nbo_==0 )
            return 0;
        size_t res = 0;
        for ( size_t i = 0; i < nbo_; ++i )
            if ( val_[i] != v ) ++res;
        return res;
    }

    
    /**
     Remove all entries which are equal to `zero`, and pack array by shuffling values around.
     The order of the elements is not preserved, and copy operations are minimized
     */
    template <typename T>
    static T * remove_pack(T * s, T * e, T const& zero)
    {
        if ( e <= s )
            return s;
        --e;
        while ( s < e )
        {
            // find the next `zero` going upward:
            while ( *s != zero )
            {
                ++s;
                if ( e <= s )
                    return e + ( *e != zero );
            }
            // going downward, skip `zero` values:
            while ( *e == zero )
            {
                --e;
                if ( e <= s )
                    return e;
            }
            // flip the two values:
            *s = *e;
            *e = zero; // maybe not necessary
            ++s;
            --e;
        }
        return e + ( *e != zero );
    }
    
    
    

    /// Remove all entries which are equal to `zero`, and pack array
    void remove_pack(VAL const& zero)
    {
        assert_true( val_ || nbo_==0 );
        nbo_ = remove_pack(val_, val_+nbo_, zero) - val_;
    }
    
    
    /// Sort array using `std::qsort()` and the provided comparison function
    void quick_sort(int (*comp)(const void *, const void *))
    {
        assert_true( val_ || nbo_==0 );
        qsort(val_, nbo_, sizeof(VAL), comp);
    }
    
    void sort()
    {
        //std::sort using a lambada function
        std::sort(val_, val_+nbo_, [](VAL const& a, VAL const& b) { return a < b; });
    }
    
    /// Return one of the value in the array, chosen randomly
    VAL& pick_one()
    {
        assert_true( nbo_ > 0 );
        assert_true( nbo_ <= UINT32_MAX );
        return val_[RNG.pint32(static_cast<uint32_t>(nbo_))];
    }

    /// Move the last Object on top, pushing all other values down by one slot
    void rotate()
    {
        if ( nbo_ > 1 )
        {
            assert_true(val_);
        
            VAL * tmp = val_[0];
            for ( size_t i = 0; i < nbo_-1; ++i )
                val_[i] = val_[i+1];
            val_[nbo_-1] = tmp;
        }
    }
    
    /// exchange the values of `a` and `b`
    static void swap(VAL* a, VAL* b)
    {
        VAL tmp = *a;
        *a = *b;
        *b = tmp;
    }
    
    /// Swap two random values in the array
    void permute()
    {
        assert_true(val_);
        assert_true( nbo_ <= UINT32_MAX );
        size_t ii = RNG.pint32(static_cast<uint32_t>(nbo_));
        size_t jj = RNG.pint32(static_cast<uint32_t>(nbo_));
        if ( ii != jj )
            swap(val_+ii, val_+jj);
    }
    
    
    /// Randomly permutes all objects in the array
    /**
     Fisher-Yates shuffle
     This produces uniform shuffling in linear time.
     see Knuth's The Art of Programming, Vol 2 chp. 3.4.2 
     */
    void shuffle32()
    {
        assert_true( nbo_ <= UINT32_MAX );
        assert_true( val_ || nbo_==0 );
        uint32_t jj = (uint32_t)nbo_, kk;
        while ( jj > 1 )
        {
            kk = RNG.pint32(jj);  // 32 bits in [0, j-1]
            --jj;
            swap(val_+jj, val_+kk);
        }
    }
    
    /// Randomly permutes all objects in the array
    /**
     Fisher-Yates shuffle
     This produces uniform shuffling in linear time.
     see Knuth's The Art of Programming, Vol 2 chp. 3.4.2
     */
    void shuffle64()
    {
        assert_true( val_ || nbo_==0 );
        uint64_t jj = nbo_, kk;
        while ( jj > 1 )
        {
            kk = RNG.pint64(jj);  // 64 bits in [0, j-1]
            --jj;
            swap(val_+jj, val_+kk);
        }
    }

    /// Randomly permutes all objects in the array
    void shuffle()
    {
        shuffle32();
    }
    
    /// possibly reduce array to a maximum of 'cnt' objects, chosen randomly
    void shuffle_truncate(size_t cnt)
    {
        if ( cnt < size() )
        {
            if ( cnt == 1 )
                val_[0] = pick_one();
            else
                shuffle();
            truncate(cnt);
        }
    }
};


#endif
