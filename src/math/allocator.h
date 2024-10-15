// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef ALLOCATOR_H
#define ALLOCATOR_H

#include "real.h"
#include "assert_macro.h"
#include <cstdio>

/// Iterative methods to solve a system of linear equations
namespace LinearSolvers
{
    /// allocates vectors of real
    /** An class to allocate and keep track of memory used for vectors */
    class Allocator
    {
    private:
        
        /// length of vectors to be allocated
        size_t siz_;
        
        /// number of vectors allocated
        size_t alc_;
        
        /// memory
        real * mem_;
        
    public:
        
        /// initialize
        Allocator() { siz_ = 0; alc_ = 0; mem_ = nullptr; }

        /// calls release()
        ~Allocator() { deallocate(); }
        
        /// allocate n vectors of size `s`
        void allocate(size_t s, size_t n)
        {
            // Keep memory aligned
            // pad with 4 doubles to allow SIMD instructions overspill:
            siz_ = chunk_real(s) + 4;
            size_t a = siz_ * n;
            if ( a > alc_ )
            {
                free_real(mem_);
                mem_ = new_real(a);
                alc_ = a;
                //fprintf(stdout, "Allocator::allocate %lu bytes\n", a*sizeof(real));
            }
        }
        
        /// release memory
        void deallocate()
        {
            free_real(mem_);
            //fprintf(stdout, "Allocator::releases %lu bytes\n", alc*sizeof(real));
            mem_ = nullptr;
            alc_ = 0;
            siz_ = 0;
        }
        
        /// called to declare that memory is not needed anymore
        void release()
        {
            //deallocate();
        }
        
        /// return pointer to the memory allocated for i-th vector
        real * bind(size_t i) const
        {
            if ( !mem_ )
            {
                fprintf(stderr, "Error: Allocator::allocate() not called\n");
                return nullptr;
            }
            if ( (i+1)*siz_ > alc_ )
            {
                fprintf(stderr, "Error: inconsistent declaration in Allocator\n");
                return nullptr;
            }
            real * ptr = mem_ + i * siz_;
            //fprintf(stdout, "Allocator::bind(%u) at %p\n", i, ptr);
            zero_real(siz_, ptr);
            return ptr;
        }
    };

}

#endif

