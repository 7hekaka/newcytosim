// Cytosim was created by Francois Nedelec. Copyright 2019 Cambridge University

// replace global new & delete operators to control memory alignment

#include <new>
#include <cstdio>
#include <cstdlib>

/**
 Cytosim's new() operator should replace the default operator.
 It returns memory aligned to a 32-bytes boundary.
 */
void* operator new(std::size_t size)
{
    //printf("new(%lu)\n", size);
    void * ptr = nullptr;
#if ( 0 )
    constexpr std::size_t sup = 1 << 30;
    if ( size > sup )
    {
        std::printf("Error: excessive memory requested (%5zu)\n", size);
        throw std::bad_alloc();
    }
#endif
#if ( 0 )
    // get memory aligned to 32 bytes
    if ( posix_memalign(&ptr, 32, size) )
        throw std::bad_alloc();
#else
    // system's default allocation, not necessarily aligned
    ptr = std::malloc(size);
    //if ( size > 2048 ) std::printf("Cytosim:new(%5zu) %p %+li\n", size, ptr, ((uintptr_t)ptr & 63));
#endif
    if ( !ptr )
        throw std::bad_alloc();
#if ( 0 )
    static char * old = nullptr;
    std::printf("Cytosim:new(%5zu) %p %+li\n", size, ptr, (char*)ptr-old);
    old = (char*)ptr;
#endif
    return ptr;
}


void operator delete(void * ptr) throw()
{
    //std::printf("Cytosim:delete(%p)\n", ptr);
    std::free(ptr);
}


/*
void* operator new[](std::size_t s) throw(std::bad_alloc)
{
    std::printf("Cytosim new[] %5zu\n", s);
    return ::operator new(s);
}

void operator delete[](void *ptr) throw()
{
    std::printf("Cytosim delete[]    %p\n", ptr);
    ::operator delete(ptr);
}
*/
