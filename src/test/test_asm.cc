

#include <cstddef>
#include <cstdio>

const size_t SIZE = 8;

/**
 This is code by Katsura, from
 https://katsurashareware.com/forum/viewtopic.php?f=9&t=16&sid=29d42f2f4ceb346c30efb844d6656a4a
 
 https://developer.arm.com/documentation
 Non-temporal loads
 In AArch64 state, cache requests made by a non-temporal load instruction (LDNP) are treated as normal loads.
 */
void *ntbread(void *dst, const void *src, size_t bs)
{
    __asm__ __volatile__ (
        "BREAD1_%=:\n"
        "ldnp   q0, q1, [%0]\n"
        "ldnp   q2, q3, [%0, #0x20]\n"
        "add    %0, %0, #0x40\n"
        "subs   %1, %1, #0x40\n"
        "b.hi   BREAD1_%=\n"
        : "=r"(src) /* %0 */, "=r"(bs) /* %1 */ /* output */
        : "0"(src) /* %0 */, "1"(bs) /* %1 */ /* input */
        : "v0", "v1", "v2", "v3" /* clobbered */
    );
    return dst;
}


void test()
{
    int src[SIZE] = { 0 };
    int dst[SIZE] = { 0 };
    
    for ( int i = 0; i < SIZE; ++i )
        src[i] = i;
    
    ntbread(dst, src, SIZE*sizeof(int));
    
    for ( int i = 0; i < SIZE; ++i )
        printf("%i ", src[i]);
    printf("\n");
    
    for ( int i = 0; i < SIZE; ++i )
        printf("%i ", dst[i]);
    printf("\n");
}


int main(int argc, char* argv[])
{
    test();
}
