// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#include <cstdio>

#include "array.h"
#include "random.h"


int compare(const void * A, const void * B)
{
    int a = *static_cast<const int*>(A);
    int b = *static_cast<const int*>(B);
    return ( a > b ) - ( a < b );
}


int main(int argc, char* argv[])
{
    RNG.seed();
    Array<int> a;
    
    for( size_t cnt = 0; cnt < 10; ++cnt )
    {
        a.clear();
        size_t n = RNG.poisson(8);
        for( size_t i=0; i < n; ++i )
            a.push_back(RNG.pint32(3));
        
        printf("\nsize %lu", a.size());
        {
            Array<int> b;
            b = a;
            a.deallocate();
            
            printf("\n   copy %2lu :", b.size());
            for( int i : b )
                printf(" %i", i);

            b.remove_pack(0);
            
            printf("\n   pack %2lu :", b.size());
            for( int i : b )
                printf(" %i", i);
            
            b.quick_sort(compare);
            
            printf("\n   sort %2lu :", b.size());
            for( int i : b )
                printf(" %i", i);
        }
    }
    
    printf("\ndone\n");
}
