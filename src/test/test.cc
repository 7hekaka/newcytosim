// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <cstdio>
#include <sys/types.h>

#include "array.h"
#include "random.h"


int compare(const void * a, const void * b)
{
    int av = *static_cast<const int*>(a);
    int bv = *static_cast<const int*>(b);
    return ( av > bv ) - ( av < bv );
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
            a.push_back(RNG.pint(2));
        
        printf("\nsize %lu", a.size());
        {
            Array<int> b;
            b = a;
            a.deallocate();
            
            printf("\n   copy %2lu :", b.size());
            for( size_t i=0; i < b.size(); ++i )
                printf(" %i", b[i]);

            b.remove_pack(0);
            
            printf("\n   pack %2lu :", b.size());
            for( size_t i=0; i < b.size(); ++i )
                printf(" %i", b[i]);
            
            b.sort(compare);
            
            printf("\n   sort %2lu :", b.size());
            for( size_t i=0; i < b.size(); ++i )
                printf(" %i", b[i]);
        }
    }
    
    printf("\ndone\n");
    return 0;
}
