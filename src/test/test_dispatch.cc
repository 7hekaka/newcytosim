// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
// [Grand Central Dispatch](https://developer.apple.com/documentation/DISPATCH)

#include <dispatch/dispatch.h>
#include <stdio.h>
#include <cmath>


const size_t CNT = 1024;

double results[CNT];

void work(size_t i)
{
    results[i] = std::sin(i);
}


int main(int argc, char * argv[])
{
    //dispatch_queue_t queue = dispatch_queue_create(nullptr, nullptr)
    dispatch_queue_t const& queue = dispatch_get_global_queue(0, 0);
    
    dispatch_apply(CNT, queue, ^(size_t i){ work(i); });
    
    dispatch_sync(queue, ^{ printf("Hello, from a dispatch queue!\n"); });

    //dispatch_release(queue);
    for ( int i = 0; i < 10; ++i )
        printf("%i %9.3f\n", i, results[i]);
}

