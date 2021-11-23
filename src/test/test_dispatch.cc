// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
// [Grand Central Dispatch](https://developer.apple.com/documentation/dispatch)
// [doc](https://developer.apple.com/documentation/dispatch)


#include <dispatch/dispatch.h>
#include <stdio.h>
#include <cmath>
#include <sys/time.h>


const size_t CNT = 128;

double results[CNT];

void work(void* context, size_t i)
{
    fprintf(stderr, "%lX ", i&15);
    results[i] = i * i;
    usleep((1+(i&3))*200000);
}

void hello(void * context)
{
    usleep(100000);
    const char * str = dispatch_queue_get_label(DISPATCH_CURRENT_QUEUE_LABEL);
    printf("Hello, from `%s'!\n", str);
}


int main(int argc, char * argv[])
{
    time_t sec = time(nullptr);
    
    dispatch_queue_t queue = dispatch_queue_create("QUEUE", DISPATCH_QUEUE_CONCURRENT);
    
    dispatch_async_f(queue, nullptr, hello);
    dispatch_apply_f(CNT, queue, nullptr, work);
    dispatch_release(queue);
    
    for ( int i = 0; i < 10; ++i )
        printf("\n%i %9.3f", i, results[i]);
    
    printf("\n%lu seconds!\n", time(nullptr)-sec);
}

