// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include <cstdlib>
#include <cstdio>
#include <pthread.h>

/*
 This is an elementary test where the execution of child threads is controlled
 from the their parent thread, using a condition signal.
 Using the C POSIX thread library
 */

pthread_mutex_t mutex;
pthread_cond_t condition;

thread_local void* thread_id = 0;


void* loop(void *arg)
{
    thread_id = arg;
    pthread_mutex_lock(&mutex);
    pthread_t self = pthread_self();
    for ( int cnt = 1; cnt <= 10; ++cnt )
    {
        fprintf(stderr, "-- thread %p @ %i -- ", thread_id, cnt);
        pthread_cond_wait(&condition, &mutex);
    }
    pthread_mutex_unlock(&mutex);
    return (void*)self;
}


int main(int argc, char *argv[]) 
{
    fprintf(stderr, "type `return` to signa and `q` to quit\n");
    
    pthread_mutex_init(&mutex, 0);
    pthread_cond_init(&condition, 0);
    
    pthread_t slave, child;
    pthread_create(&slave, 0, &loop, (void*)1);
    pthread_create(&child, 0, &loop, (void*)2);

    char key[32];
    while ( fgets(key, sizeof(key), stdin) )
    {
        if ( key[0] == 'q' ) break;
        if ( key[0] == '\n' )
            pthread_cond_signal(&condition);
    }
    
    if ( 1 )
    {
        pthread_cancel(slave);
        pthread_join(slave, 0);
    }
    
    pthread_cond_destroy(&condition);
    pthread_mutex_destroy(&mutex);
    fprintf(stderr, "goodbye!\n");
}
