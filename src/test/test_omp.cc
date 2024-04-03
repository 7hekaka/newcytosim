// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University
// F. Nedelec, 22.04.2018

/*
 To use OpenMP on MacOS:
 
 1. Use GCC provided by Homebrew.
 
 2. Use clang, with libomp (brew install libomp)

      clang -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp src/test/test_omp.cc
 
 */


#include <cstdio>
#include "omp.h"

typedef double real;

const int max = 1024;
const int num = 4;

void process(int val)
{
    printf("%i %i\n", omp_get_thread_num(), val);
}

/*
void pointers(real * vec, real * end)
{
#pragma omp parallel
    {
        #pragma omp single private(p)
        {
            real * p = vec;
            while ( p < end )
            {
                #pragma omp task
                process(p);
                p += 10;
            }
        }
    }
}
*/

void pooh(int off, real * vec, int inc)
{
    printf("thread %i\n", off);
    for ( int i = off; i < max; i += inc)
        vec[i] = 0;
}


int main(int argc, char* argv[])
{
    double vec[max];
    
    omp_set_num_threads(num);
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        pooh(id, vec, num);
    }
 
    #pragma omp parallel num_threads(6)
    {
        int inc = omp_get_num_threads();
        int id = omp_get_thread_num();
        pooh(id, vec, inc);
    }
    
    #pragma omp parallel for num_threads(4)
    for ( int i = 0; i < 11; ++i )
        process(i);
 
    printf("done\n");
}
