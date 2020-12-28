// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include <cstdio>
#include <sys/types.h>

#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <cstdint>
#include <cmath>

#define PRINT(arg) printf("sizeof %16s   %lu bytes\n", #arg, sizeof(arg));


void print_types()
{
    PRINT(bool);
    PRINT(char);
    PRINT(short);
    PRINT(int);
    PRINT(unsigned);
    PRINT(long int);
    PRINT(long long int);
    
    PRINT(int_fast8_t);
    PRINT(int_fast16_t);
    PRINT(int_fast32_t);
    PRINT(int_fast64_t);

    PRINT(float);
    PRINT(double);
    PRINT(long double);
    
    PRINT(void *);
    PRINT(fpos_t);
    PRINT(off_t);
    PRINT(size_t);
}

void print_sizes()
{
    printf("Size of Data Types:\n");

    printf("bool          = %lu bytes\n", sizeof(bool) );
    printf("char          = %lu bytes\n", sizeof(char) );
    printf("short         = %lu bytes\n", sizeof(short) );
    printf("int           = %lu bytes\n", sizeof(int) );
    printf("unsigned      = %lu bytes\n", sizeof(unsigned) );
    printf("long int      = %lu bytes\n", sizeof(long int) );
    printf("long long int = %lu bytes\n", sizeof(long long int) );
    printf("float         = %lu bytes\n", sizeof(float) );
    printf("double        = %lu bytes\n", sizeof(double) );
    printf("long double   = %lu bytes\n", sizeof(long double) );
    printf("void *        = %lu bytes\n", sizeof(void *) );
    printf("fpos_t        = %lu bytes\n", sizeof(fpos_t) );
    printf("off_t         = %lu bytes\n", sizeof(off_t) );
    printf("size_t        = %lu bytes\n", sizeof(size_t) );
}

int main()
{
    print_types();
    //printf("\n");
    //print_sizes();
    printf("done\n");
}
