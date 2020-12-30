// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University
// Created by Francois Nedelec on 24/04/2010.


#include "backtrace.h"

// enable/disable backtrace with the '#if' below:
#if 0

#include <execinfo.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <cxxabi.h>

/**
 * print the current call-stack using functions from the GNU C Library:
 * - backtrace()
 * - backtrace_symbols()
 * .
 * provided by <execinfo.h>
 */
void print_backtrace(int out)
{
    void* buffer[128];
    size_t n_ptr = 512;
    int size = backtrace(buffer, 128);
    if ( size < 2 )
    {
        //(void) write(out, "Empty execution stack!\n", 23);
        return;
    }
#if ( 1 )
    char** buf = backtrace_symbols(buffer, size);
    char * ptr = (char*)malloc(n_ptr);

    (void) write(out, "Cytosim execution stack:\n", 25);
    for ( int i = 1; i < size; ++i )
    {
        int status = -1;
        //printf("%i: %s\n", i, buf[i]);
        char* str = buf[i];
        while ( *str )
        {
            // find start of C++ mangled name
            if ( *str == '_' && *(str+1) == 'Z' )
                break;
            ++str;
        }
        if ( *str )
        {
            char* end = str;
            // find end of string
            while ( *end && *end != ' ' )
                ++end;
            *end = 0;
            ptr = abi::__cxa_demangle(str, ptr, &n_ptr, &status);
        }
        if ( status == 0 )
        {
            (void) write(out, buf[i], str-buf[i]);
            (void) write(out, ptr, strlen(ptr));
        }
        else
            (void) write(out, buf[i], strlen(buf[i]));
        (void) write(out, "\n", 1);
    }
    free(ptr);
    free(buf);
#else
    backtrace_symbols_fd(buffer, size, out);
#endif
}

#else

void print_backtrace(int out)
{
}

#endif

