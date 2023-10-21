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
 *
 * Note that you may need to compile with `-g -rdynamic` to enable this
 */
void print_backtrace(int out)
{
    void* buffer[128];
    int size = backtrace(buffer, 128);
    if ( size < 2 )
    {
        //(void) write(out, "Empty execution stack!\n", 23);
        return;
    }
#if ( 1 )
    size_t len = 512;
    char * ptr = (char*)malloc(len);
    char** buf = backtrace_symbols(buffer, size);

    ssize_t __attribute__((unused)) u;
    u = write(out, "Cytosim execution stack:\n", 25);
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
        char* end = str;
        if ( *str )
        {
            // find end of string
            while ( *end && *end != ' ' )
                ++end;
            *end = 0;
            ptr = abi::__cxa_demangle(str, ptr, &len, &status);
            *end = ' ';
        }
        if ( status == 0 )
        {
            u = write(out, buf[i], str-buf[i]);
            u = write(out, ptr, strlen(ptr));
            u = write(out, end, strlen(end));
        }
        else
            u = write(out, buf[i], strlen(buf[i]));
        u = write(out, "\n", 1);
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

