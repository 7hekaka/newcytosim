// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include <cstdlib>
#include <csignal>
#include <unistd.h>
#include <exception>
#include <new>

#include "signal_handlers.h"
#include "backtrace.h"

//---------------------------- Signal Handlers --------------------------------

static void out_of_memory_handler()
{
    ssize_t __attribute__((unused)) u;
    u = write(STDERR_FILENO, "\n* * * * *\n", 11);
    u = write(STDERR_FILENO, "Cytosim: memory allocation failed", 33);
    u = write(STDERR_FILENO, "\n* * * * *\n", 11);
    print_backtrace();
    _exit(1);
}

static void termination_handler()
{
    ssize_t __attribute__((unused)) u;
    u = write(STDERR_FILENO, "\n* * * * *\n", 11);
    u = write(STDERR_FILENO, "Cytosim: uncaught exception", 27);
    u = write(STDERR_FILENO, "\n* * * * *\n", 11);
    print_backtrace();
    abort();
}

static void signal_handler(int sig)
{
    ssize_t __attribute__((unused)) u;
    u = write(STDERR_FILENO, "\n*  *  *  *  *\n", 15);
    psignal(sig, "Cytosim");
    u = write(STDERR_FILENO, "*  *  *  *  *\n", 14);
    print_backtrace();
    _exit(sig);
}


void register_signal_handlers()
{
    // Register a function to be called if operator new fails:
    std::set_new_handler(out_of_memory_handler);
    
    // Register a function to be called upon abortion:
    std::set_terminate(termination_handler);
    
    // Register a function to be called for Floating point exceptions:
    if ( signal(SIGFPE, signal_handler) == SIG_ERR )
        write(STDERR_FILENO, "no  SIGFPE handler\n", 19);
#if 1
    if ( signal(SIGSEGV, signal_handler) )
        write(STDERR_FILENO, "no SIGSEGV handler\n", 19);
              
    if ( signal(SIGILL,  signal_handler) )
        write(STDERR_FILENO, "no  SIGILL handler\n", 19);
              
    if ( signal(SIGABRT, signal_handler) )
        write(STDERR_FILENO, "no SIGABRT handler\n", 19);
#endif

}
