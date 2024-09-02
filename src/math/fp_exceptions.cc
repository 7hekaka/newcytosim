
#include <stdlib.h>
#include <fenv.h>
#include <signal.h>
#include <stdio.h>


#if defined(__APPLE__)

/*
 Intel's OSX implementation adapted from the Corsika project
 https://gitlab.ikp.kit.edu/AirShowerPhysics/corsika
 by Lukas Nellen
*/

extern "C"
{
    static void fpe_signal_handler(int sig)
    {
        psignal(sig, "FPE Exception");
        exit(sig);
    }
    
#if defined(__arm64__)
    
    int enable_floating_point_exceptions()
    {
        signal(SIGFPE, fpe_signal_handler);
        fenv_t env;
        if (fegetenv(&env))
            return -1;
        env.__fpcr = env.__fpcr | __fpcr_trap_invalid;
        return fesetenv(&env);
    }
    
#elif defined(__x86_64__)
    
    int enable_floating_point_exceptions()
    {
        signal(SIGFPE, fpe_signal_handler);
        fenv_t env;
        if (fegetenv(&env))
            return -1;
        int old = env.__control & FE_ALL_EXCEPT;
        // unmask
        env.__control &= ~FE_ALL_EXCEPT;
        env.__mxcsr &= ~(FE_ALL_EXCEPT << 7);

        return fesetenv(&env) ? -1 : old;
    }

    int enable_floating_point_exceptions_except(int arg)
    {
        fenv_t env;
        int val = arg & FE_ALL_EXCEPT;
        int old;
        
        if (fegetenv(&env))
            return -1;
        old = env.__control & FE_ALL_EXCEPT;
        
        // unmask
        env.__control &= ~val;
        env.__mxcsr &= ~(val << 7);
        
        return fesetenv(&env) ? -1 : old;
    }
    
    
    int disable_floating_point_exceptions_except(int arg)
    {
        fenv_t env;
        int val = arg & FE_ALL_EXCEPT;
        int old;
        
        if (fegetenv(&env))
            return -1;
        old = env.__control & FE_ALL_EXCEPT;
        
        // mask
        env.__control |= val;
        env.__mxcsr |= val << 7;
        
        return fesetenv(&env) ? -1 : old;
    }
#endif
    
    #pragma STDC_FENV_ACCESS on
    void print_floating_point_exceptions(const char* str, FILE * out = stdout)
    {
        int n = fetestexcept(FE_ALL_EXCEPT);
        n &= ~FE_INEXACT;
        if ( n )
        {
            fprintf(out, " %s(", str);
            if (n&FE_DIVBYZERO)  fprintf(out, "DIVBYZERO ");
            if (n&FE_INEXACT)    fprintf(out, "INEXACT ");
            if (n&FE_INVALID)    fprintf(out, "INVALID ");
            if (n&FE_OVERFLOW)   fprintf(out, "OVERFLOW ");
            if (n&FE_UNDERFLOW)  fprintf(out, "UNDERFLOW ");
            if (n&FE_ALL_EXCEPT) fprintf(out, "UNKNOWN ");
            if ( feclearexcept(FE_ALL_EXCEPT) )
                fprintf(out, "unclear ");
            fprintf(out, "\b)");
        }
    }

}

#endif
