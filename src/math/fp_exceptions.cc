
#include <fenv.h>

#if defined(__APPLE__) && defined(__x86_64__)

/*
 OSX implementation adapted from the Corsika project
 https://gitlab.ikp.kit.edu/AirShowerPhysics/corsika
 by Lukas Nellen
*/

extern "C"
{
    int feenableexcept(int arg)
    {
        fenv_t fenv;
        int val = arg & FE_ALL_EXCEPT;
        int old;
        
        if (fegetenv(&fenv))
            return -1;
        old = fenv.__control & FE_ALL_EXCEPT;
        
        // unmask
        fenv.__control &= ~val;
        fenv.__mxcsr &= ~(val << 7);
        
        return fesetenv(&fenv) ? -1 : old;
    }
    
    
    int fedisableexcept(int arg)
    {
        fenv_t fenv;
        int val = arg & FE_ALL_EXCEPT;
        int old;
        
        if (fegetenv(&fenv))
            return -1;
        old = fenv.__control & FE_ALL_EXCEPT;
        
        // mask
        fenv.__control |= val;
        fenv.__mxcsr |= val << 7;
        
        return fesetenv(&fenv) ? -1 : old;
    }
    
}

#endif
