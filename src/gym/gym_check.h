// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "assert_macro.h"
#include <cstdio>


namespace gym
{
    
    /// convert OpenGL error code to string
    const char* errorString(unsigned code);
    
    /// check and print OpenGL error(s)
    void reportErrors(FILE*, const char* msg);
    
    /// print current color properties of OpenGL context
    void printColors(FILE*);
    
    /// print some info for debugging purpose
    void printCaps(const char[]);
    
    /// print OpenGL matrices for debugging purpose
    void printMatrices(FILE*);

};

#ifdef NDEBUG
#  define CHECK_GL_ERROR(ARG) ((void) 0)
#  define assert_enabled(ARG) ((void) 0)
#else
#  define CHECK_GL_ERROR(ARG) gym::reportErrors(stderr, ARG)
#  define assert_enabled(CAP) { if (!glIsEnabled(CAP))\
 { fprintf(stderr, "%s is not enabled in `%s`,  %s:%d\n", #CAP, SFUNC, SFILE, __LINE__); }}
#endif
